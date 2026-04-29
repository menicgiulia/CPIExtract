'''Gene identifier harmonization via MyGene.'''

import mygene
import pandas as pd
import time
from ..utils.helper import *

class MyGeneServer(metaclass=Singleton):
    '''
    Gene identifier harmonization via the MyGene API.

    Replacement for BiomartServer using the same search() and subset_search() interface.

    Supported input_type values:
        entrezgene_id, uniprotswissprot, ensembl_peptide_id,
        ensembl_gene_id, hgnc_id, chembl

    Note on ChEMBL target IDs
    -------------------------
    MyGene does not natively index ChEMBL target IDs (e.g. CHEMBL1977) as a queryable scope. 
    BioMart does have this natively, can flip to biomart if necessary.
    '''

    MAX_RETRIES  = 3
    BACKOFF_BASE = 2    # seconds: waits 2, 4, 8 between attempts
    CHUNK_SIZE   = 1000 # MyGene supports up to 1000 IDs per batch

    # Map BioMart filter names to MyGene query scopes
    SCOPE_MAP = {
        'entrezgene_id':      'entrezgene',
        'uniprotswissprot':   'uniprot.Swiss-Prot',
        'ensembl_peptide_id': 'ensembl.protein',
        'ensembl_gene_id':    'ensembl.gene',
        'hgnc_id':            'HGNC',
        'hgnc_symbol':        'symbol',
        'chembl':             'chembl',
    }

    # Map caller-defined column names to MyGene response fields.
    # Handles the naming variation across DB classes (e.g. 'entrez' vs 'entrezgene_id').
    COLUMN_FIELD_MAP = {
        'entrez':        'entrezgene',
        'entrezgene_id': 'entrezgene',
        'gene_type':     'type_of_gene',
        'gene_biotype':  'type_of_gene',
        'hgnc_symbol':   'symbol',
        'description':   'name',
    }

    def __init__(self) -> None:
        self.mg = mygene.MyGeneInfo()

    def _normalize_ids(self, input_type: str, input_ids: list) -> list:
        '''Strip known prefixes that MyGene.info doesn't expect.'''
        if input_type == 'hgnc_id':
            # BioMart passes 'HGNC:12679'; MyGene expects '12679'
            return [i.replace('HGNC:', '') if isinstance(i, str) else i for i in input_ids]
        return list(input_ids)

    def _hits_to_dataframe(self, hits: list, input_type: str,
                           original_ids: list, column_names: list) -> pd.DataFrame:
        '''
        Convert MyGene query many results to a DataFrame that matches column_names.
        The first column of column_names is always the input ID column — its values come from hit['query']. 
        All remaining columns are resolved via COLUMN_FIELD_MAP against MyGene response fields.
        '''
        records = []
        for hit in hits:
            if hit.get('notfound'):
                continue

            record = {}

            # First column: the input ID (query value as sent to MyGene)
            record[column_names[0]] = hit.get('query')

            # Remaining columns: map via COLUMN_FIELD_MAP
            for col in column_names[1:]:
                mygene_field = self.COLUMN_FIELD_MAP.get(col)
                if mygene_field is None:
                    record[col] = None
                    continue

                val = hit.get(mygene_field)

                # Normalize gene type: MyGene uses 'protein-coding', BioMart used 'protein_coding'
                if mygene_field == 'type_of_gene' and val == 'protein-coding':
                    val = 'protein_coding'
                record[col] = val

            records.append(record)

        if not records:
            return pd.DataFrame(columns=column_names)

        df = pd.DataFrame(records, columns=column_names)
        return df

    def search(self, input_type: str, input_id: list | str | int,
               attributes: list[str], column_names: list[str]) -> pd.DataFrame:
        '''
        Query MyGene with retry.
        The attributes argument is accepted for interface compatibility with BioMart.
        '''
        if isinstance(input_id, (str, int)):
            input_id = [input_id]
        input_id = list(input_id)

        scope = self.SCOPE_MAP.get(input_type)
        if scope is None:
            return pd.DataFrame(columns=column_names)

        normalized = self._normalize_ids(input_type, input_id)

        for attempt in range(1, self.MAX_RETRIES + 1):
            try:
                hits = self.mg.querymany(
                    normalized,
                    scopes=scope,
                    fields='entrezgene,type_of_gene,symbol,name',
                    species='human',
                    returnall=False,
                    verbose=False,
                )

                df = self._hits_to_dataframe(hits, input_type, normalized, column_names)

                # Retry if 0 rows returned for a non-empty input —
                # true biological negatives will consistently return 0 across all retries
                if len(df) == 0 and len(normalized) > 0:
                    raise ValueError(
                        f"MyGene returned 0 results for {len(normalized)} "
                        f"'{input_type}' IDs (attempt {attempt}/{self.MAX_RETRIES})."
                    )

                return df

            except Exception as e:
                if attempt < self.MAX_RETRIES:
                    wait = self.BACKOFF_BASE ** attempt
                    time.sleep(wait)

        return pd.DataFrame(columns=column_names)

    def subset_search(self, input_type: str, input_ids: list[int | str],
                      attributes: list[str], column_names: list[str]) -> pd.DataFrame:
        '''
        Chunk size is 1000 (MyGene batch limit)
        '''
        targets = pd.DataFrame(columns=column_names)

        if len(input_ids) > self.CHUNK_SIZE:
            for start, end in generate_subsets(len(input_ids), self.CHUNK_SIZE):
                sub = self.search(input_type, input_ids[start:end], attributes, column_names)
                targets = pd.concat([targets, sub], ignore_index=True)
        else:
            targets = self.search(input_type, input_ids, attributes, column_names)

        return targets.reset_index(drop=True)