'''Loading,searching,filtering and preprocessing data from PubChem.'''

import requests
import json
import time
import pandas as pd
import numpy as np
import pubchempy as pcp
import duckdb
import os

from ..utils.typing import Connection
from ..servers.PubChemServer import PubChemServer
from ..servers.BiomartServer import BiomartServer
from ..servers.MyGeneServer import MyGeneServer
from ..data_manager import *
from ..utils.helper import generate_subsets
from .Database import Database

class PubChem(Database):

    def __init__(self, connection: Connection| None=None, database: pd.DataFrame|None=None, 
                bioact_file=None, gene_server=None, server_select='mygene'):

        if gene_server is not None:
            self.gene_server = gene_server
        elif server_select == 'mygene':
            self.gene_server = MyGeneServer()
        elif server_select == 'biomart':
            self.gene_server = BiomartServer()
        else:
            raise ValueError(f"server_select must be 'mygene' or 'biomart', got '{server_select}'")

        if connection is not None or database is not None:
            raise ValueError('No SQL connection or database expected for Pubchem')
    
        self.db_file = None
        self.use_local = False
    
        if bioact_file:
            # Check if it's already a .duckdb file
            if bioact_file.endswith('.duckdb'):
                self.db_file = bioact_file
            else:
                self.db_file = None
        
            self.use_local = os.path.exists(self.db_file)
    
        funcs = {
            'comps' : self._retrieve_compounds,
            'targets' : self._retrieve_targets,
            'proteins' : self._retrieve_proteins
        }
        self.data_manager = APIManager(funcs)

    # Recognized concentration units in PubChem's Activity Unit field
    UNIT_TO_M = {'uM': 1e-6, 'nM': 1e-9}

    # activity names that are convertible to pChEMBL
    TRUSTED_ACTIVITY_NAMES = {'Ki', 'Kd', 'IC50', 'EC50', 'AC50', 'Potency', 'GI50'}

    def _filter_database(self, pubchem_raw: pd.DataFrame, identifier: str, pChEMBL_thres: float) -> pd.DataFrame:
        """
        Filters PubChem bioactivity data and derives a pChEMBL-equivalent value.

        Constraints
        -----------
            - Activity Outcome must not be 'Unspecified' or 'Inconclusive'. 'Inactive' is kept as negative evidence
              setting pChEMBL_eq = 0 only if there's no numeric information at all.
            - Activity Name != 'CD' (no potency measure)
            - Activity Unit is 'uM' or 'nM'
            - Activity Qualifier is used to find and control for relational operaters
              '>' -> weaker/upper bound, '<' -> stronger/lower bound. A '>'-censored
              row (or an 'Inactive') only counts as a confirmed negative if the derived value
              is <= pChEMBL_thres. A '<'-censored row only counts as a confirmed positive if it
              clears pChEMBL_thres.

        Parameters
        ----------
        pubchem_raw : DataFrame
            Raw PubChem bioactivity data
        identifier : str
            Column name identifying the compound or protein (e.g. 'CID', 'Target GeneID')
        pChEMBL_thres : float
            pChEMBL value used to classify a resolved interaction as positive (above) or
            negative (at or below)

        Returns
        -------
        DataFrame
            One row per activity record that passed filtering, with pChEMBL_eq / pChEMBL_lt /
            pChEMBL_gt and standard_type (Ki/IC50/Kd/EC50/Potency/...) columns.
        """

        # Handle different column names between API and bulk file
        activity_col = 'Activity Value' if 'Activity Value' in pubchem_raw.columns else 'Activity Value [uM]'

        # Drop invalid/empty identifiers
        pubchem_filt = pubchem_raw.dropna(subset=[identifier])
        pubchem_filt = pubchem_filt[~pubchem_filt[identifier].eq('')]
        # Remove records with ambiguous evidence
        pubchem_filt = pubchem_filt[(~pubchem_filt['Activity Outcome'].isin(['Unspecified', 'Inconclusive'])) &
                                    # Not a potency measure
                                    (pubchem_filt['Activity Name']!='CD')]

        # Find negative associations
        comment_negative_raw = pubchem_filt['Activity Outcome'] == 'Inactive'
        # Find associations without pChEMBL values with only eligible when Activity Name
        trusted_type_raw = pubchem_filt['Activity Name'].isin(self.TRUSTED_ACTIVITY_NAMES)
        bare_positive_raw = pubchem_filt['Activity Outcome'].isin(['Active', 'Probe']) & trusted_type_raw
        has_value_raw = pubchem_filt[activity_col].notna() & ~pubchem_filt[activity_col].astype(str).eq('')
        pubchem_filt = pubchem_filt[has_value_raw | comment_negative_raw | bare_positive_raw]

        # Only trust the value/qualifier for measurement types in TRUSTED_ACTIVITY_NAMES
        untrusted = ~pubchem_filt['Activity Name'].isin(self.TRUSTED_ACTIVITY_NAMES)
        pubchem_filt.loc[untrusted, activity_col] = np.nan

        # Remove duplicates with same identifier
        pubchem_filt = pubchem_filt.drop_duplicates(subset=[identifier, activity_col, 'Activity Name']).copy()

        # Get the relation activities
        if 'Activity Qualifier' in pubchem_filt.columns:
            relation = pubchem_filt['Activity Qualifier'].fillna('=').replace('', '=')
            numeric_str = pubchem_filt[activity_col].astype(str).str.replace(r'[^0-9.]', '', regex=True)
        else:
            raw_val = pubchem_filt[activity_col].astype(str)
            relation = raw_val.str.extract(r'^\s*(>=|<=|>|<)?')[0].fillna('=')
            numeric_str = raw_val.str.replace(r'[^0-9.]', '', regex=True)
        pubchem_filt[activity_col] = pd.to_numeric(numeric_str, errors='coerce')
        pubchem_filt['_relation'] = relation.values

        # Find activity unit
        if 'Activity Unit' in pubchem_filt.columns:
            no_value = pubchem_filt[activity_col].isnull()
            pubchem_filt = pubchem_filt[pubchem_filt['Activity Unit'].isin(self.UNIT_TO_M) | no_value]
            unit_to_M = pubchem_filt['Activity Unit'].map(self.UNIT_TO_M)
        else:
            # API data has no separate unit column
            unit_to_M = pd.Series(self.UNIT_TO_M['uM'], index=pubchem_filt.index)

        if len(pubchem_filt) == 0:
            return pubchem_filt.assign(
                pChEMBL_eq=pd.Series(dtype=float), pChEMBL_lt=pd.Series(dtype=float),
                pChEMBL_gt=pd.Series(dtype=float), standard_type=pd.Series(dtype=object))

        has_value = pubchem_filt[activity_col].notnull() & (pubchem_filt[activity_col] > 0)

        # Convert to M for the pChEMBL transform. NaN/non-positive  values propagate as NaN
        derived = -np.log10(pubchem_filt[activity_col].where(has_value) * unit_to_M)
        is_gt = pubchem_filt['_relation'].isin(['>', '>=', '>>']) & has_value
        is_lt = pubchem_filt['_relation'].isin(['<', '<=', '<<']) & has_value
        is_eq = (pubchem_filt['_relation'] == '=') & has_value

        pubchem_filt['pChEMBL_eq'] = derived.where(is_eq)
        pubchem_filt['pChEMBL_lt'] = derived.where(is_gt)
        pubchem_filt['pChEMBL_gt'] = derived.where(is_lt)
        pubchem_filt['standard_type'] = pubchem_filt['Activity Name']
        pubchem_filt = pubchem_filt.drop(columns=['_relation'])

        # Only a trusted potency measure is eligible to be treated as presence/absence evidence
        trusted_type = pubchem_filt['Activity Name'].isin(self.TRUSTED_ACTIVITY_NAMES)

        comment_negative = (pubchem_filt['Activity Outcome'] == 'Inactive') & trusted_type
        # A '>'-censored row is negative evidence if below pChEMBL_thres
        true_negatives = comment_negative | (pubchem_filt['pChEMBL_lt'].notnull() & (pubchem_filt['pChEMBL_lt'] <= pChEMBL_thres))
        # A '<'-censored row is positive evidence if above pChEMBL_thres
        censored_positive = pubchem_filt['pChEMBL_gt'].notnull() & (pubchem_filt['pChEMBL_gt'] > pChEMBL_thres)
        # An exact value is retained
        exact_value = pubchem_filt['pChEMBL_eq'].notnull()
        # A bare qualitative 'Active' call with no derivable value at all is kept as binary information
        bare_positive = (pubchem_filt['Activity Outcome'].isin(['Active', 'Probe']) &
            pubchem_filt['pChEMBL_eq'].isnull() &
            pubchem_filt['pChEMBL_lt'].isnull() &
            pubchem_filt['pChEMBL_gt'].isnull() &
            trusted_type)

        pubchem_filt = pubchem_filt.loc[exact_value | true_negatives | censored_positive | bare_positive].reset_index(drop=True)

        # Bare 'Inactive' calls with no derivable bound at all get an explicit 0. (Any surviving
        # 'Inactive' row here already passed the trusted_type check above via
        # comment_negative/true_negatives, so no need to re-check it here.)
        missing_val_negative = (
            (pubchem_filt['Activity Outcome'] == 'Inactive') &
            pubchem_filt['pChEMBL_eq'].isnull() &
            pubchem_filt['pChEMBL_lt'].isnull()
        )
        pubchem_filt.loc[missing_val_negative, 'pChEMBL_eq'] = 0

        return pubchem_filt

    def compounds(self, input_comp: pd.DataFrame, pChEMBL_thres: float=3.0, 
                    verbose: bool=False) -> tuple[pd.DataFrame, str, pd.DataFrame]:
        """
        Retrieves proteins from pubchem database interacting with compound passed as input.

        Steps
        -----
        - Filters database to obtain proteins interacting with input compound \\
        Constraints:
            - Only Homo Sapiens interactions
            - Activity Outcome must be specified
            - Activity unit is convertible to pchembl value.
        
        Parameters
        ----------
        input_comp : dictionary
            dictionary of input compound data from which interacting proteins are found
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid
        verbose : bool
            states whether API or Local file is in use
            
        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values: \\
            entrez, gene_type, hgnc_symbol, description, datasource (pc), pChEMBL_eq/lt/gt, standard_type
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all PubChem info about the input compound
        """

        # Print database info if verbose
        if verbose:
            if self.use_local:
                print(f"Using PubChem database: {self.db_file}")
            else:
                print("Using PubChem API")

        # On self for internal functions
        self._verbose = verbose

        columns = ['entrez','gene_type','hgnc_symbol','description','pChEMBL_eq','pChEMBL_lt','pChEMBL_gt','standard_type','datasource','inchikey']
        pubchem_act = pd.DataFrame(columns=columns)
        pubchem_raw = pd.DataFrame()
    
        # Determine CID list
        input_comp = input_comp.dropna(subset=['inchikey']).reset_index(drop=True)
        if len(input_comp) == 0:
            return pubchem_act, 'Input compound does not contain inchikey', pubchem_raw
        
        firstblock = input_comp['inchikey_fb'][0]
        
        # Try local database first
        if self.use_local:
            cid_list = self._get_cids_from_firstblock(firstblock)
            if not cid_list:  # Database empty, use API
                try:
                    cid_list = pcp.get_cids(firstblock, namespace='inchikey', searchtype=None)
                    if not cid_list:
                        return pubchem_act, 'No CIDs found for first block inchikey', pubchem_raw
                except Exception as e:
                    return pubchem_act, f'Error retrieving CIDs: {str(e)}', pubchem_raw
        else:  # No local database, use API
            try:
                cid_list = pcp.get_cids(firstblock, namespace='inchikey', searchtype=None)
                if not cid_list:
                    return pubchem_act, 'No CIDs found for first block inchikey', pubchem_raw
            except Exception as e:
                return pubchem_act, f'Error retrieving CIDs: {str(e)}', pubchem_raw
        
        # Retrieve bioactivities - local or API
        if self.use_local:
            # Use local files with DuckDB
            pubchem_act = self._retrieve_from_local(cid_list, pChEMBL_thres)
            
            if len(pubchem_act) == 0:
                return pubchem_act, 'No interaction data in local files', pubchem_raw
            
            statement = 'completed (local)'
            
        else: # Use API
            raw_cid_data = []
            for cid in cid_list:
                try:
                    raw_data = self.data_manager.retrieve_raw_data('comps', cid)
                    if len(raw_data) > 0:
                        raw_cid_data.append(raw_data)
                    time.sleep(0.5)
                except Exception as e:
                    continue
            
            if len(raw_cid_data) == 0:
                return pubchem_act, 'Failed to retrieve data from PubChem API', pubchem_raw
            
            pubchem_raw = pd.concat(raw_cid_data, ignore_index=True)
            
            if len(pubchem_raw) == 0:
                return pubchem_act, 'No interaction data', pubchem_raw
            
            # Filter database
            pubchem_act = self._filter_database(pubchem_raw, 'Target GeneID', pChEMBL_thres)
            
            # Get taxonomy - use it if SDQ already supplied it in _retrieve_compounds
            if 'TaxonomyID' not in pubchem_act.columns:
                pubchem_act['TaxonomyID'] = np.nan

            missing_genes = list(pubchem_act.loc[pubchem_act['TaxonomyID'].isnull(), 'Target GeneID'].unique())
            if len(missing_genes) > 0:
                gene_tax = self.data_manager.retrieve_raw_data('targets', missing_genes)
                gene_tax = gene_tax.rename(columns={'GeneID': 'Target GeneID', 'TaxonomyID': '_fetched_TaxonomyID'})
                pubchem_act = pubchem_act.merge(gene_tax, on='Target GeneID', how='left')
                pubchem_act['TaxonomyID'] = pubchem_act['TaxonomyID'].fillna(pubchem_act['_fetched_TaxonomyID'])
                pubchem_act = pubchem_act.drop(columns=['_fetched_TaxonomyID'])

            pubchem_act = pubchem_act.loc[pubchem_act['TaxonomyID']==9606].reset_index(drop=True)
            
            if len(pubchem_act) == 0:
                return pubchem_act, 'Filter reduced interactions to 0', pubchem_raw
            
            statement = 'completed (API)'
        
        # Get InChIKeys for CIDs
        unique_cids = pubchem_act['CID'].dropna().unique().tolist()
        if len(unique_cids) > 0:
            pc = PubChemServer()
            inchikey_columns = ['inchikey']
            selected_columns = pc.get_columns(inchikey_columns)

            try:
                cid_inchikey_map = pd.DataFrame()

                if len(unique_cids) > 1000:
                    for start, end in generate_subsets(len(unique_cids), 1000):
                        subset = [str(cid) for cid in unique_cids[start:end]]
                        batch_result = pc.get_compounds(subset, selected_columns, namespace='cid')
                        cid_inchikey_map = pd.concat([cid_inchikey_map, batch_result])
                        time.sleep(0.5)
                else:
                    cid_list_str = [str(cid) for cid in unique_cids]
                    cid_inchikey_map = pc.get_compounds(cid_list_str, selected_columns, namespace='cid')

                if len(cid_inchikey_map) > 0:
                    cid_inchikey_map = cid_inchikey_map.rename(columns=pc.properties)
                    pubchem_act['CID'] = pubchem_act['CID'].astype(str)
                    cid_inchikey_map['CID'] = cid_inchikey_map['CID'].astype(str)
                    pubchem_act = pubchem_act.merge(cid_inchikey_map[['CID', 'inchikey']], 
                                                    on='CID', how='left')
                else:
                    pubchem_act['inchikey'] = None
            except Exception as e:
                pubchem_act['inchikey'] = None

        # Unify gene identifiers by Harmonizing IDs
        ensembl = self.gene_server
        input_type = 'entrezgene_id' 
        attributes = ['entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
        names = ['entrez','gene_type','hgnc_symbol','description']

        genelist = list(pubchem_act['Target GeneID'])
        input_genes = [int(i) for i in genelist]

        pubchem_targets = ensembl.subset_search(input_type, input_genes, attributes, names)
        
        # For each compound, assign protein query column values to the ones from the original database
        pubchem_targets['entrez']=pubchem_targets['entrez'].astype(str)
        pubchem_act['Target GeneID']=pubchem_act['Target GeneID'].astype(str)
        for index, row in pubchem_act.iterrows():
            S1 = pubchem_targets.loc[pubchem_targets['entrez']==row['Target GeneID']]
            if len(S1) > 0:
                pubchem_act.loc[index,'entrez'] = S1['entrez'].iloc[0]
                pubchem_act.loc[index,'gene_type'] = S1['gene_type'].iloc[0]
                pubchem_act.loc[index,'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                pubchem_act.loc[index,'description'] = S1['description'].iloc[0]
                pubchem_act.loc[index, 'note'] ='Harmonized gene ID'
            else:
                pubchem_act.loc[index,'entrez'] = None
                pubchem_act.loc[index,'gene_type'] = None
                pubchem_act.loc[index,'hgnc_symbol'] = None
                pubchem_act.loc[index,'description'] = None
                pubchem_act.loc[index, 'note'] ='Failed to harmonize gene ID'
        pubchem_act['datasource'] = 'PubChem'

        return pubchem_act, statement, pubchem_raw

    def _retrieve_compounds(self, input_comp_id):
        # Uses PubChem's SDQ (Structured Data Query) agent instead of the PUG-REST to get activity qualifier.
        pubchem_raw = pd.DataFrame()
        if input_comp_id == '':
            return pubchem_raw

        sdq_url = 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi'
        query_template = {
            "select": ["*"],
            "collection": "bioactivity",
            "where": {"ands": [{"cid": str(input_comp_id)}]},
        }

        all_rows = []
        start = 1
        page_size = 10000
        while True:
            query = dict(query_template, start=start, limit=page_size)
            params = {"infmt": "json", "outfmt": "json", "query": json.dumps(query)}
 
            page_rows = None
            last_error = None
            for attempt in range(2):
                try:
                    response = requests.get(sdq_url, params=params, timeout=60)
                    time.sleep(0.5)
                    data = response.json()
                    page_rows = data['SDQOutputSet'][0]['rows']
                    break
                except Exception as e:
                    last_error = e
                    if attempt == 0:
                        time.sleep(1)
 
            if page_rows is None and getattr(self, '_verbose', False):
                # Both attempts failed. Stop paging rather than looping forever
                print(f"Warning: PubChem SDQ query failed for CID {input_comp_id} at page "
                      f"start={start} after 2 attempts ({last_error}). Results for this compound may be incomplete.")
            if page_rows is None:
                break
 
            all_rows.extend(page_rows)
            if len(page_rows) < page_size:
                break
            start += page_size

        if len(all_rows) == 0:
            return pubchem_raw

        raw = pd.DataFrame(all_rows)

        pubchem_raw = pd.DataFrame({
            'CID': raw.get('cid'),
            # SDQ returns geneid as a JSON float (e.g. 23097.0); normalize to the same clean integer-string
            'Target GeneID': raw.get('geneid').apply(lambda x: str(int(x)) if pd.notnull(x) else None) if 'geneid' in raw.columns else None,
            'Target Accession': raw.get('protacxn'),
            'Activity Outcome': raw.get('activity'),
            'Activity Name': raw.get('acname'),
            'Activity Qualifier': raw.get('acqualifier'),
            'Activity Value [uM]': raw.get('acvalue'),
            'Assay Name': raw.get('aidname'),
            'Assay Type': raw.get('aidtype'),
            'PubMed ID': raw.get('pmid'),
            'RNAi': raw.get('rnai'),
            'TaxonomyID': raw.get('taxid'),
        })

        return pubchem_raw
    
    def _retrieve_targets(self, gene_list) -> pd.DataFrame:
        # Use the PUG REST Pubchem API to assign the tax_id for each gene
        url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/gene/geneid/summary/JSON'

        def _fetch(subset, attempts=2):
            payload = {"geneid": ",".join(map(str, subset))}
            last_error = None
            for attempt in range(attempts):
                try:
                    response = requests.post(url, data=payload)
                    data = response.json()
                    gene_summary_list = data['GeneSummaries']['GeneSummary']
                    return pd.DataFrame(gene_summary_list, columns=['GeneID', 'TaxonomyID'])
                except Exception as e:
                    last_error = e
                    if attempt < attempts - 1:
                        time.sleep(1)
            # Both attempts failed, verbose print
            if getattr(self, '_verbose', False):
                print(f"Warning: PubChem gene taxonomy lookup failed for {len(subset)} gene(s) "
                      f"after {attempts} attempt(s) ({last_error}). These genes will be excluded.")
            return pd.DataFrame(columns=['GeneID', 'TaxonomyID'])

        if len(gene_list) > 1000:
            gene_tax = pd.DataFrame(columns=['GeneID', 'TaxonomyID'])
            for start, end in generate_subsets(len(gene_list), 1000):
                subset = gene_list[start:end]
                gene_tax = pd.concat([gene_tax, _fetch(subset)])
                time.sleep(0.5)
        else:
            gene_tax = _fetch(gene_list)
        gene_tax['GeneID'] = gene_tax['GeneID'].astype(str)
        return gene_tax
    
    def _retrieve_proteins(self, input_protein_id) -> pd.DataFrame:
        pubchem_raw = pd.DataFrame()
        # Use PUG REST to get the information of gene activity from PubChem
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/gene/geneid/%s/concise/JSON' %(str(input_protein_id))
 
        data = None
        last_error = None
        attempts = 2
        for attempt in range(attempts):
            try:
                response = requests.get(url)
                data = response.json()
                break
            except Exception as e:
                last_error = e
                if attempt < attempts - 1:
                    time.sleep(1)
 
        if data is None:
            # Both attempts failed
            if getattr(self, '_verbose', False):
                print(f"Warning: PubChem gene concise lookup failed for gene {input_protein_id} "
                      f"after {attempts} attempt(s) ({last_error}). Returning no data for this protein.")
            return pubchem_raw
        
        try:
            Des=data['Table']
        except:
            Des='No Data'

        if Des != 'No Data':
            cols=data['Table']['Columns']['Column']
            rows=[]
            for row in data['Table']['Row']:
                rows.append(row['Cell'])
            # Create dataframe from pubchem information
            pubchem_raw = pd.DataFrame(rows, columns=cols)

        return pubchem_raw

    def proteins(self, input_protein: pd.DataFrame, pChEMBL_thres: float=3.0, 
                verbose: bool=False):
        """
        Retrieves compounds from pubchem database interacting with proteins passed as input.

        Steps
        -----
        - Filters database to obtain proteins interacting with input compound \\
        Constraints:
            - Only Homo Sapiens interactions
            - Activity Outcome must be specified
            - Activity unit is convertible to pchembl value.
        
        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compound are found
        pChEMBL_thres : float
            pChEMBL value used to classify a resolved interaction as positive (above) or
            negative (at or below); censored bounds are only kept when they're tight enough to
            confidently fall on one side of this value.

        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values: \\
            inchi, inchikey, smiles, iupac_name, datasource (pc), pChEMBL_eq/lt/gt, standard_type
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all PubChem info about the protein
        """
        
        # Print database info if verbose
        if verbose:
            if self.use_local:
                print(f"Using PubChem database: {self.db_file}")
            else:
                print("Using PubChem API")

        # On self for internal functions
        self._verbose = verbose

        columns = ['inchi','inchikey','CID','smiles','connectivity_smiles','iupac_name','datasource','pChEMBL_eq','pChEMBL_lt','pChEMBL_gt','standard_type']
        pubchem_c1 = pd.DataFrame(columns=columns)
        pubchem_raw = pd.DataFrame()
        
        # Drop duplicated and null values of entrez
        input_protein = input_protein.dropna(subset=['entrez'])
        input_protein = input_protein.drop_duplicates(subset='entrez').reset_index(drop=True)
        
        if len(input_protein) == 0:
            pubchem_c1 = pubchem_c1[columns]
            return pubchem_c1, 'Input protein does not contain entrez', pubchem_raw
        
        # Get entrez ID
        input_protein_id = int(input_protein['entrez'][0])
        
        # Try local database first
        if self.use_local:
            try:
                con = duckdb.connect(self.db_file, read_only=True)
                con.execute("SET enable_progress_bar=false")
                
                # Query compounds that interact with this protein
                query = f"""
                    SELECT DISTINCT c.CID, c.InChI, c.InChIKey, 
                        b."Activity Value", b."Activity Name", b."Activity Outcome", b."Activity Unit", b."Activity Qualifier"
                    FROM bioactivities b
                    INNER JOIN cid_inchikey c ON b.CID = c.CID
                    WHERE TRY_CAST(b."Gene ID" AS INTEGER) = {input_protein_id}
                    AND b."Activity Outcome" NOT IN ('Unspecified', 'Inconclusive')
                    AND b."Activity Name" != 'CD'
                    AND (b."Activity Unit" IN ('uM', 'nM') OR b."Activity Value" IS NULL)
                """
                result = con.execute(query).df()
                con.close()
                
                if len(result) > 0:
                    # Rename columns
                    result = result.rename(columns={
                        'InChI': 'inchi',
                        'InChIKey': 'inchikey'
                    })
                    
                    # Apply the same relation-aware filtering/derivation as every other path
                    result = self._filter_database(result, 'CID', pChEMBL_thres)
                    
                    if len(result) > 0:
                        # Create output with all required columns
                        pubchem_c1 = result[['inchi', 'inchikey', 'CID', 'pChEMBL_eq', 'pChEMBL_lt', 'pChEMBL_gt', 'standard_type']].drop_duplicates()
                        pubchem_c1['smiles'] = None  # Will be filled in postprocessing
                        pubchem_c1['iupac_name'] = None  # Will be filled in postprocessing
                        pubchem_c1['datasource'] = 'PubChem'
                        
                        # Reorder columns
                        pubchem_c1 = pubchem_c1[columns]
                        
                        statement = 'completed (local)'
                    else:
                        statement = 'Filter reduced interactions to 0'
                else:
                    statement = 'No interaction data in local database'
                    
            except Exception as e:
                if verbose:
                    print(f"Local database query failed: {e}, falling back to API")
                self.use_local = False  # Fall through to API
        
        # Use API if local failed or not available
        if not self.use_local or len(pubchem_c1) == 0:
            pubchem_raw = self.data_manager.retrieve_raw_data('proteins', input_protein_id)

            if len(pubchem_raw) > 0:
                # Filter for high quality interactions
                pubchem_act = self._filter_database(pubchem_raw, 'CID', pChEMBL_thres)

                if len(pubchem_act) > 0:
                    # Get unique CIDs
                    unique_cids = pubchem_act['CID'].dropna().unique().tolist()
                    
                    if len(unique_cids) > 0:
                        pc = PubChemServer()
                        try:
                            cid_list = [str(cid) for cid in unique_cids]
                            selected_columns = pc.get_columns(['inchi', 'inchikey', 'smiles', 'iupac_name'])
                            
                            if len(cid_list) > 1000:
                                compounds = pd.DataFrame()
                                for start, end in generate_subsets(len(cid_list), 1000):
                                    subset = cid_list[start:end]
                                    batch = pc.get_compounds(subset, selected_columns, namespace='cid')
                                    compounds = pd.concat([compounds, batch])
                                    time.sleep(0.5)
                            else:
                                compounds = pc.get_compounds(cid_list, selected_columns, namespace='cid')
                            
                            compounds = compounds.rename(columns=pc.properties)
                            
                            if len(compounds) > 0:
                                pubchem_act['CID'] = pubchem_act['CID'].astype(int)
                                compounds['CID'] = compounds['CID'].astype(int)
                                
                                pubchem_info = pubchem_act[['CID', 'pChEMBL_eq', 'pChEMBL_lt', 'pChEMBL_gt', 'standard_type']].drop_duplicates()
                                
                                pubchem_c1 = pd.merge(compounds, pubchem_info, on='CID', how='left')
                                pubchem_c1['datasource'] = 'PubChem'
                                
                                # Ensure all columns exist
                                for col in columns:
                                    if col not in pubchem_c1.columns:
                                        pubchem_c1[col] = None
                                
                                pubchem_c1 = pubchem_c1[columns]
                                statement = 'completed (API)'
                            else:
                                statement = 'Compounds not found using PubChem'
                        except Exception as e:
                            statement = f'Error retrieving compound info: {str(e)}'
                    else:
                        statement = 'No valid CIDs found'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        
        pubchem_c1 = pubchem_c1[columns]
        return pubchem_c1, statement, pubchem_raw

    def _get_cids_from_firstblock(self, firstblock):
        """Get all CIDs matching a firstblock from database"""
    
        if not self.use_local:
            return None
        con = None
    
        try:
            con = duckdb.connect(self.db_file, read_only=True)
            con.execute("SET enable_progress_bar=false")

            result = con.execute(f"""
                SELECT DISTINCT CID 
                FROM cid_inchikey 
                WHERE firstblock = '{firstblock}'
            """).df()
            con.close()
        
            if len(result) == 0:
                return []
            return result['CID'].tolist()
        
        #except Exception as e:
        #    return None
        finally:
            if con is not None:
                con.close()

    def _retrieve_from_local(self, cid_list, pChEMBL_thres):
        """Retrieve bioactivities from database with SQL JOIN for human genes"""

        if not cid_list: # do not check duckdb if cid_list is empty
            return pd.DataFrame()

        con=None
        try:
            con = duckdb.connect(self.db_file, read_only=True)
            con.execute("SET enable_progress_bar=false")

            cid_str = ','.join(map(str, cid_list))
    
            pubchem_raw = con.execute(f"""
                SELECT b.* EXCLUDE ("Target TaxID"), g.TaxonomyID
                FROM bioactivities b
                INNER JOIN gene_info g ON TRY_CAST(b."Gene ID" AS INTEGER) = g.GeneID
                WHERE b.CID IN ({cid_str})
                AND g.TaxonomyID = 9606
                AND b."Activity Outcome" NOT IN ('Unspecified', 'Inconclusive')
                AND b."Activity Name" != 'CD'
                AND b."Gene ID" IS NOT NULL
                AND (b."Activity Unit" IN ('uM', 'nM') OR b."Activity Value" IS NULL)
            """).df()

            if len(pubchem_raw) == 0:
                return pd.DataFrame()
    
            # Rename columns to match API format
            pubchem_raw['Gene ID'] = pubchem_raw['Gene ID'].astype(int)
            pubchem_raw = pubchem_raw.rename(columns={
                'Gene ID': 'Target GeneID',
                'Protein Accession': 'Target Accession'
            })
        
            # Apply pChEMBL filtering
            pubchem_filt = self._filter_database(pubchem_raw, 'Target GeneID', pChEMBL_thres)
            return pubchem_filt
        
        #except Exception as e:
        #    return pd.DataFrame()
        finally:
            if con is not None:
                con.close()