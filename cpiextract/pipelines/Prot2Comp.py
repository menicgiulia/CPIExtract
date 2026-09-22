'''Retrieve compound interactions for proteins provided as input.'''

import numpy as np
import pandas as pd
from ..databases import *
from ..utils.identifiers import protein_identifiers
from ..servers.PubChemServer import PubChemServer
from .Pipeline import Pipeline

class Prot2Comp(Pipeline):
    '''Retrieve compound interactions for proteins provided as input.'''

    def _update_args(self, pChEMBL_thres: float, dtc_mutated: bool, dc_extra: bool,
                     verbose: bool, experimental_thres: float) -> None:
        
        # Keyword-argument dicts, each database's compounds() call binds arguments by name
        self.database_args = {
            'pc':     {'pChEMBL_thres': pChEMBL_thres, 'verbose': verbose},
            'chembl': {'pChEMBL_thres': pChEMBL_thres},
            'bdb':    {'pChEMBL_thres': pChEMBL_thres},
            'stitch': {'experimental_thres': experimental_thres},
            'ctd':    {},
            'dtc':    {'dtc_mutated': dtc_mutated, 'pChEMBL_thres': pChEMBL_thres, 'verbose': verbose},
            'otp':    {},
            'dc':     {'dc_extra': dc_extra, 'pChEMBL_thres': pChEMBL_thres},
            'db':     {},
        }      

    # Calls functions to collect data and merges all the data from the various sources together
    # The parameters are:
    #    - input_id - the protein id
    #    - pChEMBL_thresh - the minimum interaction pChEMBL value required to be added to the output file
    #    - dtc_mutated - to select whether also to consider interactions with mutated target proteins from DTC
    #    - dc_extra - to select whether to include possibly non-Homo sapiens interactions
    #    - pchembl_grouping - how to compute the average pChEMBL for a pair: 'all' (combined),
    #      'type_group' (K-types vs C50-types separately), or 'unique' (one average per exact type)
    #    - experimental_thres - minimum STITCH/STRING 'experimental' confidence score (0-999) required
    #    - strong_positive_thres - pchembl_eq/pchembl_gt average above this is "strong positive",
    #      above pChEMBL_thres but at or below this is "weak positive"

    def prot_interactions(self, input_id: int|str, pChEMBL_thres: float=3.0, dtc_mutated: bool=False, dc_extra: bool=False, 
                          verbose: bool=False, pchembl_grouping: str='all',
                          experimental_thres: float=400, strong_positive_thres: float=6.0) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:

        # Run interaction select with all databases selected
        prot_comp, states, raw_comp = self.prot_interactions_select(input_id, pChEMBL_thres=pChEMBL_thres, dtc_mutated=dtc_mutated, 
                                                        dc_extra=dc_extra,verbose=verbose,
                                                        pchembl_grouping=pchembl_grouping,
                                                        experimental_thres=experimental_thres,
                                                        strong_positive_thres=strong_positive_thres)
        
        return prot_comp, states, raw_comp

    def prot_interactions_select(self, input_id: int|str, selected_dbs: str='pc_chembl_bdb_stitch_ctd_dtc_otp_dc_db', 
                                 pChEMBL_thres: float=3.0, dtc_mutated: bool=False, dc_extra: bool=False,
                                 verbose: bool=False, pchembl_grouping: str='all', experimental_thres: float=400,
                                 strong_positive_thres: float=6.0) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:

        self._update_args(pChEMBL_thres, dtc_mutated, dc_extra, verbose, experimental_thres)

        prot_ids = protein_identifiers(input_id,gene_server=self.gene_server)

        unified_columns = ['inchikey', 'CID', 'smiles', 'connectivity_smiles', 'iupac_name', 'synonyms',
                            'entrez', 'hgnc_symbol', 'gene_type', 'description',
                            'pchembl_eq', 'pchembl_lt', 'pchembl_gt', 'standard_type', 'datasource']

        prot_comp = pd.DataFrame(columns=unified_columns)
        raw_comp = pd.DataFrame()

        # Outputs the brief statement regarding the processing of each source 
        states=pd.DataFrame(columns=['prot'])
        states.loc[0, 'prot'] = input_id   
        
        if len(prot_ids) > 0:
            comp_all_list = []
            comp_all_raw_list = []
            db_list = selected_dbs.split('_')
            for name, db in self.databases.items():
                if name in db_list:
                    # Perform database search, protein-input direction
                    result, state, _ = db.proteins(prot_ids, **self.database_args.get(name, {}))
                    if len(result) > 0:
                        # Keep each database's own native, filtered output as-is
                        comp_all_raw_list.append(result.copy())

                        # Harmonize this database's native column names for aggregation
                        result = self._harmonize_columns(result)
                        # Append to list
                        comp_all_list.append(result)
                    if verbose:
                        print(f'{name} done!')
                else:
                    state = f'Did not run {name}'
                # Add statement for this source
                states.loc[0, name] = state
        
            if len(comp_all_raw_list) > 0:
                raw_comp = pd.concat(comp_all_raw_list, ignore_index=True)

            # Concatenate once AFTER loop
            if len(comp_all_list) > 0:
                comp_all = pd.concat(comp_all_list, ignore_index=True)
            else:
                comp_all = pd.DataFrame(columns=unified_columns)

            prot_comp = self._postprocess_databases(comp_all, pChEMBL_thres, pchembl_grouping, strong_positive_thres)
                        
            if len(prot_comp) > 0:
                # Add std compound ids to output
                std_ids = prot_ids[['entrez', 'gene_type', 'hgnc_symbol', 'description']]
                # Select first nonNa value from each column
                std_ids = std_ids.apply(lambda col: std_ids[col.name].dropna().iloc[0] if not std_ids[col.name].dropna().empty else None)
                # Assign to each row the standard ids
                prot_comp = prot_comp.assign(**std_ids)
                # Add input id to the results
                prot_comp['input_id'] = input_id

                # Reorder columns
                fixed_columns = ['input_id', 'entrez', 'gene_type', 'hgnc_symbol', 'description',
                                    'inchi', 'inchikey', 'inchikey_fb', 'CID', 'smiles','connectivity_smiles', 'iupac_name',
                                    'synonyms']
                pchembl_columns = [c for c in prot_comp.columns if c.startswith(('pchembl_count', 'ave_pchembl', 'std_pchembl'))]
                trailing_columns = ['interaction_class', 'src_count', 'pubchem', 'chembl', 'bindingdb', 'stitch',
                                    'ctd', 'dtc', 'otp', 'drugcentral', 'drugbank']
                prot_comp = prot_comp[fixed_columns + pchembl_columns + trailing_columns]

        return prot_comp, states, raw_comp
            

    def _postprocess_databases(self, comp_all: pd.DataFrame, pChEMBL_thres: float, pchembl_grouping: str,
                               strong_positive_thres: float = 6.0) -> pd.DataFrame:
        # Remove compounds without InChI
        comp_all = comp_all.dropna(subset=['inchikey']).reset_index(drop=True)
    
        # Add first-block column for reference
        pc = PubChemServer()
        comp_all['inchikey_fb'] = comp_all['inchikey'].apply(
            lambda x: pc.get_inchikey_first_block(x) if pd.notna(x) else None
        )
    
        # Deduplicate by full InChIKey (not by first-block)
        comp_list = comp_all['inchikey'].unique()

        tar_comp = pd.DataFrame(columns=['inchi', 'inchikey', 'inchikey_fb', 'CID', 'smiles', 'connectivity_smiles','iupac_name',
                                     'synonyms', 'src_count'] + 
                                     [source.lower() for source in self.sources])
    
        for index, inchikey in enumerate(comp_list):
            comp = comp_all.loc[comp_all['inchikey'] == inchikey]
        
            if len(comp) > 0:
                # Take first occurrence for compound identifiers
                tar_comp.loc[index, 'inchi'] = comp['inchi'].iloc[0] if 'inchi' in comp.columns else None
                tar_comp.loc[index, 'inchikey'] = comp['inchikey'].iloc[0]
                tar_comp.loc[index, 'inchikey_fb'] = comp['inchikey_fb'].iloc[0]
                tar_comp.loc[index, 'CID'] = comp['CID'].iloc[0] if 'CID' in comp.columns else None
                tar_comp.loc[index, 'smiles'] = comp['smiles'].iloc[0] if 'smiles' in comp.columns else None
                tar_comp.loc[index, 'connectivity_smiles'] = comp['connectivity_smiles'].iloc[0] if 'connectivity_smiles' in comp.columns else None
                tar_comp.loc[index, 'iupac_name'] = comp['iupac_name'].iloc[0] if 'iupac_name' in comp.columns else None
                tar_comp.loc[index, 'synonyms'] = self._top_synonyms(comp['synonyms']) if 'synonyms' in comp.columns else None

                # Aggregate pchembl values
                tar_comp = self._aggregate_pchembl(tar_comp, index, comp, pChEMBL_thres, pchembl_grouping, strong_positive_thres)
            
                # Count sources
                tar_comp.loc[index, 'src_count'] = len(comp['datasource'].unique())
            
                # Create source binary matrix
                for source in self.sources:
                    tar_comp.loc[index, source.lower()] = 1 if source in comp['datasource'].unique() else 0
    
        return tar_comp