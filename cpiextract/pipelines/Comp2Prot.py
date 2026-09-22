'''Retrieve proteins interacting with the small molecule passed as input.'''

import numpy as np
import pandas as pd
from ..utils.identifiers import compound_identifiers
from ..databases import *
from .Pipeline import Pipeline

class Comp2Prot(Pipeline):
    '''Retrieve proteins interacting with the small molecule passed as input.'''

    def _update_args(self, pChEMBL_thres: float, dtc_mutated: bool, dc_extra: bool,
                     verbose: bool, experimental_thres: float,protein_types: set|None):

        # Keyword-argument dicts, not positional tuples
        self.database_args = {
            'pc':     {'pChEMBL_thres': pChEMBL_thres, 'verbose': verbose},
            'chembl': {'pChEMBL_thres': pChEMBL_thres},
            'bdb':    {'pChEMBL_thres': pChEMBL_thres},
            'stitch': {'experimental_thres': experimental_thres},
            'ctd':    {},
            'dtc':    {'dtc_mutated': dtc_mutated, 'pChEMBL_thres': pChEMBL_thres, 'verbose': verbose},
            'otp':    {},
            'dc':     {'dc_extra': dc_extra, 'pChEMBL_thres': pChEMBL_thres},
            'db':     {'protein_types': protein_types},
        }       

    # Calls functions to collect data and merges all the data from the various sources together  
    # The parameters are:
    #    - input_id - the compound id
    #    - pChEMBL_thresh - the minimum interaction pChEMBL value required to be added to the output file
    #    - dtc_mutated - to select whether also to consider interactions with mutated target proteins from DTC
    #    - dc_extra - to select whether to include possibly non-Homo sapiens interactions
    #    - pchembl_grouping - how to compute the average pChEMBL for a pair: 'all' (combined),
    #      'type_group' (K-types vs C50-types separately), or 'unique' (one average per exact type)
    #    - experimental_thres - minimum STITCH/STRING 'experimental' confidence score (0-999) required
    #    - protein_types - which DrugBank protein_type categories to include (target/enzyme/carrier/
    #      transporter); defaults to all four if not specified
    #    - strong_positive_thres - pchembl_eq/pchembl_gt average above this is "strong positive",
    #      above pChEMBL_thres but at or below this is "weak positive"

    def comp_interactions(self, input_id: int|str, pChEMBL_thres: float=3.0, 
                    dtc_mutated: bool=False, dc_extra: bool=False,
                    verbose: bool=False, pchembl_grouping: str='all', experimental_thres: float=400,
                    protein_types: set|None=None, strong_positive_thres: float=6.0,
                    prebuilt_comp_ids: pd.DataFrame|None=None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    
        # Run interaction select with all databases selected
        comp_tar, states, raw_tar = self.comp_interactions_select(input_id, pChEMBL_thres=pChEMBL_thres, 
                                                dtc_mutated=dtc_mutated, dc_extra=dc_extra, 
                                                verbose=verbose,
                                                pchembl_grouping=pchembl_grouping,
                                                experimental_thres=experimental_thres,
                                                protein_types=protein_types,
                                                strong_positive_thres=strong_positive_thres,
                                                prebuilt_comp_ids=prebuilt_comp_ids)
        return comp_tar, states, raw_tar


    def comp_interactions_select(self, input_id: int|str, selected_dbs: str='pc_chembl_bdb_stitch_ctd_dtc_otp_dc_db', 
                             pChEMBL_thres: float=3.0, dtc_mutated: bool=False, dc_extra: bool=False, 
                             verbose: bool=False, pchembl_grouping: str='all',
                             experimental_thres: float=400, protein_types: set|None=None,
                             strong_positive_thres: float=6.0,
                             prebuilt_comp_ids: pd.DataFrame|None=None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:

        self._update_args(pChEMBL_thres, dtc_mutated, dc_extra, verbose,
                          experimental_thres, protein_types)

        unified_columns = ['inchikey', 'CID', 'smiles', 'connectivity_smiles', 'iupac_name', 'synonyms',
                            'entrez', 'hgnc_symbol', 'gene_type', 'description',
                            'pchembl_eq', 'pchembl_lt', 'pchembl_gt', 'standard_type', 'datasource']
        comp_tar = pd.DataFrame(columns=unified_columns)
        raw_tar = pd.DataFrame()

        # Use prebuilt comp_ids if provided, otherwise call compound_identifiers
        if prebuilt_comp_ids is not None:
            comp_ids = prebuilt_comp_ids
        else:
            comp_ids = compound_identifiers(input_id)

        states = pd.DataFrame(columns=['comp'])
        states.loc[0, 'comp'] = input_id

        if len(comp_ids) > 0:
            tar_all_list = []
            tar_all_raw_list = []
            db_list = selected_dbs.split('_')
            for name, db in self.databases.items():
                if name in db_list:
                    # Perform database search
                    result, state, _ = db.compounds(comp_ids, **self.database_args.get(name, {}))
                    if len(result) > 0:
                        # Keep each database's own native, filtered output as-is
                        tar_all_raw_list.append(result.copy())

                        # Harmonize this database's native column names for aggregated
                        result = self._harmonize_columns(result)
                        # Append to list
                        tar_all_list.append(result)
                    if verbose:
                        print(f'{name} done!')
                else:
                    state = f'Did not run {name}'
                # Add statement for this source
                states.loc[0, name] = state
        
            if len(tar_all_raw_list) > 0:
                raw_tar = pd.concat(tar_all_raw_list, ignore_index=True)

            if len(tar_all_list) > 0: # Concatenate
                # Filter out empty DataFrames before concatenating
                tar_all_list = [df for df in tar_all_list if len(df) > 0]
                if len(tar_all_list) > 0:
                    tar_all = pd.concat(tar_all_list, ignore_index=True)
                else:
                    tar_all = pd.DataFrame(columns=unified_columns)
            else:
                tar_all = pd.DataFrame(columns=unified_columns)
            
            tar_all=tar_all.rename(columns={'inchikey':'db_inchikey'})
            comp_tar = self._postprocess_databases(tar_all, pChEMBL_thres, pchembl_grouping, strong_positive_thres)

            if len(comp_tar) > 0:
                # Add std compound ids to output
                comp_ids=comp_ids.rename(columns={'smiles':'pc_iso_smiles','connectivity_smiles':'pc_canonical_smiles','inchi':'pc_inchi',
                                                  'inchikey':'pc_inchikey','inchikey_fb':'pc_firstblock',
                                                  'iupac_name':'pc_iupac_name','CID':'pc_cid'})
                std_ids = comp_ids[['pc_iso_smiles', 'pc_canonical_smiles', 'pc_inchi', 'pc_inchikey', 'pc_firstblock', 'pc_iupac_name', 'pc_cid']]
                # Select first nonNa value from each column
                std_ids = std_ids.apply(lambda col: std_ids[col.name].dropna().iloc[0] if not std_ids[col.name].dropna().empty else None)
                # Assign to each row the standard ids
                comp_tar = comp_tar.assign(**std_ids)       
                # Add input id to the results
                comp_tar['input_id'] = input_id

                if 'synonyms' in comp_ids.columns:
                    comp_tar['synonyms'] = self._top_synonyms(comp_ids['synonyms'])

                def _structure_match(row):
                    db_ik, pc_ik = row.get('db_inchikey'), row.get('pc_inchikey')
                    if pd.isna(db_ik) or pd.isna(pc_ik):
                        return None
                    db_blocks, pc_blocks = str(db_ik).split('-'), str(pc_ik).split('-')
                    if len(db_blocks) < 2 or len(pc_blocks) < 2:
                        return None
                    if db_blocks[0] != pc_blocks[0]:
                        return None
                    return 'stereochemical' if db_blocks[1] == pc_blocks[1] else 'scaffold'

                comp_tar['structure_match'] = comp_tar.apply(_structure_match, axis=1)

                # Reorder columns
                fixed_columns = ['input_id', 'pc_inchi', 'pc_inchikey', 'pc_firstblock', 'pc_iso_smiles', 'pc_canonical_smiles', 
                                     'pc_iupac_name','pc_cid', 'db_inchikey', 'structure_match', 'entrez', 'gene_type', 
                                     'hgnc_symbol', 'description', 'synonyms']
                pchembl_columns = [c for c in comp_tar.columns if c.startswith(('pchembl_count', 'ave_pchembl', 'std_pchembl'))]
                trailing_columns = ['interaction_class', 'src_count', 'pubchem',
                                    'chembl', 'bindingdb', 'stitch', 'ctd', 'dtc', 'otp', 'drugcentral',
                                    'drugbank']
                comp_tar = comp_tar[fixed_columns + pchembl_columns + trailing_columns]

        return comp_tar, states, raw_tar
    

    def _postprocess_databases(self, tar_all: pd.DataFrame, pChEMBL_thres: float, pchembl_grouping: str,
                               strong_positive_thres: float = 6.0) -> pd.DataFrame:
        # Remove non-protein coding interactions
        tar_all = tar_all[tar_all['gene_type']=='protein_coding']
        # Remove proteins with no Symbol or Entrez
        tar_all = tar_all.dropna(subset=['hgnc_symbol', 'entrez']).reset_index(drop=True)  
    
        # Group by BOTH protein AND compound (to keep stereoisomers separate)
        tar_list = tar_all[['hgnc_symbol', 'db_inchikey']].drop_duplicates().values.tolist()
    
        # Create output dataframe
        comp_tar = pd.DataFrame(columns=['entrez','hgnc_symbol','description','gene_type','db_inchikey',
                                   'synonyms', 'src_count'] + 
                                   [source.lower() for source in self.sources])
    
        for index, (target, inchikey) in enumerate(tar_list):
            # Select rows for this protein-compound pair
            tar = tar_all.loc[(tar_all['hgnc_symbol'] == target) & (tar_all['db_inchikey'] == inchikey)]
        
            # Copy the protein ID data to the output dataframe
            comp_tar.loc[index,'entrez'] = int(tar['entrez'].iloc[0])
            comp_tar.loc[index,'hgnc_symbol'] = tar['hgnc_symbol'].iloc[0]
            comp_tar.loc[index,'description'] = tar['description'].iloc[0]
            comp_tar.loc[index,'gene_type'] = tar['gene_type'].iloc[0]
            comp_tar.loc[index,'db_inchikey'] = inchikey  # Add the compound inchikey
            comp_tar.loc[index,'synonyms'] = self._top_synonyms(tar['synonyms'])

            comp_tar = self._aggregate_pchembl(comp_tar, index, tar, pChEMBL_thres, pchembl_grouping, strong_positive_thres)
        
            # Count sources for this protein-compound pair
            comp_tar.loc[index, 'src_count'] = len(tar['datasource'].unique())
        
            # Create source binary matrix
            for source in self.sources:
                comp_tar.loc[index, source.lower()] = 1 if source in tar['datasource'].unique() else 0

        return comp_tar