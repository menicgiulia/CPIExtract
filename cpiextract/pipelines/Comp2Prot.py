'''Retrieve proteins interacting with the small molecule passed as input.'''

import numpy as np
import pandas as pd
from ..utils.identifiers import compound_identifiers
from ..databases import *
from .Pipeline import Pipeline

class Comp2Prot(Pipeline):
    '''Retrieve proteins interacting with the small molecule passed as input.'''

    def _update_args(self, pChEMBL_thres: float, dtc_mutated: bool, dc_extra: bool, chembl_ids: list[str], merge_stereoisomers: bool, verbose: bool):
        self.database_args = {
            'pc': (pChEMBL_thres,merge_stereoisomers,),
            'chembl': (chembl_ids, pChEMBL_thres,merge_stereoisomers,),
            'bdb': (pChEMBL_thres,merge_stereoisomers,),
            'stitch': (merge_stereoisomers,),
            'ctd': (merge_stereoisomers,),
            'dtc': (chembl_ids, dtc_mutated, pChEMBL_thres,merge_stereoisomers,),
            'otp': (chembl_ids,merge_stereoisomers,),
            'dc': (dc_extra, pChEMBL_thres,merge_stereoisomers,),
            'db': (merge_stereoisomers,),
        }       


    # Calls functions to collect data and merges all the data from the various sources together  
    # The parameters are:
    #    - input_id - the compound id
    #    - pChEMBL_thresh - the minimum interaction pChEMBL value required to be added to the output file
    #    - dtc_mutated - to select whether also to consider interactions with mutated target proteins from DTC
    #    - dc_extra - to select whether to include possibly non-Homo sapiens interactions
    #    - merge_stereoisomers - to select whether data collected is stereo-specific or generic to a structure

    def comp_interactions(self, input_id: int|str, pChEMBL_thres: float=0, 
                        dtc_mutated: bool=False, dc_extra: bool=False, merge_stereoisomers: bool=False, verbose: bool=False) -> tuple[pd.DataFrame, pd.DataFrame]:
        
        # Run interaction select with all databases selected
        comp_tar, state = self.comp_interactions_select(input_id, pChEMBL_thres=pChEMBL_thres, 
                                                    dtc_mutated=dtc_mutated, dc_extra=dc_extra, merge_stereoisomers=merge_stereoisomers, verbose=verbose)
        return comp_tar, state

    # Calls functions to collect data and merges all the data from the selected sources together
    # The parameters are:
    #    - input_id - the compound id
    #    - selected_dbs - underscore-separated string containing the names of the databases which are selected
    #    - pChEMBL_thresh - the minimum interaction pChEMBL value required to be added to the output file
    #    - dtc_mutated - to select whether also to consider interactions with mutated target proteins from DTC
    #    - dc_extra - to select whether to include possibly non-Homo sapiens interactions
    #    - merge_stereoisomers - to select whether data collected is stereo-specific or generic to a structure

    def comp_interactions_select(self, input_id: int|str, selected_dbs: str='pc_chembl_bdb_stitch_ctd_dtc_otp_dc_db', 
                                 pChEMBL_thres: float=0, dtc_mutated: bool=False, dc_extra: bool=False, merge_stereoisomers: bool=False, verbose: bool=False) -> tuple[pd.DataFrame, pd.DataFrame]:
        #pChEMBL_thres default 0: Will only interactions without activity data in sources with activity data present
        #merge_stereoisomers default False: False=selects interactions for a structure without stereo-specificity

        chembl_ids: list[str] = []
        self._update_args(pChEMBL_thres, dtc_mutated, dc_extra, chembl_ids, merge_stereoisomers,verbose)

        main_columns = ['entrez','gene_type','hgnc_symbol','description','pchembl_value','datasource','inchikey']
        comp_tar = pd.DataFrame(columns=main_columns)
        #Collect all the data through the other functions
        comp_ids = compound_identifiers(input_id) #find comp_ids
        #Outputs the brief statement regarding the processing of each source 
        states = pd.DataFrame(columns=['comp'])
        states.loc[0, 'comp'] = input_id

        if len(comp_ids) > 0:
            tar_all_list = []
            db_list = selected_dbs.split('_')
            for name, db in self.databases.items():
                if name in db_list:
                    # Perform database search
                    result, state, _ = db.interactions(comp_ids, *self.database_args.get(name, ()))
                    if len(result) > 0:
                        # Keep only wanted columns
                        result = result.loc[:, main_columns]
                        # Append to list
                        tar_all_list.append(result)
                    if verbose:
                        print(f'{name} done!')
                else:
                    state = f'Did not run {name}'
                # Add statement for this source
                states.loc[0, name] = state
        
            if len(tar_all_list) > 0: # Concatenate
                # Filter out empty DataFrames before concatenating
                tar_all_list = [df for df in tar_all_list if len(df) > 0]
                if len(tar_all_list) > 0:
                    tar_all = pd.concat(tar_all_list, ignore_index=True)
                else:
                    tar_all = pd.DataFrame(columns=main_columns)
            else:
                tar_all = pd.DataFrame(columns=main_columns)
            
            tar_all=tar_all.rename(columns={'inchikey':'db_inchikey'})
            comp_tar = self._postprocess_databases(tar_all)

            if len(comp_tar) > 0:
                # Add std compound ids to output
                comp_ids=comp_ids.rename(columns={'isomeric_smiles':'pc_iso_smiles','inchi':'pc_inchi',
                                                  'inchikey':'pc_inchikey','inchikey_fb':'pc_firstblock',
                                                  'iupac_name':'pc_iupac_name','cid':'pc_cid'})
                std_ids = comp_ids[['pc_iso_smiles', 'pc_inchi', 'pc_inchikey', 'pc_firstblock', 'pc_iupac_name', 'pc_cid']]
                # Select first nonNa value from each column
                std_ids = std_ids.apply(lambda col: std_ids[col.name].dropna().iloc[0] if not std_ids[col.name].dropna().empty else None)
                # Assign to each row the standard ids
                comp_tar = comp_tar.assign(**std_ids)       
                # Add input id to the results
                comp_tar['input_id'] = input_id

                # Reorder columns
                comp_tar = comp_tar[['input_id', 'pc_inchi', 'pc_inchikey', 'pc_firstblock', 'pc_iso_smiles', 
                                     'pc_iupac_name','pc_cid', 'db_inchikey', 'entrez', 'gene_type', 
                                     'hgnc_symbol', 'description',
                                    'pchembl_count', 'ave_pchembl', 'std_pchembl', 'src_count', 'pubchem',
                                    'chembl', 'bindingdb', 'stitch', 'ctd', 'dtc', 'otp', 'drugcentral',
                                    'drugbank']]

        return comp_tar, states
    

    def _postprocess_databases(self, tar_all) -> pd.DataFrame:
        # Remove non-protein coding interactions
        tar_all = tar_all[tar_all['gene_type']=='protein_coding']
        # Remove proteins with no Symbol or Entrez
        tar_all = tar_all.dropna(subset=['hgnc_symbol', 'entrez']).reset_index(drop=True)  
    
        # Group by BOTH protein AND compound (to keep stereoisomers separate)
        tar_list = tar_all[['hgnc_symbol', 'db_inchikey']].drop_duplicates().values.tolist()
    
        # Create output dataframe
        comp_tar = pd.DataFrame(columns=['entrez','hgnc_symbol','description','gene_type','db_inchikey',
                                   'pchembl_count','ave_pchembl','std_pchembl','src_count'] + 
                                   [source.lower() for source in self.sources])
    
        for index, (target, inchikey) in enumerate(tar_list):
            # Select rows for this protein-compound pair
            tar = tar_all.loc[(tar_all['hgnc_symbol'] == target) & (tar_all['db_inchikey'] == inchikey)]
        
            # Copy the biomart data to the output dataframe
            comp_tar.loc[index,'entrez'] = int(tar['entrez'].iloc[0])
            comp_tar.loc[index,'hgnc_symbol'] = tar['hgnc_symbol'].iloc[0]
            comp_tar.loc[index,'description'] = tar['description'].iloc[0]
            comp_tar.loc[index,'gene_type'] = tar['gene_type'].iloc[0]
            comp_tar.loc[index,'db_inchikey'] = inchikey  # Add the compound inchikey
        
            comp_tar = self._aggregate_pchembl(comp_tar, index, tar)
        
            # Count sources for this protein-compound pair
            comp_tar.loc[index, 'src_count'] = len(tar['datasource'].unique())
        
            # Create source binary matrix
            for source in self.sources:
                comp_tar.loc[index, source.lower()] = 1 if source in tar['datasource'].unique() else 0

        return comp_tar
