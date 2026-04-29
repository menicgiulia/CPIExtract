'''Retrieve compound interactions for proteins provided as input.'''

import numpy as np
import pandas as pd
from ..databases import *
from ..utils.identifiers import protein_identifiers
from ..servers.PubChemServer import PubChemServer
from .Pipeline import Pipeline

class Prot2Comp(Pipeline):
    '''Retrieve compound interactions for proteins provided as input.'''

    def _update_args(self, pChEMBL_thres: float, dtc_mutated: bool, dc_extra: bool, merge_stereoisomers: bool, verbose: bool) -> None:
        self.database_args = {
            'pc': (pChEMBL_thres, merge_stereoisomers,),
            'chembl': (pChEMBL_thres, merge_stereoisomers,),
            'bdb': (pChEMBL_thres, merge_stereoisomers,),
            'stitch': (merge_stereoisomers,),
            'ctd': (merge_stereoisomers,),
            'dtc': (dtc_mutated, pChEMBL_thres, merge_stereoisomers,),
            'otp': (merge_stereoisomers,),
            'dc': (dc_extra, pChEMBL_thres, merge_stereoisomers,),
            'db': (merge_stereoisomers,),
        }      

    # Calls functions to collect data and merges all the data from the various sources together
    # The parameters are:
    #    - input_id - the protein id
    #    - pChEMBL_thresh - the minimum interaction pChEMBL value required to be added to the output file
    #    - dtc_mutated - to select whether also to consider interactions with mutated target proteins from DTC
    #    - dc_extra - to select whether to include possibly non-Homo sapiens interactions
    #    - merge_stereoisomers - to select whether data collected is stereo-specific or generic to a structure

    def prot_interactions(self, input_id: int|str, pChEMBL_thres: float=0, dtc_mutated: bool=False, dc_extra: bool=False, 
                          merge_stereoisomers: bool=False,verbose: bool=False) -> tuple[pd.DataFrame, pd.DataFrame]:

        # Run interaction select with all databases selected
        prot_comp, state = self.prot_interactions_select(input_id, pChEMBL_thres=pChEMBL_thres, dtc_mutated=dtc_mutated, 
                                                        dc_extra=dc_extra, merge_stereoisomers=merge_stereoisomers,verbose=verbose)
        
        return prot_comp, state
            

    # Calls functions to collect data and merges all the data from the selected sources together
    # The parameters are:
    #    - input_id - the protein id
    #    - selected_dbs - underscore-separated string containing the names of the databases which are selected
    #    - pChEMBL_thresh - the minimum interaction pChEMBL value required to be added to the output file
    #    - dtc_mutated - to select whether also to consider interactions with mutated target proteins from DTC
    #    - dc_extra - to select whether to include possibly non-Homo sapiens interactions
    #    - merge_stereoisomers - to select whether data collected is stereo-specific or generic to a structure

    def prot_interactions_select(self, input_id: int|str, selected_dbs: str='pc_chembl_bdb_stitch_ctd_dtc_otp_dc_db', 
                                 pChEMBL_thres: float=0, dtc_mutated: bool=False, dc_extra: bool=False, merge_stereoisomers: bool=False,
                                 verbose: bool=False) -> tuple[pd.DataFrame, pd.DataFrame]:

        self._update_args(pChEMBL_thres, dtc_mutated, dc_extra, merge_stereoisomers,verbose)

        prot_ids = protein_identifiers(input_id,gene_server=self.gene_server)

        main_columns=['inchi','inchikey','isomeric_smiles','iupac_name','pchembl_value','datasource']    

        prot_comp = pd.DataFrame(columns=main_columns)

        # Outputs the brief statement regarding the processing of each source 
        states=pd.DataFrame(columns=['prot'])
        states.loc[0, 'prot'] = input_id   
        
        if len(prot_ids) > 0:
            comp_all_list = []
            db_list = selected_dbs.split('_')
            for name, db in self.databases.items():
                if name in db_list:
                    # Perform database search
                    result, state, _ = db.compounds(prot_ids, *self.database_args.get(name, ()))
                    if len(result) > 0:
                        # Ensure all required columns exist before selecting
                        for col in main_columns:
                            if col not in result.columns:
                                result[col] = None
                        
                        # Keep only wanted columns
                        result = result.loc[:, main_columns]
                        # Append to list
                        comp_all_list.append(result)
                    if verbose:
                        print(f'{name} done!')
                else:
                    state = f'Did not run {name}'
                # Add statement for this source
                states.loc[0, name] = state
        
            # Concatenate once AFTER loop (OUTSIDE for loop - notice indentation)
            if len(comp_all_list) > 0:
                comp_all = pd.concat(comp_all_list, ignore_index=True)
            else:
                comp_all = pd.DataFrame(columns=main_columns)

            prot_comp = self._postprocess_databases(comp_all)
                        
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
                prot_comp = prot_comp[['input_id', 'entrez', 'gene_type', 'hgnc_symbol', 'description',
                                    'inchi', 'inchikey', 'inchikey_fb', 'isomeric_smiles', 'iupac_name', 'pchembl_count',
                                    'ave_pchembl', 'std_pchembl', 'src_count', 'pubchem', 'chembl', 'bindingdb', 'stitch',
                                    'ctd', 'dtc', 'otp', 'drugcentral', 'drugbank']]

        return prot_comp, states
            

    def _postprocess_databases(self, comp_all: pd.DataFrame) -> pd.DataFrame:
        # Remove compounds without InChI
        comp_all = comp_all.dropna(subset=['inchi']).reset_index(drop=True)
    
        # Add first-block column for reference
        pc = PubChemServer()
        comp_all['inchikey_fb'] = comp_all['inchikey'].apply(
            lambda x: pc.get_inchikey_first_block(x) if pd.notna(x) else None
        )
    
        # Deduplicate by full InChIKey (not by first-block)
        comp_list = comp_all['inchikey'].unique()

        tar_comp = pd.DataFrame(columns=['inchi', 'inchikey', 'inchikey_fb', 'isomeric_smiles', 'iupac_name',
                                     'pchembl_count', 'ave_pchembl', 'std_pchembl', 'src_count'] + 
                                     [source.lower() for source in self.sources])
    
        for index, inchikey in enumerate(comp_list):
            comp = comp_all.loc[comp_all['inchikey'] == inchikey]
        
            if len(comp) > 0:
                # Take first occurrence for compound identifiers
                tar_comp.loc[index, 'inchi'] = comp['inchi'].iloc[0]
                tar_comp.loc[index, 'inchikey'] = comp['inchikey'].iloc[0]
                tar_comp.loc[index, 'inchikey_fb'] = comp['inchikey_fb'].iloc[0]
                tar_comp.loc[index, 'isomeric_smiles'] = comp['isomeric_smiles'].iloc[0]
                tar_comp.loc[index, 'iupac_name'] = comp['iupac_name'].iloc[0]
            
                # Aggregate pchembl values
                tar_comp = self._aggregate_pchembl(tar_comp, index, comp)
            
                # Count sources
                tar_comp.loc[index, 'src_count'] = len(comp['datasource'].unique())
            
                # Create source binary matrix
                for source in self.sources:
                    tar_comp.loc[index, source.lower()] = 1 if source in comp['datasource'].unique() else 0
    
        return tar_comp
    
