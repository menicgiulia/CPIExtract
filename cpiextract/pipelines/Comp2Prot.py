import numpy as np
import pandas as pd
from ..utils.identifiers import compound_identifiers
from ..databases import *
from .Pipeline import Pipeline

class Comp2Prot(Pipeline):

    def _update_args(self, pChEMBL_thres, stitch_stereo, otp_biblio, dtc_mutated, dc_extra, chembl_ids):
        self.database_args = {
            'pc': (pChEMBL_thres,),
            'chembl': (chembl_ids, pChEMBL_thres,),
            'bdb': (pChEMBL_thres,),
            'stitch': (stitch_stereo,),
            'ctd': (),
            'dtc': (chembl_ids, dtc_mutated, pChEMBL_thres,),
            'otp': (chembl_ids, otp_biblio,),
            'dc': (dc_extra, pChEMBL_thres,),
            'db': (),
        }       


    # Calls functions to collect data and merges all the data from the various sources together  
    # The parameters are:
    #    - input_id - the compound id
    #    - pChEMBL_thresh - the minimum interaction pChEMBL value required to be added to the output file
    #    - stitch_stereo - to select whether to consider the specific compound stereochemistry or group all stereoisomers interactions from STITCH
    #    - otp_biblio - to select whether to include the *bibliography* data from OTP. This parameter is only available for Comp2Prot as OTP provides only known drug interactions for proteins
    #    - dtc_mutated - to select whether also to consider interactions with mutated target proteins from DTC
    #    - dc_extra - to select whether to include possibly non-Homo sapiens interactions

    def comp_interactions(self, input_id, pChEMBL_thres=0, stitch_stereo=True, 
                          otp_biblio=False, dtc_mutated=False, dc_extra=False):
        
        # Run interaction select with all databases selected
        comp_tar, state = self.comp_interactions_select(input_id, pChEMBL_thres=pChEMBL_thres, 
                                                    stitch_stereo=stitch_stereo, otp_biblio=otp_biblio, 
                                                    dtc_mutated=dtc_mutated, dc_extra=dc_extra)
        return comp_tar, state

    # Calls functions to collect data and merges all the data from the selected sources together
    # The parameters are:
    #    - input_id - the compound id
    #    - selected_dbs - underscore-separated string containing the names of the databases which are selected
    #    - pChEMBL_thresh - the minimum interaction pChEMBL value required to be added to the output file
    #    - stitch_stereo - to select whether to consider the specific compound stereochemistry or group all stereoisomers interactions from STITCH
    #    - otp_biblio - to select whether to include the *bibliography* data from OTP. This parameter is only available for Comp2Prot as OTP provides only known drug interactions for proteins
    #    - dtc_mutated - to select whether also to consider interactions with mutated target proteins from DTC
    #    - dc_extra - to select whether to include possibly non-Homo sapiens interactions

    def comp_interactions_select(self, input_id, selected_dbs='pc_chembl_bdb_stitch_ctd_dtc_otp_dc_db', 
                                 pChEMBL_thres=0, stitch_stereo=True, otp_biblio=False, dtc_mutated=False, dc_extra=False):
        #pChEMBL_thres default 0: Will only interactions without activity data in sources with activity data present
        #stitch_stereo default True: True=stereo specific False=non-specific stereochemistry (True ensures exact match to input_id)
        #otp_biblio default False: False=mechanism data only True=mechanism and bibliography (False ensures real interactions, avoids keyword search matches within abstracts that can lead to false associations that are not real interactions)
        
        chembl_ids = []
        self._update_args(pChEMBL_thres, stitch_stereo, otp_biblio, dtc_mutated, dc_extra, chembl_ids)

        main_columns = ['entrez','gene_type','hgnc_symbol','description','pchembl_value','datasource']
        comp_tar = pd.DataFrame(columns=main_columns)
        #Collect all the data through the other functions
        comp_ids = compound_identifiers(input_id) #find comp_ids
        #Outputs the brief statement regarding the processing of each source 
        states = pd.DataFrame(columns=['comp'])
        states.loc[0, 'comp'] = input_id

        if len(comp_ids) > 0:
            tar_all = pd.DataFrame(columns=main_columns)
            db_list = selected_dbs.split('_')
            for name, db in self.databases.items():
                if name in db_list:
                    # Perform database search
                    result, state, _ = db.interactions(comp_ids, *self.database_args.get(name, ()))
                    if len(result) > 0:
                        # Keep only wanted columns
                        result = result.loc[:, main_columns]
                        # Concatenate the source
                        tar_all = pd.concat([tar_all, result])
                    print(f'{name} done!')
                else:
                    state = f'Did not run {name}'
                # Add statement for this source
                states.loc[0, name] = state
            
            comp_tar = self._postprocess_databases(tar_all)

            if len(comp_tar) > 0:
                # Add std compound ids to output
                std_ids = comp_ids[['isomeric_smiles', 'inchi', 'inchikey', 'iupac_name']]
                # Select first nonNa value from each column
                std_ids = std_ids.apply(lambda col: std_ids[col.name].dropna().iloc[0] if not std_ids[col.name].dropna().empty else None)
                # Assign to each row the standard ids
                comp_tar = comp_tar.assign(**std_ids)       
                # Add input id to the results
                comp_tar['input_id'] = input_id

                # Reorder columns
                comp_tar = comp_tar[['input_id', 'inchi', 'inchikey', 'isomeric_smiles', 'iupac_name', 
                                    'entrez', 'gene_type', 'hgnc_symbol', 'description',
                                    'pchembl_count', 'ave_pchembl', 'std_pchembl', 'src_count', 'pubchem',
                                    'chembl', 'bindingdb', 'stitch', 'ctd', 'dtc', 'otp', 'drugcentral',
                                    'drugbank']]

        return comp_tar, states
    

    def _postprocess_databases(self, tar_all):

        # Remove non-protein coding interactions
        tar_all = tar_all[tar_all['gene_type']=='protein_coding']
        # Remove proteins with no Symbol or Entrez
        tar_all = tar_all.dropna(subset=['hgnc_symbol', 'entrez']).reset_index(drop=True)  
        # Make a list of the unique targets for the compound
        tar_list = tar_all['hgnc_symbol'].unique() 
        
        #Create an output dataframe 
        comp_tar=pd.DataFrame(columns=['entrez','hgnc_symbol','description','gene_type','pchembl_count',
                                       'ave_pchembl','src_count'] + 
                                       [source.lower() for source in self.sources])
        
        for index, target in enumerate(tar_list):
            tar = tar_all.loc[tar_all['hgnc_symbol'] == target] #Select target
            
            # Copy the biomart data to the output dataframe
            comp_tar.loc[index,'entrez'] = int(tar['entrez'][tar.index[0]])
            comp_tar.loc[index,'hgnc_symbol'] = tar['hgnc_symbol'][tar.index[0]]
            comp_tar.loc[index,'description'] = tar['description'][tar.index[0]]
            comp_tar.loc[index,'gene_type'] = tar['gene_type'][tar.index[0]]
            
            comp_tar = self._aggregate_pchembl(comp_tar, index, tar)
            
            # Sets a source binary matrix to the dataframe to easily filter output later
            comp_tar.loc[index, 'src_count'] = len(tar['datasource'].unique())
            # Sets a source binary matrix to the dataframe to easily filter output later
            for source in self.sources:
                comp_tar.loc[index, source.lower()] = 1 if source in tar['datasource'].unique() else 0


        return comp_tar
