import numpy as np
import pandas as pd
from ..databases import *
from ..utils.identifiers import protein_identifiers
from .Pipeline import Pipeline

class Prot2Comp(Pipeline):

    def _update_args(self, pChEMBL_thres, stitch_stereo, dtc_mutated, dc_extra):
        self.database_args = {
            'pc': (pChEMBL_thres,),
            'chembl': (pChEMBL_thres,),
            'bdb': (pChEMBL_thres,),
            'stitch': (stitch_stereo,),
            'ctd': (),
            'dtc': (dtc_mutated, pChEMBL_thres,),
            'otp': (),
            'dc': (dc_extra, pChEMBL_thres,),
            'db': (),
        }      

    # Calls functions to collect data and merges all the data from the various sources together
    # The parameters are:
    #    - input_id - the protein id
    #    - pChEMBL_thresh - the minimum interaction pChEMBL value required to be added to the output file
    #    - stitch_stereo - to select whether to consider the specific compound stereochemistry or group all stereoisomers interactions from STITCH
    #    - dtc_mutated - to select whether also to consider interactions with mutated target proteins from DTC
    #    - dc_extra - to select whether to include possibly non-Homo sapiens interactions

    def prot_interactions(self, input_id, pChEMBL_thres=0, stitch_stereo=True, 
                          dtc_mutated=False, dc_extra=False):

        # Run interaction select with all databases selected
        prot_comp, state = self.prot_interactions_select(input_id, pChEMBL_thres=pChEMBL_thres, stitch_stereo=stitch_stereo, 
                                                         dtc_mutated=dtc_mutated, dc_extra=dc_extra)
        
        return prot_comp, state
            

    # Calls functions to collect data and merges all the data from the selected sources together
    # The parameters are:
    #    - input_id - the protein id
    #    - selected_dbs - underscore-separated string containing the names of the databases which are selected
    #    - pChEMBL_thresh - the minimum interaction pChEMBL value required to be added to the output file
    #    - stitch_stereo - to select whether to consider the specific compound stereochemistry or group all stereoisomers interactions from STITCH
    #    - dtc_mutated - to select whether also to consider interactions with mutated target proteins from DTC
    #    - dc_extra - to select whether to include possibly non-Homo sapiens interactions

    def prot_interactions_select(self, input_id, selected_dbs='pc_chembl_bdb_stitch_ctd_dtc_otp_dc_db', 
                                 pChEMBL_thres=0, stitch_stereo=True, dtc_mutated=False, dc_extra=False):

        self._update_args(pChEMBL_thres, stitch_stereo, dtc_mutated, dc_extra)

        prot_ids = protein_identifiers(input_id)

        main_columns=['inchi','inchikey','isomeric_smiles','iupac_name','pchembl_value','datasource']    

        prot_comp = pd.DataFrame(columns=main_columns)

        # Outputs the brief statement regarding the processing of each source 
        states=pd.DataFrame(columns=['prot'])
        states.loc[0, 'prot'] = input_id   
        
        if len(prot_ids) > 0:
            comp_all = pd.DataFrame(columns=main_columns)
            db_list = selected_dbs.split('_')
            for name, db in self.databases.items():      
                if name in db_list:
                    # Perform database search
                    db_comp, state, _ = db.compounds(prot_ids, *self.database_args.get(name, ()))
                    if len(db_comp) > 0:
                        # Keep only wanted columns
                        db_comp = db_comp.loc[:, main_columns]
                        # Concatenate the source
                        comp_all = pd.concat([comp_all, db_comp])
                    print(f'{name} done!')
                else:
                    state = f'Did not run {name}'
                # Add statement for this source
                states.loc[0, name] = state

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
                                    'inchi', 'inchikey', 'isomeric_smiles', 'iupac_name', 'pchembl_count',
                                    'ave_pchembl', 'std_pchembl', 'src_count', 'pubchem', 'chembl', 'bindingdb', 'stitch',
                                    'ctd', 'dtc', 'otp', 'drugcentral', 'drugbank']]

        return prot_comp, states
            

    def _postprocess_databases(self, comp_all):
        comp_all = comp_all.dropna(subset=['inchi']).reset_index(drop=True)

        # Find the compound with same inchikey in 1st and 2nd block. 
        comp_all['12Key'] = comp_all['inchikey'].apply(lambda x: x[:-1]) 
        comp_all['3Key'] = comp_all['inchikey'].apply(lambda x: x[-1:])

        # Make a list of the unique targets for the protein
        comp_list = comp_all['12Key'].unique()
        tar_comp=pd.DataFrame(columns=['inchi','inchikey','isomeric_smiles','iupac_name',
                                    'pchembl_count','ave_pchembl','src_count'] + 
                                    [source.lower() for source in self.sources]
                                    +['notes'])

        for index, compound in enumerate(comp_list):

            comp = comp_all.loc[comp_all['12Key'] == compound]
            comp1 = comp.drop_duplicates()
            if len(comp) > 0:
                # Make sure the neutral one is in the final output
                is_N_present = comp1['3Key'].str.contains('N') 
                indexes_with_N = comp1.loc[is_N_present].index
                if len(indexes_with_N) > 0:
                    tar_comp.loc[index,'inchi'] = comp1['inchi'][indexes_with_N[0]]
                    tar_comp.loc[index,'inchikey'] = comp1['inchikey'][indexes_with_N[0]]
                    tar_comp.loc[index,'isomeric_smiles'] = comp1['isomeric_smiles'][indexes_with_N[0]]
                    tar_comp.loc[index,'iupac_name'] = comp1['iupac_name'][indexes_with_N[0]]
                    # tar_comp.loc[index,'pchembl_value'] = comp1['pchembl_value'][indexes_with_N[0]]
                else:
                    tar_comp.loc[index,'inchi'] = comp1['inchi'][comp.index[0]]
                    tar_comp.loc[index,'inchikey'] = comp1['inchikey'][comp.index[0]]
                    tar_comp.loc[index,'isomeric_smiles'] = comp1['isomeric_smiles'][comp.index[0]]
                    tar_comp.loc[index,'iupac_name'] = comp1['iupac_name'][comp.index[0]]
                    # tar_comp.loc[index,'pchembl_value'] = comp1['pchembl_value'][comp.index[0]]
            
            tar_comp = self._aggregate_pchembl(tar_comp, index, comp)

            tar_comp.loc[index, 'src_count'] = len(comp['datasource'].unique())
            # Sets a source binary matrix to the dataframe to easily filter output later
            for source in self.sources:
                tar_comp.loc[index, source.lower()] = 1 if source in comp['datasource'].unique() else 0

        return tar_comp
    
