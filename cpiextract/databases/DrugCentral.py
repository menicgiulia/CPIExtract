import numpy as np
import pandas as pd
from ..servers.BiomartServer import BiomartServer
from ..servers.PubchemServer import PubChemServer
from .Database import Database
from ..data_manager import *
import time


class DrugCentral(Database):

    def __init__(self, connection=None, database=None):
        # if not connection and not database:
        #     raise ValueError('Either SQL connection or database should be not None')
        if database is not None:
            self.data_manager = LocalManager(database)
        else:
            self.data_manager = SQLManager(connection, 'DC')

    def _filter_database(self, dc_raw: pd.DataFrame, dc_extra: bool):

        dc_filt = dc_raw.dropna(subset=['ACTION_TYPE'])
        dc_filt = dc_filt.loc[(dc_filt['ORGANISM']=='Homo sapiens') &
                             (dc_filt['ACT_TYPE'].isin(['IC50', 'Ki', 'EC50', 'Kd', 'AC50'])) &
                             (dc_filt['ACT_SOURCE'] != 'UNKNOWN')]\
                                .drop_duplicates(ignore_index=True).copy()
        if not dc_extra:
            invalid_action_types = ['PHARMACOLOGICAL CHAPERONE', 'RELEASING AGENT']
            invalid_classes = ['CD molecules', 'RNA', 'Unclassified', 'Viral envelope protein', 'Polyprotein']
            dc_filt = dc_filt.loc[(~dc_filt['ACTION_TYPE'].isin(invalid_action_types)) & 
                                  (~dc_filt['TARGET_CLASS'].isin(invalid_classes))].copy()
        return dc_filt
    

    def _compute_pchembl(self, dc_dat: pd.DataFrame, pChEMBL_thres: float):

        dc_dat['pchembl_value'] = dc_dat['ACT_VALUE'].apply(lambda x: -np.log10(x / 1e9))
        dc_act = dc_dat.rename(columns={'ACT_TYPE': 'notes'})

        # Filter interactions by pChEMBL value
        dc_act = dc_act.loc[dc_act['pchembl_value'] > pChEMBL_thres].reset_index(drop=True) 

        return dc_act

    def interactions(self, input_comp: pd.DataFrame, dc_extra: bool=False, pChEMBL_thres: float=0):
        """
        Retrieves proteins from DrugCentral database interacting with compound passed as input.

        Steps
        -----
        - Filters dc database to obtain proteins interacting with input compound \\
        Constraints:
            - Only Homo sapiens interactions
            - Activity Type must be one of the following: IC50, Ki, EC50, Kd, AC50
            - Action Type must not be: PHARMACOLOGICAL CHAPERONE, RELEASING AGENT
            - Target Class must not be: CD molecules, RNA, Unclassified, Viral envelope protein, Polyprotein
        - Uses activity column to compute pchembl values.
        - Uses Biomart to obtain proteins info (modified with data from original dc database) to return,
        searching with UniProt ID.
        
        Parameters
        ----------
        input_comp : DataFrame
            Dataframe of input compound data from which interacting proteins are found
        dc_extra: bool
            bool to select whether to include possibly non-Homo sapiens interactions
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid
            
        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values: \\
            entrez, gene_type, hgnc_symbol, description, datasource (DrugCentral), pchembl_value
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all DrugCentral info about the input compound
        """

        columns = ['entrez','gene_type','hgnc_symbol','description','pchembl_value','datasource']
        # Create an empty DataFrame with the specified columns
        dc_act = pd.DataFrame(columns=columns)
        dc_raw = pd.DataFrame()
        # Drop null values of inchikey, to make sure the first line of the dataframe has the inchi key 
        input_comp = input_comp.dropna(subset=['inchikey']).reset_index(drop=True)
        # Check if there are any input compounds remaining  
        if len(input_comp) > 0:
            # Select compound from BindingDB data based on inchi key
            input_comp_id = input_comp['inchikey'][0]
            dc_raw = self.data_manager.retrieve_raw_data('InChIKey', input_comp_id)
            if len(dc_raw) > 0:
                # Filter database
                dc_filt = self._filter_database(dc_raw, dc_extra)
                dc_filt.drop(['SMILES','InChI'], axis=1, inplace=True)
                # Compute pchembl values
                dc_act = self._compute_pchembl(dc_filt, pChEMBL_thres)
                
                if len(dc_act) > 0:
                    # Unify gene identifiers by Converting BindingDB with biomart
                    ensembl = BiomartServer()
                    # BindingDB uses Uniprot IDs
                    input_type='uniprotswissprot' 
                    attributes = ['uniprotswissprot', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                    names = ['uniprot','entrez','gene_type','hgnc_symbol','description']
                    dc_targets = pd.DataFrame(columns=names)
                    
                    input_genes = list(dc_act['ACCESSION'])
                    
                    dc_targets = ensembl.subset_search(input_type, input_genes, attributes, names)
    
                    # For each protein, assign specific biomart column values to the ones from the original dc database
                    for index, row in dc_act.iterrows():
                        # Find the compound using uniprot
                        S1 = dc_targets.loc[dc_targets['uniprot'] == row['ACCESSION']]
                        if len(S1) > 0:
                            dc_act.loc[index, 'entrez'] = S1['entrez'].iloc[0]
                            dc_act.loc[index, 'gene_type'] = S1['gene_type'].iloc[0]
                            dc_act.loc[index, 'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                            dc_act.loc[index, 'description'] = S1['description'].iloc[0]
                        else:
                            dc_act.loc[index, 'entrez'] = None
                            dc_act.loc[index, 'gene_type'] = None
                            dc_act.loc[index, 'hgnc_symbol'] = None
                            dc_act.loc[index, 'description'] = None
                            Des='Failed to convert the gene ID'
                    dc_act['datasource'] = 'DrugCentral'
                    statement = 'completed'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input compound doesn\'t contain inchi key'

        return dc_act, statement, dc_raw
    

    def compounds(self, input_protein: pd.DataFrame, dc_extra: bool=False, pChEMBL_thres: float=0):
        """
        Retrieves compounds from DrugCentral database interacting with proteins passed as input.

        Steps
        -----
        - Filters dc database to obtain compounds interacting with input proteins:
        Constraints:
            - Only Homo sapiens interactions
            - Activity Type must be one of the following: IC50, Ki, EC50, Kd, AC50
            - Action Type must not be: PHARMACOLOGICAL CHAPERONE, RELEASING AGENT
            - Target Class must not be: CD molecules, RNA, Unclassified, Viral envelope protein, Polyprotein
        - Uses activity column to compute pchembl values. \\
        - Uses Pubchempy to obtain compounds info (modified with data from original dc database) to return,
        searching with InChiKey
        
        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compound are found
        dc_extra: bool
            bool to select whether to include possibly non-Homo sapiens interactions
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid       

        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values: \\
            inchi, inchikey, isomeric_smiles, iupac_name, datasource (DrugCentral), pchembl_value, notes (activity type)
        """

        columns = ['inchi','inchikey','isomeric_smiles','iupac_name','datasource','pchembl_value']
        # Create an empty DataFrame with the specified columns
        dc_c1 = pd.DataFrame(columns=columns)
        dc_raw = pd.DataFrame()
        # Drop duplicated and null values of uniprot, to make sure the first line of the dataframe has the uniprot ID
        input_protein= input_protein.dropna(subset=['uniprot'])
        input_protein= input_protein.drop_duplicates(subset='uniprot').reset_index(drop=True) 
        # Check if there are any input proteins remaining       
        if len(input_protein) > 0:
            input_protein_id=input_protein['uniprot'].iloc[0]
            # Select only compounds interacting with the input protein
            dc_raw = self.data_manager.retrieve_raw_data('ACCESSION', input_protein_id)
            if len(dc_raw):
                # Filter database
                dc_c = self._filter_database(dc_raw, dc_extra)
                
                # Compute pchembl values
                dc_act = self._compute_pchembl(dc_c, pChEMBL_thres)
                # Check if at least one interacting compound has been found
                if len(dc_act) > 0:
                    pc = PubChemServer()

                    if 'CID' in dc_act.columns:
                        # Retrieve compounds using cids
                        compounds = self._pubchem_search_cid(dc_act, columns, pc)
                        
                        if len(compounds) > 0:
                            # Add additional values from activity dataframe
                            dc_info = dc_act[['CID', 'ACCESSION', 'TARGET_CLASS', 'ACT_COMMENT',
                                                'pchembl_value', 'notes']].\
                                        rename(columns={'ACCESSION': 'uniprotid'}).\
                                                        drop_duplicates()
                            
                            dc_c1 = pd.merge(compounds, dc_info, on='CID', how='left')
                    else:
                        # Filter only columns from pubchempy to return
                        selected_columns = pc.get_columns(columns[:-2])
                        compounds = pd.DataFrame()
                        ids = list(dc_act['InChIKey'].unique())
                        for id in ids:
                            try:
                                comp = pc.get_compounds(id, selected_columns, namespace='inchikey')
                                compounds = pd.concat([compounds, comp])
                            except:
                                None
                            time.sleep(0.5)
                        # Check if at least a match has been found
                        if len(compounds) > 0:
                            dc_c1 = compounds.rename(columns=pc.properties)

                            # Add additional values from activity dataframe
                            dc_info = dc_act[['InChIKey', 'ACCESSION', 'TARGET_CLASS', 'ACT_COMMENT',
                                                'pchembl_value', 'notes']].\
                                        rename(columns={'InChIKey': 'inchikey', 
                                                        'ACCESSION': 'uniprotid'}).\
                                                        drop_duplicates()
                            
                            dc_c1 = pd.merge(dc_c1, dc_info, on='inchikey', how='left')

                    if len(dc_c1) > 0:
                        dc_c1['datasource'] = 'DrugCentral'
                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using PubChem'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein doesn\'t contain UniProt ID'
        return dc_c1, statement, dc_raw
