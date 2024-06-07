import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client
from ..servers.BiomartServer import BiomartServer
from ..servers.PubchemServer import PubChemServer 
from ..servers.ChEMBLServer import ChEMBLServer as chembl
from .Database import Database
from ..data_manager import *
import re

class ChEMBL(Database):

    def __init__(self, connection=None, database=None):
        if database is not None:
            self.data_manager = LocalManager(database)
        elif connection is not None:
            self.data_manager = SQLManager(connection, 'CHEMBL')
        else:
            funcs = {
                'Molecule ChEMBL ID' : self._retrieve_compounds,
                'Target ChEMBL ID' : self._retrieve_proteins
            }
            self.data_manager = APIManager(funcs) 

        # Mapping for API and local database column names
        self.column_names = {
            'molecule_chembl_id': 'Molecule ChEMBL ID',
            'activity_comment': 'Comment',
            'data_validity_comment': 'Data Validity Comment',
            'target_tax_id': 'Target Organism',
            'src_id': 'Source ID',
            'pchembl_value': 'pChEMBL Value',
            'target_chembl_id': 'Target ChEMBL ID'
        }

    def _filter_database(self, chembl_raw: pd.DataFrame, pChEMBL_thres: float):

        invalid_activities = ['Not Determined', 'Not Active', 'Not Evaluated', 'inactive', 'inconclusive', 'undetermined', 'No data']

        chembl_act = chembl_raw.loc[
                    # Remove all non-human
                    (chembl_raw['Target Organism'].isin(['9606', 'Homo sapiens'])) &  
                    # Remove various activity comments
                    (~chembl_raw['Comment'].isin(invalid_activities)) &  
                    # Only use interactions extracted from scientific literature
                    (chembl_raw['Source ID'] == 1) &  
                    # Only use interactions with no data validity comments
                    (chembl_raw['Data Validity Comment'].isnull()) &  
                    # Only use interactions with pchembl value
                    (chembl_raw['pChEMBL Value'].notnull())  
                ].drop_duplicates(subset=['Comment', 'Molecule ChEMBL ID', 'pChEMBL Value']).copy()
        chembl_act['pChEMBL Value'] = chembl_act['pChEMBL Value'].apply(lambda x: re.sub(r'[^0-9.]', '', x) if type(x) is str else x)
        chembl_act['pChEMBL Value'] = chembl_act['pChEMBL Value'].apply(pd.to_numeric, errors='coerce') 
        # Filter interactions by pChEMBL value
        chembl_act = chembl_act.loc[chembl_act['pChEMBL Value'] > pChEMBL_thres].reset_index(drop=True) 

        return chembl_act

    def interactions(self, input_comp: pd.DataFrame, chembl_ids: list=[], pChEMBL_thres=0):
        """
        Retrieves proteins from chembl API interacting with compound passed as input.

        Steps
        -----
        - Finds chembl id of all synonyms for the compound passed as input \\
        Constraints:
            - Only Homo Sapiens interactions
            - Remove various activity comments (Not Determined, Not Active, Not Evaluated, inactive, inconclusive, undetermined, No data)
            - Only interactions extracted from scientific literature
            - Only interactions with no data validity comments
            - Only interactions with pchembl value
        - Identifies target proteins using chembl API.
        - Uses Biomart to obtain proteins info (modified with data from first API search) to return,
        searching with chembl id.
        
        Parameters
        ----------
        input_comp : DataFrame
            Dataframe of input compounds from which interacting proteins are found
        chembl_ids : list
            Empty list that will contain all ChEMBL ids found for the input compound
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid

        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values:
            - entrez 
            - gene_type 
            - hgnc_symbol 
            - description 
            - datasource (ChEMBL) 
            - pchembl_value

        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all ChEMBL info about the input compound
        """
        
        columns = ['entrez','gene_type','hgnc_symbol','description','datasource', 'pchembl_value']
        # Create an empty DataFrame with the specified columns
        chembl_act = pd.DataFrame(columns=columns)
        chembl_raw = pd.DataFrame()
        # Find chembl ids from input compound
        chembl_ids.extend(chembl.identify_chembl_ids(input_comp))
        # Check if at least one id has been found
        if len(chembl_ids) > 0 and None not in chembl_ids:
 
            chembl_raw = self.data_manager.retrieve_raw_data('Molecule ChEMBL ID', chembl_ids)
            try:
                chembl_raw.reset_index()
            except:
                chembl_raw = pd.DataFrame()
            
            if len(chembl_raw) > 0:
                # Filter ChEMBL for high quality interactions. This follows: Bosc, N. et al. J. cheminformatics 11, 1–16 (2019)
                chembl_act = self._filter_database(chembl_raw, pChEMBL_thres)

                if len(chembl_act) > 0:
                    # Retrieve target type if not already available
                    if 'Target Type' not in chembl_act.columns:
                        chembl_act = self._retrieve_target_type(chembl_act)

                    # Keep only relevant columns
                    chembl_act = chembl_act[['Molecule ChEMBL ID', 'Comment', 'Data Validity Comment', 'Target Organism',
                                             'Source ID', 'pChEMBL Value', 'Target ChEMBL ID', 'Target Type']]
                    # Only use int. w/proteins, removes cell lines and RNA
                    chembl_act = chembl_act.loc[chembl_act['Target Type']=='SINGLE PROTEIN'] 
                    chembl_act = chembl_act.rename(columns={'pChEMBL Value': 'pchembl_value'}).reset_index(drop=True)
                    
                    if len(chembl_act) > 0:
                        # Unify gene identifiers by Converting ChEMBL with biomart
                        ensembl = BiomartServer()
                        # ChEMBL uses chembl ids
                        input_type='chembl' 
                        attributes = ['chembl', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                        names = ['chembl','entrez','gene_type','hgnc_symbol','description']
                        chembl_targets=pd.DataFrame(columns=names)
                        
                        input_genes=list(chembl_act['Target ChEMBL ID'])
                        
                        chembl_targets=ensembl.subset_search(input_type, input_genes, attributes, names)

                        # For each compound, assign specific biomart column values to the ones from the original chembl database
                        for index, row in chembl_act.iterrows():
                            S1 = chembl_targets.loc[chembl_targets['chembl']==row['Target ChEMBL ID']]
                            if len(S1) > 0:
                                chembl_act.loc[index, 'entrez'] = S1['entrez'].iloc[0]
                                chembl_act.loc[index, 'gene_type'] = S1['gene_type'].iloc[0]
                                chembl_act.loc[index, 'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                                chembl_act.loc[index, 'description'] = S1['description'].iloc[0]
                            else:
                                chembl_act.loc[index, 'entrez'] = None
                                chembl_act.loc[index, 'gene_type'] = None
                                chembl_act.loc[index, 'hgnc_symbol'] = None
                                chembl_act.loc[index, 'description'] = None  
                                Des='Failed to convert the gene ID'
                        chembl_act['datasource'] = 'ChEMBL'
                        statement = 'completed'
                    else:
                        statement = 'Filter reduced interactions to 0'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'No ChEMBL ids found'

        return chembl_act, statement, chembl_raw
    
    def compounds(self, input_protein: pd.DataFrame, pChEMBL_thres: float=0):
        """
        Retrieves compounds from chEMBL database interacting with proteins passed as input.

        Steps
        -----
        - Makes an API call to chEMBL API to retrieve compounds interacting with input protein \\
        Constraints:
            - Only Homo Sapiens interactions
            - Remove various activity comments (Not Determined, Not Active, Not Evaluated, inactive, inconclusive, undetermined, No data)
            - Only interactions extracted from scientific literature
            - Only interactions with no data validity comments
            - Only interactions with pchembl value
        - Uses Pubchempy to retrieve the compound info (with additional data from chembl) to return, 
        searching with chEMBL ID or with inchi (obtained from chembl API) if the first search is unsuccessful.
        
        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compound are found
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid

        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values:
            - inchi
            - inchikey
            - isomeric smiles
            - iupac name
            - datasource (ChEMBL)
            - pchembl value
            - notes (NaN)
        """

        columns = ['inchi','inchikey','isomeric_smiles','iupac_name','datasource','pchembl_value']
        # Create an empty DataFrame with the specified columns
        chembl_c1 = pd.DataFrame(columns=columns) 
        chembl_raw = pd.DataFrame()
        # Drop duplicated and null value of chembl_id, to make sure the first line of the dataset has the chembl ID
        input_protein = input_protein.dropna(subset=['chembl'])
        input_protein = input_protein.drop_duplicates(subset='chembl').reset_index(drop=True)        
        # Check if there are any input proteins remaining
        if len(input_protein) > 0:
            input_protein_id = input_protein['chembl'].iloc[0]

            chembl_raw = self.data_manager.retrieve_raw_data('Target ChEMBL ID', input_protein_id)
            # Check if the API returned any interacting compound
            if len(chembl_raw) > 0:
                # Filter ChEMBL for high quality interactions.
                chembl_act = self._filter_database(chembl_raw, pChEMBL_thres)

                if len(chembl_act) > 0:

                    pc = PubChemServer()

                    if 'CID' in chembl_act.columns:
                        # Retrieve compounds using cids
                        compounds = self._pubchem_search_cid(chembl_act, columns, pc)
                        
                        if len(compounds) > 0:
                            chembl_info = chembl_act[['CID', 'Molecule ChEMBL ID', 'pChEMBL Value', 'Target ChEMBL ID']].\
                                            rename(columns={'Molecule ChEMBL ID': 'id', 'pChEMBL Value': 'pchembl_value', 
                                                            'Target ChEMBL ID': 'target_chembl_id'}).drop_duplicates()
                            chembl_c1 = pd.merge(compounds, chembl_info, on='CID', how='left')
                        
                    else:
                        # Retrieve compounds using chembl ids and inchis if formers miss
                        compounds = self._pubchem_search_chembl(chembl_act, 'Molecule ChEMBL ID', columns, pc)

                        if len(compounds) > 0:
                            # Update pchembl value with that of filtered compounds
                            chembl_info = chembl_act[['Molecule ChEMBL ID', 'pChEMBL Value', 'Target ChEMBL ID']].\
                                            rename(columns={'Molecule ChEMBL ID': 'id', 'pChEMBL Value': 'pchembl_value', 
                                                            'Target ChEMBL ID': 'target_chembl_id'}).drop_duplicates()
                            chembl_c1 = pd.merge(compounds, chembl_info, on='id', how='left')

                    # Check if at least one compound has been found in total
                    if len(chembl_c1) > 0:
                        # Store into dataframe additional data
                        chembl_c1.loc[:, 'datasource'] = 'ChEMBL'
                        chembl_c1.loc[:, 'notes'] = np.nan
                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using Pubchem and ChEMBL'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein doesn\'t contain chembl id'

        return chembl_c1, statement, chembl_raw
    

    def _retrieve_compounds(self, chembl_ids):
            chembl_raw = pd.DataFrame(columns=['entrez','gene_type','hgnc_symbol','description','datasource', 'pchembl_value'])
            # API to Obtain all the targets for the list of chembl ids and their reported activities
            for id in chembl_ids: 
                # Use the chembl webresource client to find interactions
                activities = new_client.activity.filter(molecule_chembl_id=id)\
                .only(['molecule_chembl_id','activity_comment','data_validity_comment','target_tax_id',\
                       'src_id','pchembl_value','target_chembl_id'])
                if len(activities) != 0:
                    mol_act = pd.DataFrame(activities)
                    chembl_raw = pd.concat([chembl_raw, mol_act])
            
            chembl_raw = chembl_raw.rename(columns=self.column_names)

            return chembl_raw
    
    def _retrieve_target_type(self, chembl_act):
        # Filter ChEMBL based on target_type, looking to only use proteins
        target = new_client.target
        dat = []
        # Use the chembl webresource client to find target_type
        for _, row in chembl_act.iterrows(): 
            res = target.filter(target_chembl_id=row['Target ChEMBL ID']).only(['target_type'])
            dat.append(res[0]['target_type'])
        chembl_act['Target Type'] = dat

        return chembl_act

    def _retrieve_proteins(self, input_protein_id):
        # Use package of Chembl API to get the compound that have an activity with the protein
        activities = new_client.activity.filter(target_chembl_id=input_protein_id)\
                    .only(['molecule_chembl_id','activity_comment','data_validity_comment','target_tax_id',\
                            'src_id','pchembl_value','target_chembl_id'])
        # Create dataframe from API output
        chembl_raw = pd.DataFrame.from_dict(activities)
        chembl_raw = chembl_raw.rename(columns=self.column_names)

        return chembl_raw