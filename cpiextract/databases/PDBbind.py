import pandas as pd
import numpy as np
from ..servers.PubchemServer import PubChemServer
from ..servers.BiomartServer import BiomartServer
from .Database import Database
from ..data_manager import *
import requests

class PDBbind(Database):

    def __init__(self, connection=None, database=None):
        # if not connection and not database:
        #     raise ValueError('Either SQL connection or database should be not None')
        if database is not None:
            self.data_manager = LocalManager(database)
        else:
            self.data_manager = SQLManager(connection, 'PDBbind')

        # Set the endpoint URL
        self.api_url = 'https://search.rcsb.org/rcsbsearch/v2/query'
        self.api_query = {
            "query": {
                "type": "terminal",
                "service": "chemical",
                "parameters": {
                "value": "",
                "type": "descriptor",
                "descriptor_type": "InChI",
                "match_type": "graph-exact"
                }
            },
            "return_type": "mol_definition"
            }
        
        # Set the GraphQL endpoint URL
        self.graphql_url = 'https://data.rcsb.org/graphql'


    def _filter_database(self, PDB_raw: pd.DataFrame, pChEMBL_thres: float) -> pd.DataFrame:
                                # Filter to Homo Sapiens only interactions
        PDB_act = PDB_raw.loc[(PDB_raw['organism'] == 'Homo sapiens') & 
                                # Filter interactions by pChEMBL value
                                PDB_raw['pchembl'] > pChEMBL_thres].reset_index(drop=True)
    
        return PDB_act

    def interactions(self, input_comp: pd.DataFrame, pChEMBL_thres: float):
        """
        Retrieves proteins from PDB database interacting with compound passed as input.

        Steps
        -----
        - Finds matches for compounds with the same inchi as the one passed as input
        - Finds all proteins interacting with input compound \\
        Constraints:
            - Only Homo Sapiens interactions
        - Uses Biomart to obtain proteins info (modified with data from from original PDB database) to return,
        searching with Uniprot ID.
        
        Parameters
        ----------
        input_comp : DataFrame
            Dataframe of input compounds from which interacting proteins are found
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid

        Returns
        -------
        DataFrame
        DataFrame
            Dataframe of interacting proteins, containing the following values:
            - entrez 
            - gene_type 
            - hgnc_symbol 
            - description 
            - datasource (PDB) 
            - pchembl_value
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all PDB info about the input compound
        """

        columns = ['entrez','gene_type','hgnc_symbol','description','pchembl_value','datasource']
        # Create an empty DataFrame with the specified columns
        PDB_act = pd.DataFrame(columns=columns)
        PDB_raw = pd.DataFrame()
        input_comp = input_comp.dropna(subset=['inchi'])
        if len(input_comp) > 0:
            input_comp_id = input_comp['inchi'][0]

            pdb_id = self._retrieve_comp_id(input_comp_id)

            PDB_raw = self.data_manager.retrieve_raw_data('ligand ID', pdb_id)

            # Check if at least one match has been found
            if len(PDB_raw) > 0:
                # Extract target organism data
                PDB_act = self._retrieve_prot_info(PDB_raw)
                # Filter database for human target interactions
                PDB_act = self._filter_database(PDB_act, pChEMBL_thres)
                
                # Check if at least one interaction has been found
                if len(PDB_act) > 0:
                    # Unify gene identifiers by Converting PDB with Biomart
                    ensembl = BiomartServer()
                    attributes = ['uniprotswissprot', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                    names = ['uniprot', 'entrez','gene_type','hgnc_symbol','description']
                    # Use entrez for conversion
                    input_type='uniprotswissprot' 
                    
                    PDB_targets=pd.DataFrame(columns=names)
                    
                    input_genes=list(PDB_act['Uniprot ID'])
                    
                    PDB_targets=ensembl.subset_search(input_type, input_genes, attributes, names)

                    # For each compound, assign specific biomart column values to the ones from the original PDB database
                    for index, row in PDB_act.iterrows():
                        S1 = PDB_targets.loc[PDB_targets['uniprot']==row['Uniprot ID']]
                        if len(S1) > 0:
                            PDB_act.loc[index,'entrez'] = S1['entrez'].iloc[0]
                            PDB_act.loc[index,'gene_type'] = S1['gene_type'].iloc[0]
                            PDB_act.loc[index,'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                            PDB_act.loc[index,'description'] = S1['description'].iloc[0]
                        else:
                            PDB_act.loc[index, 'entrez'] = None
                            PDB_act.loc[index, 'gene_type'] = None
                            PDB_act.loc[index, 'hgnc_symbol'] = None
                            PDB_act.loc[index, 'description'] = None  
                            Des='Failed to convert the Uniprot ID'

                    PDB_act['datasource'] = 'PDB'
                    # pChEMBL score already calculated, rename column
                    PDB_act.rename(columns={'pchembl': 'pchembl_value'}, inplace=True)
                    statement = 'completed'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Compound doesn\'t have any InChI'

        return PDB_act, statement, PDB_raw
    
    
    def compounds(self, input_protein: pd.DataFrame, pChEMBL_thres: float):
        """
        Retrieves compounds from PDB database interacting with proteins passed as input.

        Steps
        -----
        - Filters PDB database to obtain compounds interacting with input proteins
        - Extracts CID using PDB graphql API \\
        Constraints:
            - Organism = homo sapiens
        - Uses Pubchempy to obtain compounds info to return, searching with CID.

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
            - datasource (PDB)
            - pchembl value
            - notes (NaN)
        """

        columns = ['inchi','inchikey','isomeric_smiles','iupac_name','datasource','pchembl_value']
        # Create an empty DataFrame with the specified columns
        PDB_c1 = pd.DataFrame(columns=columns)
        PDB_raw = pd.DataFrame()

        input_protein = input_protein.dropna(subset=['uniprot'])
        input_protein = input_protein.drop_duplicates(subset='uniprot').reset_index(drop=True) 
        # Check if there are any input proteins remaining
        if len(input_protein) > 0:
            # Filter database
            input_protein_id = input_protein['uniprot'].iloc[0]
            PDB_raw = self.data_manager.retrieve_raw_data('Uniprot ID', input_protein_id)

            # Check if at least one match has been found
            if len(PDB_raw) > 0:

                PDB_act = self._retrieve_prot_info(PDB_raw)

                PDB_act = self._filter_database(PDB_act, pChEMBL_thres)

                if len(PDB_act) > 0:

                    pc = PubChemServer()
                    
                    if 'CID' not in PDB_act:
                        # Remove duplicates and null values for inchikey
                        PDB_act = PDB_act.dropna(subset=['ligand ID']).drop_duplicates(subset=['ligand ID'])\
                            .reset_index(drop=True)    
                        
                        PDB_act = self._retrieve_cids(PDB_act)

                    # Retrieve compounds using cids
                    compounds = self._pubchem_search_cid(PDB_act, columns, pc)

                    if len(compounds) > 0:
                        PDB_info = PDB_act[['CID', 'ligand ID', 'pchembl', 'affinity', 
                                            'unit', 'PDB code', 'organism', 'protein name']].drop_duplicates()
                        
                        PDB_c1 = pd.merge(compounds, PDB_info, on='CID')
                    

                    # Check if at least one compound has been found
                    if len(PDB_c1) > 0:
                        PDB_c1['datasource'] = 'PDB'
                        PDB_c1.rename(columns={'pchembl': 'pchembl_value'}, inplace=True)
                        PDB_c1.loc[:, 'notes'] = np.nan 
                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using PubChem'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein doesn\'t contain uniprot'  

        return PDB_c1, statement, PDB_raw


    def _retrieve_comp_id(self, inchi: str) -> str:

        # Update InChI value in query parameters
        self.api_query["query"]["parameters"]["value"] = inchi

        # Make the request
        response = requests.post(self.api_url, json=self.api_query)

        # Check for successful response
        if response.status_code == 200:
            # Save JSON response
            data = response.json()
        else:
            return None

        # Check for molecule with similarity score 1 for InChI and return PDB id, otherwise return none
        for mol in data['result_set']:
            if mol['score'] == 1:
                return mol['identifier']
            
        return None
    
    def _retrieve_cids(self, pdb_act: pd.DataFrame) -> pd.DataFrame:

        ids = list(pdb_act['ligand ID'].unique())
        # Define the GraphQL query
        query = """
        {
            chem_comps(comp_ids: %s) {
                rcsb_id
                comp: rcsb_chem_comp_related {
                    source: resource_name
                    comp_id: resource_accession_code
                }
            }
        }
        """ % str(ids).replace("'", '"')

        # Make the request
        response = requests.post(self.graphql_url, json={'query': query})

        # Check for successful response
        if response.status_code == 200:
            # Print the JSON response
            data = response.json()
            # Extract CID data of pubchem if available
            pdb_info = [
                (chem_comp['rcsb_id'],
                next(
                (int(rel['comp_id']) for rel in chem_comp['comp'] if rel['source'] == "PubChem"), None)
                )
                for chem_comp in data['data']['chem_comps']
            ]

            # Create DataFrame
            comp_info = pd.DataFrame(pdb_info, columns=["ligand ID", "CID"])

            pdb_act = pd.merge(pdb_act, comp_info, on='ligand ID', how='left')

        else:
            pdb_act['CID'] = None
        
        return pdb_act

    
    def _retrieve_prot_info(self, pdb_raw: pd.DataFrame) -> pd.DataFrame:

        ids = list(pdb_raw['PDB code'].unique())
        # Turn targets PDB code to uppercase and appens entity id (e.g. 1n4k -> 1N4K_1)
        ids = [id.upper() + "_1" for id in ids]
        # Define the GraphQL query
        query = """
        {
            polymer_entities(entity_ids: %s) {
                rcsb_id
                organism: rcsb_entity_source_organism {
                    scientific_name
                }
            }
        }
        """ % str(ids).replace("'", '"')

        # Make the request
        response = requests.post(self.graphql_url, json={'query': query})

        # Check for successful response
        if response.status_code == 200:
            # Save JSON response
            data = response.json()

            # Extract ID and data for each target
            pdb_info = [
                # Remove entity id from code, leaving only entry id lowercase (e.g. 1N4K_1 -> 1n4k)
                (entity['rcsb_id'][:-2].lower(), 
                # select Homo sapiens if present, otherwise select first organism
                next((org['scientific_name'] for org in entity['organism'] if org['scientific_name'] == "Homo sapiens"), entity['organism'][0]['scientific_name']))
                for entity in data['data']['polymer_entities']
            ]
            
            # Create DataFrame
            prot_info = pd.DataFrame(pdb_info, columns=["PDB code", "organism"])

            pdb_raw = pd.merge(pdb_raw, prot_info, on='PDB code', how='left')
        else:
            pdb_raw['organism'] = None

        return pdb_raw
        
