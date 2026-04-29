'''Loading,searching,filtering and preprocessing data from CTD.'''

import pandas as pd
import numpy as np

from ..utils.typing import Connection
from ..servers.PubChemServer import PubChemServer
from ..servers.BiomartServer import BiomartServer
from ..servers.MyGeneServer import MyGeneServer
from .Database import Database
from ..data_manager import *

class CTD(Database):

    def __init__(self, connection: Connection|None=None, database: pd.DataFrame|None=None, 
                 merge_stereoisomers=False, gene_server=None, server_select='mygene'):
        super().__init__(merge_stereoisomers)

        if gene_server is not None:
            self.gene_server = gene_server
        elif server_select == 'mygene':
            self.gene_server = MyGeneServer()
        elif server_select == 'biomart':
            self.gene_server = BiomartServer()
        else:
            raise ValueError(f"server_select must be 'mygene' or 'biomart', got '{server_select}'")

        # if not connection and not database:
        #     raise ValueError('Either SQL connection or database should be not None')
        if database is not None:
            self.data_manager = LocalManager(database)
        else:
            self.data_manager = SQLManager(connection, 'CTD')

    def _filter_database(self, CTD_raw: pd.DataFrame) -> pd.DataFrame:
                        # Filter to Homo Sapiens only
        CTD_act = CTD_raw.loc[(CTD_raw['OrganismID'] == 9606) &
                            # Only select ligand-protein binding
                            (CTD_raw['InteractionActions'].str.contains(r'affects\^binding')) &
                            # Only protein interactions are collected
                            (CTD_raw['GeneForms']==('protein'))].drop_duplicates().reset_index(drop=True)
        return CTD_act

    def interactions(self, input_comp: pd.DataFrame, merge_stereoisomers: bool=False) -> tuple[pd.DataFrame, str, pd.DataFrame]:
        """
        Retrieves proteins from CTD database interacting with compound passed as input.

        Steps
        -----
        - Finds matches for all synonyms of the compound passed as input
        - Finds all proteins interacting with input compound \\
        Constraints:
            - Only Homo Sapiens interactions
            - Only ligand-protein binding
            - Only protein interactions
        
        Parameters
        ----------
        CTD_data : DataFrame
            Dataframe containing all CTD database info
        input_comp : DataFrame
            Dataframe of input compounds from which interacting proteins are found
        merge_stereoisomers : bool
            determines if results respect stereochemical specificity of input compound

        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values: \\
            entrez, gene_type, hgnc_symbol, description, datasource (CTD), pchembl_value
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all CTD info about the input compound
        """

        columns = ['entrez','gene_type','hgnc_symbol','description','pchembl_value','datasource']
        # Create an empty DataFrame with the specified columns
        CTD_act = pd.DataFrame(columns=columns)
        CTD_raw = pd.DataFrame()
        input_comp = input_comp.dropna(subset=['inchikey'])
        if len(input_comp) > 0:
            if merge_stereoisomers == True: #FirstBlock only
                input_comp_id = input_comp['inchikey_fb'][0]
                CTD_raw = self.data_manager.retrieve_raw_data('FirstBlock', input_comp_id)
            else: #Full inchikey
                input_comp_id = input_comp['inchikey'][0]
                CTD_raw = self.data_manager.retrieve_raw_data('inchikey', input_comp_id)

            # Check if at least one match has been found
            if len(CTD_raw) > 0:
                # Ensure a binding type interaction occurs
                CTD_act = self._filter_database(CTD_raw)
                
                # Check if at least one interaction has been found
                if len(CTD_act) > 0:
                    # Unify gene identifiers by Harmonizing IDs
                    ensembl = self.gene_server
                    
                    attributes = ['entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                    names = ['entrez','gene_type','hgnc_symbol','description']
                    # Use entrez for conversion
                    input_type='entrezgene_id' 
                    
                    ctd_targets=pd.DataFrame(columns=names)
                    
                    input_genes=list(CTD_act['GeneID'])
                    
                    ctd_targets=ensembl.subset_search(input_type, input_genes, attributes, names)

                    # For each compound, assign protein query column values to the ones from the original database
                    ctd_targets['entrez']=ctd_targets['entrez'].astype(str)
                    CTD_act['GeneID']=CTD_act['GeneID'].astype(str)
                    for index, row in CTD_act.iterrows():
                        S1 = ctd_targets.loc[ctd_targets['entrez']==row['GeneID']]
                        if len(S1) > 0:
                            CTD_act.loc[index,'entrez'] = S1['entrez'].iloc[0]
                            CTD_act.loc[index,'gene_type'] = S1['gene_type'].iloc[0]
                            CTD_act.loc[index,'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                            CTD_act.loc[index,'description'] = S1['description'].iloc[0]
                            CTD_act.loc[index, 'note'] ='Harmonized gene ID'
                        else:
                            CTD_act.loc[index, 'entrez'] = None
                            CTD_act.loc[index, 'gene_type'] = None
                            CTD_act.loc[index, 'hgnc_symbol'] = None
                            CTD_act.loc[index, 'description'] = None
                            CTD_act.loc[index, 'note'] ='Failed to harmonize gene ID'
                    CTD_act['datasource'] = 'CTD'
                    # No activity values in database so no pChEMBL score calculated
                    CTD_act['pchembl_value'] = np.nan
                    statement = 'completed'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Compound does not have any cid'

        return CTD_act, statement, CTD_raw
    
    
    def compounds(self, input_protein: pd.DataFrame, merge_stereoisomers: bool=False) -> tuple[pd.DataFrame, str, pd.DataFrame]:
        """
        Retrieves compounds from CTD database interacting with proteins passed as input.

        Steps
        -----
        - Filters CTD database to obtain compounds interacting with input proteins \\
        Constraints:
            - Organism = homo sapiens
            - InteractionActions contains 'affects\\^binding'
            - GeneForms = protein 
        - Uses Pubchempy to obtain compounds info to return, searching with CID.

        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compound are found
        CTD_data : DataFrame
            Dataframe containing all CTD database info

        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values:
            - inchi
            - inchikey
            - isomeric smiles
            - iupac name
            - datasource (CTD)
            - pchembl value (NaN)
            - notes (NaN)
        """

        columns = ['inchi','inchikey','isomeric_smiles','iupac_name','datasource','pchembl_value']
        # Create an empty DataFrame with the specified columns
        ctd_c1 = pd.DataFrame(columns=columns)
        ctd_raw = pd.DataFrame()
        input_protein = input_protein.dropna(subset=['entrez'])
        input_protein = input_protein.drop_duplicates(subset='entrez').reset_index(drop=True) 
        # Check if there are any input proteins remaining
        if len(input_protein) > 0:
            # Filter database
            input_protein_id = int(input_protein['entrez'].iloc[0])
            ctd_raw = self.data_manager.retrieve_raw_data('GeneID', input_protein_id)

            # Check if at least one match has been found
            if len(ctd_raw) > 0:

                ctd_act = self._filter_database(ctd_raw)

                if len(ctd_act) > 0:
                    pc = PubChemServer()
                    # Retrieve compounds using cids
                    ctd_c1 = self._pubchem_search_cid(ctd_raw, columns, pc)
                    
                    if len(ctd_c1) > 0:
                        ctd_c1['datasource'] = 'CTD'
                        ctd_c1['pchembl_value'] = np.nan
                        ctd_c1.loc[:, 'notes'] = np.nan 
                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using PubChem'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein doesn\'t contain entrez'  

        return ctd_c1, statement, ctd_raw   
