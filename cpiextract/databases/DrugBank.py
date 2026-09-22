'''Loading,searching,filtering and preprocessing data from DrugBank.'''

import pandas as pd
import numpy as np
import time

from ..utils.typing import Connection
from collections import defaultdict
from ..servers.PubChemServer import PubChemServer
from ..servers.BiomartServer import BiomartServer
from ..servers.MyGeneServer import MyGeneServer
from ..data_manager import *
from .Database import Database

class DB(Database):
    '''Loading,searching,filtering and preprocessing data from DrugBank.'''

    # protein_type values in DrugBank
    DEFAULT_PROTEIN_TYPES = {'target', 'enzyme', 'carrier', 'transporter'}

    # actions/action values for small molecules interactions indicating binding
    VALID_ACTIONS = {'inhibitor', 'antagonist', 'agonist', 'modulator', 'binder', 'ligand', 'activator',
        'cofactor', 'potentiator', 'substrate', 'blocker', 'partial agonist',
        'positive allosteric modulator', 'negative modulator', 'inverse agonist',
        'allosteric modulator', 'stimulator', 'stabilization', 'binding', 'metabolizer',
        'adduct', 'neutralizer', 'degradation', 'inactivator', 'chaperone', 'cleavage',
        'cross-linking/alkylation', 'inhibitory allosteric modulator', 'aggregation inhibitor',
        'positive modulator', 'inhibition of synthesis', 'partial antagonist', 'weak inhibitor',
        'oxidizer', 'reducer', 'nucleotide exchange blocker', 'translocation inhibitor',
        'intercalation', 'incorporation into and destabilization'}

    OUTPUT_COLUMNS=['entrez','gene_type','hgnc_symbol','description','pchembl_value',
                    'standard_type','datasource']

    # Final output schema for proteins()
    PROTEIN_OUTPUT_COLUMNS = ['inchi','inchikey','smiles','connectivity_smiles','iupac_name',
                'datasource','pchembl_eq','standard_type']

    PUBCHEM_FETCH_COLUMNS = ['inchi','inchikey','smiles','connectivity_smiles','iupac_name',
                'datasource','pchembl_value']
 
    @staticmethod
    def _has_valid_action(actions_str) -> bool:
        """
        True if actions_str (pipe-joined) contains at least one VALID_ACTIONS entry, or is
        blank/null entirely (see VALID_ACTIONS comment for why blank is kept, not excluded).
        """
        if pd.isna(actions_str) or not str(actions_str).strip():
            return True
        individual_actions = {a.strip() for a in str(actions_str).split('|')}
        return len(individual_actions & DB.VALID_ACTIONS) > 0

    def __init__(self, connection:Connection|None=None, database:pd.DataFrame|None=None, 
                gene_server=None, server_select='mygene'):
        
        if gene_server is not None:
            self.gene_server = gene_server
        elif server_select == 'mygene':
            self.gene_server = MyGeneServer()
        elif server_select == 'biomart':
            self.gene_server = BiomartServer()
        else:
            raise ValueError(f"server_select must be 'mygene' or 'biomart', got '{server_select}'")

        if database is not None:
            self.data_manager = LocalManager(database)
        else:
            self.data_manager = SQLManager(connection, 'DB')


    def compounds(self, input_comp: pd.DataFrame, protein_types: set|None=None) -> tuple[pd.DataFrame, str, pd.DataFrame]:

        """
        Retrieves proteins from DrugBank database interacting with compound passed as input.

        Steps
        -----
        - Finds matches on inchikey for the compound
          passed as input.
        - Finds all proteins of the selected protein_types interacting with input compound. \\
        Constraints:
            - Only Homo Sapiens interactions \\
        
        Parameters
        ----------
        DB_data : DataFrame
            Dataframe containing all DrugBank database info
        input_comp : DataFrame
            Dataframe of input compounds from which interacting proteins are found
        protein_types : set
            Which of DrugBank's protein_type categories to include - target, enzyme, carrier,
            transporter. Defaults to all four (DEFAULT_PROTEIN_TYPES) if not specified.

        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values: \\
            entrez, gene_type, hgnc_symbol, description, datasource (DB), pchembl_value
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all DrugBank info about the input compound
        """

        if protein_types is None:
            protein_types = self.DEFAULT_PROTEIN_TYPES

        columns = self.OUTPUT_COLUMNS
        # Create an empty DataFrame with the specified columns 
        DB_act = pd.DataFrame(columns=columns) 
        # Create a drugbank activity dataframe 
        DB_raw = pd.DataFrame(columns=columns) 

        input_comp = input_comp.dropna(subset=['inchikey'])
        if len(input_comp) > 0:
            input_comp_id = input_comp['inchikey_fb'][0]
            DB_raw = self.data_manager.retrieve_raw_data('FirstBlock', input_comp_id)

            # Check if at least one match has been found
            if len(DB_raw) > 0: 
                # Remove null values for drugbank_id and HGNC
                DB_act = DB_raw.dropna(subset=['drugbank_id', 'HUGO Gene Nomenclature Committee (HGNC)'])
                # Filter to Human only interactions
                DB_act = DB_act[DB_act['organism'] == 'Humans'].reset_index(drop=True) 
                # Filter to the selected protein_type(s)
                DB_act = DB_act[DB_act['protein_type'].isin(protein_types)].reset_index(drop=True)
                # Filter to genuine direct-binding actions
                if 'actions' in DB_act.columns:
                    DB_act = DB_act[DB_act['actions'].apply(self._has_valid_action)].reset_index(drop=True)

                # Check if at least one interaction has been found
                if len(DB_act) > 0:
                    # Unify gene identifiers by Harmonizing IDs
                    ensembl = self.gene_server

                    # Search by hgnc id match to biomart
                    input_type='hgnc_id' 
                    names=['hgnc_id','entrez','gene_type','hgnc_symbol','description']
                    attributes = ['hgnc_id','entrezgene_id','gene_biotype','hgnc_symbol','description']
                    
                    db_targets = pd.DataFrame()
                    
                    input_genes = list(DB_act['HUGO Gene Nomenclature Committee (HGNC)'])
                    
                    db_targets = ensembl.subset_search(input_type, input_genes, attributes, names)
                    
                    # Fix: strip 'HGNC:' prefix to match DB_act's HGNC column format
                    db_targets['hgnc_id'] = db_targets['hgnc_id'].astype(str).str.replace('HGNC:', '', regex=False)

                    # For each compound, assign specific biomart column values to the ones from the original DrugBank database
                    db_targets['hgnc_id']=db_targets['hgnc_id'].astype(str)
                    DB_act['HUGO Gene Nomenclature Committee (HGNC)']=DB_act['HUGO Gene Nomenclature Committee (HGNC)'].str.replace('HGNC:','',regex=False)
                    for index, row in DB_act.iterrows():
                        S1 = db_targets.loc[db_targets['hgnc_id']==row['HUGO Gene Nomenclature Committee (HGNC)']]
                        if len(S1) > 0:
                            DB_act.loc[index,'entrez'] = S1['entrez'].iloc[0]
                            DB_act.loc[index,'gene_type'] = S1['gene_type'].iloc[0]
                            DB_act.loc[index,'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                            DB_act.loc[index,'description'] = S1['description'].iloc[0]
                            DB_act.loc[index, 'note'] ='Harmonized gene ID'
                        else:
                            DB_act.loc[index, 'entrez'] = None
                            DB_act.loc[index, 'gene_type'] = None
                            DB_act.loc[index, 'hgnc_symbol'] = None
                            DB_act.loc[index, 'description'] = None  
                            DB_act.loc[index, 'note'] ='Failed to harmonize gene ID'
                        DB_act['datasource'] = 'DrugBank'
                        DB_act['pchembl_eq'] = np.nan
                        DB_act['standard_type'] = np.nan
                        statement = 'completed'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement= 'Compound does not have any inchikeys'        
        return DB_act, statement, DB_raw

    def proteins(self, input_protein: pd.DataFrame) -> tuple[pd.DataFrame, str ,pd.DataFrame]:
 
        """
        Retrieves compounds from DrugBank database interacting with proteins passed as input.
 
        Steps
        -----
        - Filters DrugBank database to find compounds interacting with input proteins: \\
        Constraints:
            - organism = Humans
            - compound must have inchikey
        - Uses Pubchempy to obtain compound info to return searching with inchikey.

        
        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compound are found
        DB_data: DataFrame
            Dataframe containing all drugbank database info
 
        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values: \\
            inchi, inchikey, smiles, iupac_name, datasource (db), pchembl_value, notes (activity type)
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all DrugBank info about the protein
        """
 
        columns = self.PROTEIN_OUTPUT_COLUMNS
        # Create an empty DataFrame with the specified columns
        db_c1 = pd.DataFrame(columns=columns)
        DB_raw = pd.DataFrame()
        # Drop duplicated and null values of HGNC id
        input_protein = input_protein.dropna(subset=['hgnc_id'])
        input_protein = input_protein.drop_duplicates(subset='hgnc_id').reset_index(drop=True)
 
        # Check if there are any input proteins remaining
        if len(input_protein) > 0:
            # Search compounds interacting with input gene
            input_protein_id = input_protein['hgnc_id'][0]
 
            DB_raw = self.data_manager.retrieve_raw_data('HUGO Gene Nomenclature Committee (HGNC)', input_protein_id)
            
            # Check if at least one interaction has been found
            if len(DB_raw) > 0:
 
                # Filter to Human only interactions
                DB_act = DB_raw[DB_raw['organism'] == 'Humans'].copy().reset_index(drop=True)
                # Filter to genuine direct-binding actions
                if 'actions' in DB_act.columns:
                    DB_act = DB_act[DB_act['actions'].apply(self._has_valid_action)].reset_index(drop=True)
                                 
                if len(DB_act) > 0:
 
                    pc = PubChemServer()
 
                    if 'CID' in DB_act:
                        # Retrieve compounds using cids
                        compounds = self._pubchem_search_cid(DB_act, self.PUBCHEM_FETCH_COLUMNS, pc)
                        
                        if len(compounds) > 0:
                            db_c1 = pd.concat([db_c1, compounds])
                            db_info = DB_act[['CID', 'drugbank_id', 'name', 'HUGO Gene Nomenclature Committee (HGNC)']].\
                                rename(columns={'drugbank_id': 'compound_id', 
                                                'name': 'compound_name', 'HUGO Gene Nomenclature Committee (HGNC)': 'hgnc_id'}).drop_duplicates()
                            
                            db_c1 = pd.merge(db_c1, db_info, on='CID')
                    else:
                        # Remove duplicates and null values for inchikey
                        DB_act = DB_act.dropna(subset=['inchikey']).drop_duplicates(subset=['inchikey'])\
                            .reset_index(drop=True)    
                        compounds = pd.DataFrame()
                        selected_columns = pc.get_columns(self.PUBCHEM_FETCH_COLUMNS[:-2])
                        # Perform search using Inchi key
                        for inchi in DB_act['inchikey']:
                            try:
                                # Retrieve compound using pubchempy
                                comp = (pc.get_compounds(inchi, selected_columns, namespace='inchikey'))
                                compounds = pd.concat([compounds, comp])
                            except:
                                None
                            time.sleep(0.5)
 
                        if len(compounds) > 0:
                            db_c1 = pd.concat([db_c1, compounds.rename(columns=pc.properties)])
                            db_info = DB_act[['inchikey', 'drugbank_id', 'name', 'HUGO Gene Nomenclature Committee (HGNC)']].\
                                    rename(columns={'drugbank_id': 'compound_id', 
                                                    'name': 'compound_name', 'HUGO Gene Nomenclature Committee (HGNC)': 'hgnc_id'}).drop_duplicates()
                            # Add additional values from activity dataframe
                            db_c1 = pd.merge(db_c1, db_info, on='inchikey', how='left')
 
                    # Check if at least one compound has been found
                    if len(db_c1) > 0:
                        db_c1['datasource'] = 'DrugBank'
                        db_c1['pchembl_eq'] = np.nan
                        db_c1['standard_type'] = np.nan
                        db_c1['notes'] = np.nan
                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using Pubchem'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein does not contain HGNC'                  
 
        return db_c1, statement, DB_raw