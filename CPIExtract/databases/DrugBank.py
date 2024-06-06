import pandas as pd
import numpy as np
import time
from collections import defaultdict
from ..servers.BiomartServer import BiomartServer
from ..servers.PubchemServer import PubChemServer
from ..data_manager import *
from .Database import Database

class DB(Database):

    def __init__(self, connection=None, database=None):
        # if not connection and not database:
        #     raise ValueError('Either SQL connection or database should be not None')
        if database is not None:
            self.data_manager = LocalManager(database)
        else:
            self.data_manager = SQLManager(connection, 'DB')


    def interactions(self, input_comp: pd.DataFrame):
        """
        Retrieves proteins from DrugBank database interacting with compound passed as input.

        Steps
        -----
        - Finds matches for all synonyms of the compound passed as input.
        - Finds all proteins (targets, enzymes, carriers, transporters) interacting with input compound. \\
        Constraints:
            - Only Homo Sapiens interactions \\
        - Uses Biomart to obtain proteins info (modified with data from original DB database) to return,
        searching with GeneID.
        
        Parameters
        ----------
        DB_data : DataFrame
            Dataframe containing all DrugBank database info
            
        input_comp : DataFrame
            Dataframe of input compounds from which interacting proteins are found

        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values:
            - entrez
            - gene_type
            - hgnc_symbol
            - description
            - datasource (DrugBank)
            - pchembl_value (nan)

        String
            A statement string describing the outcome of the database search

        DataFrame
            Raw Dataframe containing all DrugBank info about the input compound
        """

        columns = ['entrez','gene_type','hgnc_symbol','description','pchembl_value','datasource']
        # Create an empty DataFrame with the specified columns 
        DB_act = pd.DataFrame(columns=columns) 
        # Create a drugbank activity dataframe 
        DB_raw = pd.DataFrame(columns=columns) 

        input_comp = input_comp.dropna(subset=['cid'])
        if len(input_comp) > 0:

            input_comp_id = input_comp['cid'][0]
            DB_raw = self.data_manager.retrieve_raw_data('CID', input_comp_id)

            # Try searching with synonyms instead
            if len(DB_raw) == 0:
                input_comp = input_comp.dropna(subset=['synonyms'])
                # Extract all synonyms into a single list
                all_synonyms = list(input_comp['synonyms'].explode().dropna().str.lower().unique())
                # Filter DB_data based on all unique synonyms
                DB_raw = self.data_manager.retrieve_raw_data('name', all_synonyms)

            # Check if at least one match has been found
            if len(DB_raw) > 0: 
 
                # Remove null values for drugbank-id and HGNC
                DB_act = DB_raw.dropna(subset=['drugbank-id', 'HGNC'])
                # Filter to Human only interactions
                DB_act = DB_act[DB_act['organism'] == 'Humans'].reset_index(drop=True) 

                # Check if at least one interaction has been found
                if len(DB_act) > 0:
                    # Unify gene identifiers by Converting DrugBank with Biomart
                    ensembl = BiomartServer()
                    # Search by hgnc id match to biomart
                    input_type='hgnc_id' 
                    names=['hgnc_id','entrez','gene_type','hgnc_symbol','description']
                    attributes = ['hgnc_id','entrezgene_id','gene_biotype','hgnc_symbol','description']
                    
                    db_targets = pd.DataFrame()
                    
                    input_genes = list(DB_act['HGNC'])
                    
                    db_targets = ensembl.subset_search(input_type, input_genes, attributes, names)
                        
                    # For each compound, assign specific biomart column values to the ones from the original DrugBank database
                    for index, row in DB_act.iterrows():
                        S1 = db_targets.loc[db_targets['hgnc_id']==row['HGNC']]
                        if len(S1) > 0:
                            DB_act.loc[index,'entrez'] = S1['entrez'].iloc[0]
                            DB_act.loc[index,'gene_type'] = S1['gene_type'].iloc[0]
                            DB_act.loc[index,'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                            DB_act.loc[index,'description'] = S1['description'].iloc[0]
                        else:
                            DB_act.loc[index, 'entrez'] = None
                            DB_act.loc[index, 'gene_type'] = None
                            DB_act.loc[index, 'hgnc_symbol'] = None
                            DB_act.loc[index, 'description'] = None  

                        DB_act['datasource'] = 'DrugBank'
                        DB_act['pchembl_value'] = np.nan
                        statement = 'completed'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement= 'Compound doesn\'t have any cids'        
        return DB_act, statement, DB_raw


    def compounds(self, input_protein: pd.DataFrame):
        """
        Retrieves compounds from db database interacting with proteins passed as input.

        Steps
        -----
        - Filters stitch database to find compounds interacting with input proteins: \\
        Constraints:
            - organism = Humans
            - compound must have inchikey
        - Uses Pubchempy to obtain compound info (modified with data from original bdb database) to return, 
        searching with inchikey.
        
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
            inchi, inchikey, isomeric_smiles, iupac_name, datasource (db), pchembl_value (nan), notes (nan)
        """

        columns = ['inchi','inchikey','isomeric_smiles','iupac_name','datasource','pchembl_value']
        # Create an empty DataFrame with the specified columns
        db_c1 = pd.DataFrame(columns=columns)
        DB_raw = pd.DataFrame()
        # Drop duplicated and null values of HGNC id
        input_protein = input_protein.dropna(subset=['hgnc_id'])
        input_protein = input_protein.drop_duplicates(subset='hgnc_id').reset_index(drop=True)

        # Check if there are any input proteins remaining
        if len(input_protein) > 0:
            # compound_data = defaultdict(list)
            # Search compounds interacting with input gene
            input_protein_id = input_protein['hgnc_id'][0]

            DB_raw = self.data_manager.retrieve_raw_data('HGNC', input_protein_id)
            
            # Check if at least one interaction has been found
            if len(DB_raw) > 0:

                # Filter to Human only interactions
                DB_act = DB_raw[DB_raw['organism'] == 'Humans'].copy().reset_index(drop=True)                    
                                 
                if len(DB_act) > 0:

                    pc = PubChemServer()

                    if 'CID' in DB_act:
                        # Retrieve compounds using cids
                        compounds = self._pubchem_search_cid(DB_act, columns, pc)
                        
                        if len(compounds) > 0:
                            db_c1 = pd.concat([db_c1, compounds])
                            db_info = DB_act[['CID', 'drugbank-id', 'name', 'HGNC']].\
                                rename(columns={'drugbank-id': 'compound_id', 
                                                'name': 'compound_name', 'HGNC': 'hgnc_id'}).drop_duplicates()
                            
                            db_c1 = pd.merge(db_c1, db_info, on='CID')
                    else:
                        # Remove duplicates and null values for inchikey
                        DB_act = DB_act.dropna(subset=['InChIKey']).drop_duplicates(subset=['InChIKey'])\
                            .reset_index(drop=True)    
                        compounds = pd.DataFrame()
                        selected_columns = pc.get_columns(columns[:-2])
                        # Perform search using Inchi key
                        for inchi in DB_act['InChIKey']:
                            try:
                                # Retrieve compound using pubchempy
                                comp = (pc.get_compounds(inchi, selected_columns, namespace='inchikey'))
                                compounds = pd.concat([compounds, comp])
                            except:
                                None
                            time.sleep(0.5)

                        if len(compounds) > 0:
                            db_c1 = pd.concat([db_c1, compounds.rename(columns=pc.properties)])
                            db_info = DB_act[['InChIKey', 'drugbank-id', 'name', 'HGNC']].\
                                    rename(columns={'InChIKey': 'inchikey', 'drugbank-id': 'compound_id', 
                                                    'name': 'compound_name', 'HGNC': 'hgnc_id'}).drop_duplicates()
                            # Add additional values from activity dataframe
                            db_c1 = pd.merge(db_c1, db_info, on='inchikey', how='left')

                    # Check if at least one compound has been found
                    if len(db_c1) > 0:
                        db_c1['datasource'] = 'DrugBank'
                        db_c1['pchembl_value'] = np.nan
                        db_c1['notes'] = np.nan
                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using Pubchem'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein doesn\'t contain HGNC'                  

        return db_c1, statement, DB_raw


    # def _find_compound_info(self, DB_act):
    #     """
    #     Extract compound info (inchikey) from DrugBank database and stores it into a new column.
    
    #     Parameters
    #     ----------
    #     DB_act : DataFrame
    #         DrugBank database

    #     Returns
    #     -------
    #     DataFrame
    #         Drugbank database with a column containing compound inchikey 
    #     """
    #     for index, row in DB_act.iterrows():
    #         if type(row['calculated-properties']) == dict:
    #             for prop in row['calculated-properties'].values():
    #                 try:
    #                     if prop['kind'][0][0] == 'InChIKey':
    #                         DB_act.loc[index, 'InChIKey'] = prop['value'][0][0]
    #                         break
    #                 except:
    #                     None
    #     return DB_act


    # def preprocess_db(self, DB_data):

    #     DB_data['name'] = DB_data['name'].apply(lambda x: x[0][0].lower() if x else None)
    #     DB_data['drugbank-id'] = DB_data['drugbank-id'].apply(lambda x: x[0][0] if x else None)

    #     # Add inchikey column
    #     DB_data['InChIKey'] = None
    #     # Extract inchikey
    #     DB_data = self._find_compound_info(DB_data)

    #     protein_data = defaultdict(list)
    #     for _, row in DB_data.iterrows():
    #         hgncs = []
    #         organisms = []
    #         protein_names = []
    #         protein_types = []
    #         for prot_type in ['targets', 'carriers', 'enzymes', 'transporters']:
    #             if row[prot_type]:
    #                 for prot in row[prot_type].values():
    #                     try:
    #                         HGNC = prot['polypeptide']['external-identifiers']['external-identifier']['identifier'][0][0]
    #                         try:
    #                             organism = prot['organism'][0][0]
    #                         except:
    #                             organism = None
    #                         protein_name = prot['name'][0][0]
    #                     except (KeyError, TypeError):
    #                         HGNC = None
    #                         organism = None
    #                         protein_name = None
    #                     hgncs.append(HGNC)
    #                     organisms.append(organism)
    #                     protein_names.append(protein_name)
    #                     protein_types.append(prot_type)
    #         protein_data['HGNCs'].append(hgncs)
    #         protein_data['organisms'].append(organism)
    #         protein_data['protein_names'].append(protein_names)
    #         protein_data['types'].append(protein_types)

    #     new_data = pd.DataFrame({'HGNC': protein_data['HGNCs'], 
    #                                     'organism': protein_data['organisms'], 
    #                                     'protein_name': protein_data['protein_names'],
    #                                     'protein_type': protein_data['types']})

    #     DB_data = pd.concat([DB_data, new_data], axis=1)
    #     DB_data = DB_data.explode(['HGNC', 'protein_name', 'protein_type']).reset_index(drop=True)
    #     DB_data = DB_data[['drugbank-id', 'name', 'InChIKey', 'HGNC', 'organism', 'protein_name', 'protein_type']]
    #     return DB_data
