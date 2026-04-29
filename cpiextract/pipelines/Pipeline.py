'''Template for a pipeline.Load,search,filter,postprocess datas from all of the databases.'''

import pandas as pd
import numpy as np
import os

from abc import ABC
from ..databases import *
from ..sql_connection import connect_to_mysql
from ..servers.BiomartServer import BiomartServer
from ..servers.MyGeneServer import MyGeneServer

class Pipeline(ABC):
    '''Template for a pipeline.Load,search,filter,postprocess datas from all of the databases.'''

    def __init__(self, execution_mode: str='local', dbs: dict|None=None, server_info: dict|None=None, 
                 pubchem_files: dict|None=None, merge_stereoisomers: bool=False, server_select: str='mygene') -> None:
        """
        Parameters
        ----------
        execution_mode : {'local', 'online', 'server'}, default 'local'
            The execution modality of the pipeline:
            - local: retrieves most databases info from local csv files, if available
            - server: retrieves most databases info from MySQL server

        dbs : dict
            Dictionary containing all the database stored into DataFrame objects

        server_info : dict
            Dictionary containing the configuration info to connect to the SQL server

        merge_stereoisomers : bool, default False
            Whether to match compounds using only the first block of InChIKey.
            - False: Exact InChIKey match (single stereoisomer)
            - True: First-block match (all stereoisomers)
        
        pubchem_files : dict, optional
            Dictionary containing path to PubChem local database for faster processing.
            Expects local in a dubckdb structure: pubchem.duckdb
            If database not found, PubChem will use API calls.

        server_select : {'mygene', 'biomart'}, default 'mygene'
            The service used for gene identifier harmonization.
            - mygene: queries MyGene.info (recommended, more stable)
            - biomart: queries Ensembl BioMart
        """

        self.cnx = None
        self.merge_stereoisomers = merge_stereoisomers

        if execution_mode not in ['local', 'server']:
            raise ValueError("Invalid execution mode provided. Must be one of the following: \
                             'local', 'server'.")
        
        if execution_mode in ['local']:
            if dbs is None:
                raise ValueError("Databases not specified.")
        
        if execution_mode in ['server']: 
            if server_info is None:
                raise ValueError("Server info not specified.")
            
            dbs = {}
            self.cnx = connect_to_mysql(config=server_info)

            if not self.cnx or not self.cnx.is_connected():
                raise ConnectionError("Couldn't connect to SQL server.")
        
        # Get PubChem database path with automatic detection
        pc_bioact = self._get_pubchem_path(pubchem_files)
        
        if server_select == 'mygene':
            self.gene_server = MyGeneServer()
        elif server_select == 'biomart':
            self.gene_server = BiomartServer()
        else:
            raise ValueError(f"server_select must be 'mygene' or 'biomart', got '{server_select}'")

        self.databases = {
            'pc': PubChem(merge_stereoisomers=self.merge_stereoisomers,bioact_file=pc_bioact,gene_server=self.gene_server),
            'chembl': ChEMBL(database=dbs.get('chembl', None), connection=self.cnx,merge_stereoisomers=self.merge_stereoisomers,gene_server=self.gene_server),
            'bdb': BindingDB(database=dbs.get('bdb', None), connection=self.cnx,merge_stereoisomers=self.merge_stereoisomers,gene_server=self.gene_server),
            'stitch': Stitch(database=dbs.get('stitch', None), connection=self.cnx,merge_stereoisomers=self.merge_stereoisomers,gene_server=self.gene_server),
            'ctd': CTD(database=dbs.get('ctd', None), connection=self.cnx,merge_stereoisomers=self.merge_stereoisomers,gene_server=self.gene_server),
            'dtc': DTC(database=dbs.get('dtc', None), connection=self.cnx,merge_stereoisomers=self.merge_stereoisomers,gene_server=self.gene_server),
            'otp': OTP(database=dbs.get('otp', None),chembl_database=dbs.get('chembl', None), merge_stereoisomers=self.merge_stereoisomers,gene_server=self.gene_server),
            'dc': DrugCentral(database=dbs.get('dc', None), connection=self.cnx,merge_stereoisomers=self.merge_stereoisomers,gene_server=self.gene_server),
            'db': DB(database=dbs.get('db', None), connection=self.cnx,merge_stereoisomers=self.merge_stereoisomers,gene_server=self.gene_server)
        }

        self.database_args = {}

        self.sources = ['PubChem', 'ChEMBL', 'BindingDB', 'Stitch', 'CTD', 'DTC', 'OTP', 'DrugCentral', 'DrugBank']

    def _get_pubchem_path(self, pubchem_files):
        """
        Get PubChem database path with automatic detection.
        """
        if pubchem_files:
            # Check for direct database file specification first
            if 'db_file' in pubchem_files:
                return pubchem_files['db_file']
    
        # Try environment variable
        pubchem_dir = os.environ.get('PUBCHEM_DATA')
    
        # Fallback to default location
        if not pubchem_dir:
            pubchem_dir = os.path.join(os.path.expanduser('~'), 'CPE_data', 'pubchem')
        if not os.path.exists(pubchem_dir):
            return None
    
        # Return path to database
        db_file = os.path.join(pubchem_dir, 'pubchem.duckdb')
        if os.path.exists(db_file):
            return db_file
    
        return None

    def _aggregate_pchembl(self, data: pd.DataFrame, index: int, comp: pd.DataFrame) -> pd.DataFrame:
        '''Calculate mean and std. of pchembl value for the target.'''
        # Calculate the pchembl value for the target, if none, states so
        pchembl = pd.to_numeric(comp['pchembl_value'], errors='coerce').dropna().unique()
        if len(pchembl) > 0:
            data.loc[index,'pchembl_count'] = len(pchembl)
            data.loc[index,'ave_pchembl'] = np.mean(pchembl)
            # Use nanstd to handle NaNs
            with np.errstate(invalid='ignore'):
                data.loc[index,'std_pchembl'] = np.std(pchembl)
        else:
            data.loc[index,'pchembl_count'] = 0
            data = data.astype({'ave_pchembl': 'object', 'std_pchembl': 'object'})
            data.loc[index,'ave_pchembl'] = 'Sources do not provide activity data'
            data.loc[index,'std_pchembl'] = 'Sources do not provide activity data'

        return data
