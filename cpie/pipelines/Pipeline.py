from abc import ABC
from ..databases import *
from ..sql_connection import connect_to_mysql
import pandas as pd
import numpy as np

class Pipeline(ABC):

    def __init__(self, execution_mode='local', dbs: dict=None, server_info: dict=None) -> None:
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
            Dictionary containing the configuration info to connect to the SQL server \\
            Example:
            {
                "host": "127.0.0.1",
                "user": "user",
                "password": "password",
                "database": "cpie",
            }
        """

        self.cnx = None

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
            
        self.databases = {
            'pc': PubChem(),
            'chembl': ChEMBL(database=dbs.get('chembl', None), connection=self.cnx),
            'bdb': BindingDB(database=dbs.get('bdb', None), connection=self.cnx),
            'stitch': Stitch(database=dbs.get('stitch', None), connection=self.cnx),
            'ctd': CTD(database=dbs.get('ctd', None), connection=self.cnx),
            'dtc': DTC(database=dbs.get('dtc', None), connection=self.cnx),
            'otp': OTP(),
            'dc': DrugCentral(database=dbs.get('dc', None), connection=self.cnx),
            'db': DB(database=dbs.get('db', None), connection=self.cnx)
        }

        self.database_args = {}

        self.sources = ['PubChem', 'ChEMBL', 'BindingDB', 'Stitch', 'CTD', 'DTC', 'OTP', 'DrugCentral', 'DrugBank']


    def _aggregate_pchembl(self, data: pd.DataFrame, index, comp: pd.DataFrame):
        # Calculate the pchembl value for the target, if none, states so
        pchembl = pd.to_numeric(comp['pchembl_value'], errors='coerce').dropna().unique()
        if len(pchembl) > 0:
            data.loc[index,'pchembl_count'] = len(pchembl)
            data.loc[index,'ave_pchembl'] = np.mean(pchembl)
            data.loc[index,'std_pchembl'] = np.std(pchembl)
        else:
            data.loc[index,'pchembl_count'] = 0
            data.loc[index,'ave_pchembl'] = 'Sources do not provide activity data'
            data.loc[index,'std_pchembl'] = 'Sources do not provide activity data'

        return data

    
