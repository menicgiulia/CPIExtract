'''Data managers for mySQL tables' data retrieving and filtering,performed by mySQL.'''

import mysql.connector
from ..utils.typing import CIDs, Connection
from .DataManager import DataManager
import pandas as pd
import mysql

class SQLManager(DataManager[Connection]):
    '''Data managers for data retrieving and filtering,performed by mySql.'''

    def retrieve_raw_data(self, filter_column: str, filter_value: CIDs, case_insensitive: bool = False, **kwargs) -> pd.DataFrame:
        '''Retrieve raw data from SQL tables.'''
        if self.data is None:
            raise ValueError('SQL connection not available')

        raw_data = pd.DataFrame()
        column_expr = f"LOWER(`{filter_column}`)" if case_insensitive else f"`{filter_column}`"

        with self.data.cursor() as cursor:
            try:
                if isinstance(filter_value, list):
                    values = [v.lower() if case_insensitive and isinstance(v, str) else v for v in filter_value]
                    placeholders = ",".join(["%s"] * len(values))
                    query = f"SELECT * FROM {self.database} WHERE {column_expr} in ({placeholders})"
                    cursor.execute(query, tuple(values))
                else:
                    value = filter_value.lower() if case_insensitive and isinstance(filter_value, str) else filter_value
                    query = f"SELECT * FROM {self.database} WHERE {column_expr} in (%s)"
                    cursor.execute(query, (value,))

                rows = cursor.fetchall()
                columns = [desc[0] for desc in cursor.description]
                raw_data = pd.DataFrame(rows, columns=columns)
            except mysql.connector.Error as err:
                print(f'Error while fetching data from SQL server: {err}')

            return raw_data