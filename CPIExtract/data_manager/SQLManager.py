from .DataManager import DataManager
import pandas as pd
import mysql

class SQLManager(DataManager):

    def retrieve_raw_data(self, filter_column, filter_value, **kwargs):

        if self.data is None:
            raise ValueError('SQL connection not available')
        
        raw_data = pd.DataFrame()
        with self.data.cursor() as cursor:
            try:
                # SQL query to select data from the table where the specified column has a certain value
                query = f"SELECT * FROM {self.database} WHERE `{filter_column}` in (%s)"
                if type(filter_value) == list:
                    values = ",".join(filter_value)
                    cursor.execute(query, (values,))
                else:
                    cursor.execute(query, (filter_value,))

                # Fetch all the rows
                rows = cursor.fetchall()

                # Create a DataFrame from the result
                columns = [desc[0] for desc in cursor.description]
                raw_data = pd.DataFrame(rows, columns=columns)

            # Connection error
            except mysql.connector.Error as err:
                print(f'Error while fetching data from SQL server: {err}')
            
            return raw_data