'''Data managers for csv/tsv data retrieving and filtering,performed by pandas.'''

from .DataManager import DataManager
import pandas as pd
from ..utils.typing import Strs

class LocalManager(DataManager[pd.DataFrame]):
    '''Data managers for csv/tsv data retrieving and filtering,performed by pandas.'''

    def retrieve_raw_data(self, filter_column: str, filter_value: Strs, **kwargs) -> pd.DataFrame:
        '''Retrieve raw data from local databases.'''
        if self.data is None:
            raise ValueError('Database not available')
        
        if isinstance(filter_value, list):
            # Filter data, create copy and reset index
            raw_data = self.data[self.data[filter_column].isin(filter_value)].copy().reset_index(drop=True)
        else:
            # Filter data, create copy and reset index
            raw_data = self.data[self.data[filter_column] == filter_value].copy().reset_index(drop=True)
        
        return raw_data