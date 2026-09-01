'''Data managers for csv/tsv data retrieving and filtering,performed by pandas.'''

from .DataManager import DataManager
import pandas as pd
from ..utils.typing import CIDs

class LocalManager(DataManager[pd.DataFrame]):
    '''Data managers for csv/tsv data retrieving and filtering,performed by pandas.'''

    def retrieve_raw_data(self, filter_column: str, filter_value: CIDs, case_insensitive: bool = False, **kwargs) -> pd.DataFrame:
        '''Retrieve raw data from local databases.'''
        if self.data is None:
            raise ValueError('Database not available')

        column = self.data[filter_column]
        value = filter_value

        if case_insensitive:
            column = column.str.lower()
            if isinstance(value, list):
                value = [v.lower() if isinstance(v, str) else v for v in value]
            elif isinstance(value, str):
                value = value.lower()

        if isinstance(value, list):
            mask = column.isin(value)
        else:
            mask = column == value

        raw_data = self.data[mask].copy().reset_index(drop=True)
        return raw_data