'''Call API function specified by filter_column.'''

from pandas import DataFrame
from .DataManager import DataManager
from ..utils.typing import Callable, Strs

class APIManager(DataManager[dict[str, Callable[..., DataFrame]]]):
    '''Call API function specified by filter_column.'''
    def retrieve_raw_data(self, filter_column: str, filter_value: Strs, **kwargs) -> DataFrame:
        '''Call API function specified by filter_column.'''
        raw_data = self.data[filter_column](filter_value)
        
        return raw_data