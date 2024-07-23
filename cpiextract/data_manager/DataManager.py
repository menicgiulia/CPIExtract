'''Base data manager for data retrieving and filtering.'''

from abc import ABC, abstractmethod
from ..utils.typing import CIDs, Generic, T


class DataManager(ABC, Generic[T]):
    '''Base data manager for data retrieving and filtering.'''

    def __init__(self, data_source: T, db_name: str| None = None):
        self.data = data_source
        self.database = db_name

    @abstractmethod
    def retrieve_raw_data(self, filter_column: str, filter_value: CIDs, **kwargs):
        '''Retrieve raw data and filter.'''
        raise NotImplementedError
