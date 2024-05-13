from .DataManager import DataManager
import pandas as pd

class LocalManager(DataManager):

    def retrieve_raw_data(self, filter_column, filter_value, **kwargs):

        if self.data is None:
            raise ValueError('Database not available')
        
        if type(filter_value) == list:
            # Filter data, create copy and reset index
            raw_data = self.data[self.data[filter_column].isin(filter_value)].copy().reset_index(drop=True)
        else:
            # Filter data, create copy and reset index
            raw_data = self.data[self.data[filter_column] == filter_value].copy().reset_index(drop=True)
        
        return raw_data