from .DataManager import DataManager

class APIManager(DataManager):

    def retrieve_raw_data(self, filter_column, filter_value, **kwargs):
        
        # Call API function specified by filter_column
        raw_data = self.data[filter_column](filter_value)
        
        return raw_data