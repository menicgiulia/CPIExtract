from abc import ABC, abstractmethod

class DataManager(ABC):

    def __init__(self, data_source, db_name=None):
        self.data = data_source
        self.database = db_name

    @abstractmethod
    def retrieve_raw_data(self, filter_column, filter_value, **kwargs):
        raise NotImplementedError
