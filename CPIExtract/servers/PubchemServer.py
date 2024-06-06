import pubchempy as pcp
import json
# import time
from typing import Union
import pandas as pd
from ..utils.helper import Singleton

class PubChemServer(metaclass=Singleton):

    def __init__(self) -> None:
        self.properties = {value: key for key, value in pcp.PROPERTY_MAP.items()}

    def get_columns(self, selected_columns: list[str]):
        return [pcp.PROPERTY_MAP[col] for col in selected_columns]

    def get_compounds(self, comp: Union[list[str], str], selected_columns: list[str], domain: str='compound', namespace: str='cid'):
        
        operation = 'property/' + ','.join(selected_columns)

        try:
            result = pcp.request(
                identifier=comp,
                domain=domain,
                namespace=namespace,
                operation=operation,
                output='JSON',
                # searchtype=None
            ).read()
        except pcp.PubChemHTTPError as e:
            raise e
        
        result = json.loads(result.decode())

        result = result['PropertyTable']['Properties']

        compounds = pd.DataFrame.from_dict(result)

        return compounds
    
    def get_synonyms(self, comp: Union[list[str], str], domain: str='compound', namespace: str='cid'):

        try:
            result = pcp.request(
                identifier=comp,
                domain=domain,
                namespace=namespace,
                operation='synonyms',
                output='JSON',
                # searchtype=None
            ).read()
        except pcp.PubChemHTTPError as e:
            raise e
        
        result = json.loads(result.decode())

        result = result['InformationList']['Information']

        compounds = pd.DataFrame.from_dict(result)
        compounds['Synonym'] = compounds['Synonym'].apply(lambda x: ','.join(x) if isinstance(x, list) else x)

        return compounds
