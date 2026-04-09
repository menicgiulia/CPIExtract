'''The server to perform pubchem requests.'''

from ..utils.typing import Strs
import pubchempy as pcp
import json
# import time
from typing import Union
import pandas as pd
from ..utils.helper import Singleton

class PubChemServer(metaclass=Singleton):
    '''The server to perform pubchem requests.'''

    def __init__(self) -> None:
        self.properties = {value: key for key, value in pcp.PROPERTY_MAP.items()}

    def get_columns(self, selected_columns: list[str]) -> list[str]:
        '''Return selected columns' original field names in pubchem.'''
        return [pcp.PROPERTY_MAP[col] for col in selected_columns]

    def get_compounds(self, comp: Strs, selected_columns: list[str], domain: str='compound', namespace: str='cid') -> pd.DataFrame:
        '''Return compounds' selected properties from pubchem.
        
        Parameters
        ----------
        comp: list[str] | str
            A (set of) compound id(s).
        selected_columns: list[str]
            Specified properties of the compound(s).
        domain: str
            Search structures in `compound` domain by default.
        namespace: str
            ID namespace of `comp`.

        Returns
        -------
        DataFrame
            A table contains selected properties of specified compound(s).
        '''
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
    
    def get_synonyms(self, comp: Strs, domain: str='compound', namespace: str='cid') -> pd.DataFrame:
        '''Return synonyms of specified compound(s).
        
        Parameters
        ----------
        comp: list[str] | str
            A (set of) compound id(s).
        domain: str
            Search structures in `compound` domain by default.
        namespace: str
            ID namespace of `comp`.

        Returns
        -------
        DataFrame
            A table contains synonyms of specified compound(s).
        '''
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
