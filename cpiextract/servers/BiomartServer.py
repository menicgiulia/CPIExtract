'''The server to perform Biomart requests.'''

import biomart
import pandas as pd
import io
from ..utils.helper import *

class BiomartServer(metaclass=Singleton):
    '''The server to perform Biomart requests.'''

    def __init__(self) -> None:
        self.mirrors = ['useast', 'www', 'asia']
        self._connect_to_biomart()
    
    def _connect_to_biomart(self) -> None:
        for mirror in self.mirrors:
            try:
                self.server = biomart.BiomartServer(f"http://{mirror}.ensembl.org/biomart")
                self.dataset = self.server.datasets['hsapiens_gene_ensembl'] # Homo Sapiens
                return 
            except:
                continue
        raise ConnectionError('Couldn\'t connect to biomart')

    def search(self, input_type: str, input_id: list[int | str] | str | int, attributes: list[str], column_names: list[str]) -> pd.DataFrame:
        '''Search and filter input ids on Biomart.'''
        try:
            bm = self.dataset.search({
            'filters': {
                input_type: input_id
            },
            'attributes': attributes
            })

            targets = pd.read_csv(io.StringIO(bm.text), sep='\t', header=None, names=column_names)
        except:
            targets = pd.DataFrame(columns=column_names)

        return targets
    

    def subset_search(self, input_type: str, input_ids: list[int | str], attributes: list[str], column_names: list[str]) -> pd.DataFrame:
        '''Search and filter input ids on Biomart by batch mode.'''
        targets = pd.DataFrame(columns=column_names)
        if len(input_ids) > 250:
            for start, end in generate_subsets(len(input_ids), 250):
                subset = input_ids[start:end]
                sub_targets = self.search(input_type, subset, attributes, column_names)
                targets = pd.concat([targets, sub_targets])
            
        else:
            sub_targets = self.search(input_type, input_ids, attributes, column_names)
            targets = pd.concat([targets, sub_targets])
        
        targets = targets.reset_index(drop=True)
        return targets
