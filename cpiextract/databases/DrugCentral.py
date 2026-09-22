'''Loading,searching,filtering and preprocessing data from DrugCentral.'''

import numpy as np
import pandas as pd
import time

from ..utils.typing import Connection
from ..servers.PubChemServer import PubChemServer
from ..servers.BiomartServer import BiomartServer
from ..servers.MyGeneServer import MyGeneServer
from .Database import Database
from ..data_manager import *

class DrugCentral(Database):

    # ACT_TYPE values within DrugCentral
    VALID_ACT_TYPES = ['IC50', 'Ki', 'EC50', 'Kd', 'AC50', 'Km', 'app Km']

    # ACTION_TYPE binding indication in DrugCentral
    INVALID_ACTION_TYPES = ['PHARMACOLOGICAL CHAPERONE', 'RELEASING AGENT', 
                            'ANTISENSE INHIBITOR', 'ANTIBODY BINDING']

    # TARGET_CLASS removals: non single protein
    INVALID_TARGET_CLASSES = ['RNA', 'Polyprotein']

    OUTPUT_COLUMNS = ['entrez','gene_type','hgnc_symbol','description',
               'pChEMBL_eq','pChEMBL_lt','pChEMBL_gt','standard_type','datasource']

    # Final output schema for proteins()
    PROTEIN_OUTPUT_COLUMNS = ['inchi','inchikey','smiles','connectivity_smiles','iupac_name','datasource',
                'pChEMBL_eq','pChEMBL_lt','pChEMBL_gt','standard_type']

    PUBCHEM_FETCH_COLUMNS = ['inchi','inchikey','smiles','connectivity_smiles','iupac_name',
                'datasource','pchembl_value']

    def __init__(self, connection: Connection|None=None, database: pd.DataFrame|None=None, 
                gene_server=None, server_select='mygene'):

        if gene_server is not None:
            self.gene_server = gene_server
        elif server_select == 'mygene':
            self.gene_server = MyGeneServer()
        elif server_select == 'biomart':
            self.gene_server = BiomartServer()
        else:
            raise ValueError(f"server_select must be 'mygene' or 'biomart', got '{server_select}'")

        if database is not None:
            self.data_manager = LocalManager(database)
        else:
            self.data_manager = SQLManager(connection, 'DC')

    def _filter_database(self, dc_raw: pd.DataFrame, dc_extra: bool) -> pd.DataFrame:

        has_valid_act_type = dc_raw['ACT_TYPE'].isin(self.VALID_ACT_TYPES)
        is_bare_positive_candidate = dc_raw['ACT_VALUE'].isna() & dc_raw['ACTION_TYPE'].notnull()
        dc_filt = dc_raw.loc[(dc_raw['ORGANISM']=='Homo sapiens') &
                             (has_valid_act_type | is_bare_positive_candidate) &
                             (dc_raw['ACT_SOURCE'] != 'UNKNOWN')]\
                                .drop_duplicates(ignore_index=True).copy()
        if not dc_extra:
            # ACTION_TYPE exclusion only applies when data is present 
            dc_filt = dc_filt.loc[(~dc_filt['ACTION_TYPE'].isin(self.INVALID_ACTION_TYPES)) & 
                                  (~dc_filt['TARGET_CLASS'].isin(self.INVALID_TARGET_CLASSES))].copy()

        # ACCESSION contains complexes, remove for single protein interactions
        dc_filt = dc_filt[~dc_filt['ACCESSION'].astype(str).str.contains('|', regex=False)].copy()

        return dc_filt


    def _compute_pchembl(self, dc_dat: pd.DataFrame, pChEMBL_thres: float) -> pd.DataFrame:
        """
        Derives pChEMBL_eq/pChEMBL_lt/pChEMBL_gt from ACT_VALUE using RELATION. 
        no directional symbols is treated as '='.

        Bare positives via statements in ACTION_TYPE. Kept as presence-only evidence
        no true negative statements exist.
        """
        has_value = dc_dat['ACT_VALUE'].notnull() & (dc_dat['ACT_VALUE'] != 0)

        relation = dc_dat['RELATION'].fillna('=').replace(['', '~', '-'], '=')
        derived = -np.log10(dc_dat['ACT_VALUE'].where(has_value) / 1e9)

        is_gt = relation.isin(['>', '>=']) & has_value
        is_lt = relation.isin(['<', '<=']) & has_value
        is_eq = (relation == '=') & has_value

        dc_dat = dc_dat.copy()
        dc_dat['pChEMBL_eq'] = derived.where(is_eq)
        dc_dat['pChEMBL_lt'] = derived.where(is_gt)
        dc_dat['pChEMBL_gt'] = derived.where(is_lt)

        dc_dat['standard_type'] = dc_dat['ACT_TYPE']

        # A '>'-censored row is negative evidence if the bound itself rules out an active result
        true_negatives = dc_dat['pChEMBL_lt'].notnull() & (dc_dat['pChEMBL_lt'] <= pChEMBL_thres)
        # A positive only counts as evidence if its lower bound clears the threshold.
        censored_positive = dc_dat['pChEMBL_gt'].notnull() & (dc_dat['pChEMBL_gt'] > pChEMBL_thres)
        # An exact value is kept regardless of which side of threshold it falls on
        exact_value = dc_dat['pChEMBL_eq'].notnull()
        # No derivable value at all, but a real ACTION_TYPE was recorded kept as binary information
        bare_positive = (~has_value) & dc_dat['ACTION_TYPE'].notnull()

        dc_act = dc_dat.loc[exact_value | true_negatives | censored_positive | bare_positive].reset_index(drop=True)

        return dc_act

    def compounds(self, input_comp: pd.DataFrame, dc_extra: bool=False, pChEMBL_thres: float=3.0, 
                ) -> tuple[pd.DataFrame, str, pd.DataFrame]:
        """
        Retrieves proteins from DrugCentral database interacting with compound passed as input.

        Steps
        -----
        - Filters dc database to obtain proteins interacting with input compound \\
        Constraints:
            - Only Homo sapiens interactions
            - Activity Type must be one of the following: IC50, Ki, EC50, Kd, AC50
            - Action Type must not be: PHARMACOLOGICAL CHAPERONE, RELEASING AGENT
            - Target Class must not be: CD molecules, RNA, Unclassified, Viral envelope protein, Polyprotein
        - Uses activity column to compute pchembl values.
        
        Parameters
        ----------
        input_comp : DataFrame
            Dataframe of input compound data from which interacting proteins are found
        dc_extra: bool
            bool to select whether to include possibly non-Homo sapiens interactions
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid

        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values: \\
            entrez, gene_type, hgnc_symbol, description, datasource (DrugCentral), pChEMBL_eq/lt/gt, standard_type
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all DrugCentral info about the input compound
        """

        columns=self.OUTPUT_COLUMNS
        # Create an empty DataFrame with the specified columns
        dc_act = pd.DataFrame(columns=columns)
        dc_raw = pd.DataFrame()
        # Drop null values of inchikey, to make sure the first line of the dataframe has the inchi key 
        input_comp = input_comp.dropna(subset=['inchikey']).reset_index(drop=True)
        # Check if there are any input compounds remaining  
        if len(input_comp) > 0:
            input_comp_id = input_comp['inchikey_fb'][0]
            dc_raw = self.data_manager.retrieve_raw_data('FirstBlock', input_comp_id)

            if len(dc_raw) > 0:
                # Filter database
                dc_filt = self._filter_database(dc_raw, dc_extra)
                # Compute pchembl values
                dc_act = self._compute_pchembl(dc_filt, pChEMBL_thres)
                
                if len(dc_act) > 0:
                    # Unify gene identifiers by Harmonizing IDs
                    ensembl = self.gene_server

                    # BindingDB uses Uniprot IDs
                    input_type='uniprotswissprot' 
                    attributes = ['uniprotswissprot', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                    names = ['uniprot','entrez','gene_type','hgnc_symbol','description']
                    dc_targets = pd.DataFrame(columns=names)
                    
                    input_genes = list(dc_act['ACCESSION'])
                    
                    dc_targets = ensembl.subset_search(input_type, input_genes, attributes, names)
    
                    # For each protein, assign specific biomart column values to the ones from the original dc database
                    for index, row in dc_act.iterrows():
                        # Find the compound using uniprot
                        S1 = dc_targets.loc[dc_targets['uniprot'] == row['ACCESSION']]
                        if len(S1) > 0:
                            dc_act.loc[index, 'entrez'] = S1['entrez'].iloc[0]
                            dc_act.loc[index, 'gene_type'] = S1['gene_type'].iloc[0]
                            dc_act.loc[index, 'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                            dc_act.loc[index, 'description'] = S1['description'].iloc[0]
                            dc_act.loc[index, 'note'] ='Harmonized gene ID'
                        else:
                            dc_act.loc[index, 'entrez'] = None
                            dc_act.loc[index, 'gene_type'] = None
                            dc_act.loc[index, 'hgnc_symbol'] = None
                            dc_act.loc[index, 'description'] = None
                            dc_act.loc[index, 'note'] ='Failed to harmonize gene ID'
                    dc_act['datasource'] = 'DrugCentral'
                    statement = 'completed'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input compound does not contain inchikey'

        return dc_act, statement, dc_raw
    

    def proteins(self, input_protein: pd.DataFrame, dc_extra: bool=False, pChEMBL_thres: float=3.0, 
                ) -> tuple[pd.DataFrame, str, pd.DataFrame]:
        """
        Retrieves compounds from DrugCentral database interacting with proteins passed as input.

        Steps
        -----
        - Filters dc database to obtain compounds interacting with input proteins:
        Constraints:
            - Only Homo sapiens interactions
            - Activity Type must be one of the following: IC50, Ki, EC50, Kd, AC50
            - Action Type must not be: PHARMACOLOGICAL CHAPERONE, RELEASING AGENT
            - Target Class must not be: CD molecules, RNA, Unclassified, Viral envelope protein, Polyprotein
        - Uses activity column to compute pchembl values. \\
        - Uses Pubchempy to obtain compounds info (modified with data from original dc database) to return,
        searching with InChiKey
        
        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compound are found
        dc_extra: bool
            bool to select whether to include possibly non-Homo sapiens interactions
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid       

        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values: \\
            entrez, gene_type, hgnc_symbol, description, datasource (DrugCentral), pChEMBL_eq/lt/gt, standard_type
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all DrugCentral info about the input compound
        """

        columns = self.PROTEIN_OUTPUT_COLUMNS
        # Create an empty DataFrame with the specified columns
        dc_c1 = pd.DataFrame(columns=columns)
        dc_raw = pd.DataFrame()
        # Drop duplicated and null values of uniprot, to make sure the first line of the dataframe has the uniprot ID
        input_protein= input_protein.dropna(subset=['uniprot'])
        input_protein= input_protein.drop_duplicates(subset='uniprot').reset_index(drop=True) 
        # Check if there are any input proteins remaining       
        if len(input_protein) > 0:
            input_protein_id=input_protein['uniprot'].iloc[0]
            # Select only compounds interacting with the input protein
            dc_raw = self.data_manager.retrieve_raw_data('ACCESSION', input_protein_id)
            if len(dc_raw):
                # Filter database
                dc_c = self._filter_database(dc_raw, dc_extra)
                
                # Compute pchembl values
                dc_act = self._compute_pchembl(dc_c, pChEMBL_thres)
                # Check if at least one interacting compound has been found
                if len(dc_act) > 0:
                    pc = PubChemServer()

                    if 'CID' in dc_act.columns:
                        # Retrieve compounds using cids
                        compounds = self._pubchem_search_cid(dc_act, self.PUBCHEM_FETCH_COLUMNS, pc)
                        
                        if len(compounds) > 0:
                            dc_info = dc_act[['CID', 'ACCESSION', 'TARGET_CLASS', 'ACT_COMMENT',
                                                'pChEMBL_eq', 'pChEMBL_lt', 'pChEMBL_gt', 'standard_type']].\
                                        rename(columns={'ACCESSION': 'uniprotid'}).\
                                                        drop_duplicates()
                            
                            dc_c1 = pd.merge(compounds, dc_info, on='CID', how='left')
                    else:
                        # Filter only columns from pubchempy to return
                        selected_columns = pc.get_columns(self.PUBCHEM_FETCH_COLUMNS[:-2])
                        compounds = pd.DataFrame()
                        ids = list(dc_act['inchikey'].unique())
                        for id in ids:
                            try:
                                comp = pc.get_compounds(id, selected_columns, namespace='inchikey')
                                compounds = pd.concat([compounds, comp])
                            except:
                                None
                            time.sleep(0.5)
                        # Check if at least a match has been found
                        if len(compounds) > 0:
                            dc_c1 = compounds.rename(columns=pc.properties)

                            # Add additional values from activity dataframe
                            dc_info = dc_act[['inchikey', 'ACCESSION', 'TARGET_CLASS', 'ACT_COMMENT',
                                                'pChEMBL_eq', 'pChEMBL_lt', 'pChEMBL_gt', 'standard_type']].\
                                        rename(columns={'ACCESSION': 'uniprotid'}).\
                                                        drop_duplicates()
                            
                            dc_c1 = pd.merge(dc_c1, dc_info, on='inchikey', how='left')

                    if len(dc_c1) > 0:
                        dc_c1['datasource'] = 'DrugCentral'
                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using PubChem'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein does not contain UniProt ID'
        return dc_c1, statement, dc_raw