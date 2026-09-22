'''Loading,searching,filtering and preprocessing data from DTC.'''

import pandas as pd
import numpy as np

from ..utils.typing import Connection
from ..servers.PubChemServer import PubChemServer
from ..servers.BiomartServer import BiomartServer
from ..servers.MyGeneServer import MyGeneServer
from .Database import Database
from ..data_manager import *

class DTC(Database):

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
            self.data_manager = SQLManager(connection, 'DTC')

    # Standard Type equivalent to pChEMBL (some in nM and some in M)
    DIRECT_STD_TYPES = ['IC50','KI','EC50','KD','Kd','ED50','XC50','AC50','Ki','CC50','AVERAGEIC50',
                         'ACTIVITYEC50','SC50','AD50','DC50','XI50','GC50','GI50',
                         'LOGEC50','LOGIC50','LOGKD','LOGKI','logIC50',
                         'LOG KD','LOG KI','LOG EC50','1/KI']
    # Inverse pChEMBL equivalent standard types (some in nM and some in M)
    INVERTED_STD_TYPES = ['-LOGIC50','PIC50','-LOGKD','PEC50','PIC50(CALC)','-LOGED50',
                           'LOG1/IC50','PED50','LOG1/KD','LOG(1/IC50)','LOG(10^6/IC50)','PXC50',
                           '-LOG IC50','PKI','PKD','PA2','PKB']
    # Log-scaled versions of pChEMBL equivalent standard types (some in nM and some in M)
    PSCALE_DIRECT_TYPES = ['-LOGIC50','PIC50','-LOGKD','PEC50','PIC50(CALC)','-LOGED50','PED50','PXC50',
                            '-LOG IC50','PKI','PKD','PA2','PKB',
                            'LOG1/IC50','LOG1/KD','LOG(1/IC50)']
    # Giving unitless rows an expected unit against real DTC data, seeing if in a physically plausible range (0-11):
    #   LOGIC50 (uM, 96.7% plausible), LOGKI (uM, 100%), LOGEC50 (uM, 96.6%),
    #   LOGKD (nM, 100%), LOG KD (Molar, 100%), LOG KI (nM, 75.2%), LOG EC50 (uM, 95.8%).
    LOG_DIRECT_OFFSETS = {'LOGIC50': 6, 'logIC50': 6, 'LOGKI': 6, 'LOGEC50': 6,
        'LOGKD': 9, 'LOG KI': 9,'LOG KD': 0,'LOG EC50': 6}
    # standard_units that are inverse concentration
    INVERTED_UNIT_TYPES = ['/NM', '10^5/M', 'M-1', '/UM', "10'-9L/MOL", "10'-10L/MOL"]
    # Whitelist of standard_units accepted for measurements
    VALID_UNIT_TYPES = ['NM','M','/NM','NMOL/L','MG.KG-1','UG.ML-1','UMOL.KG-1','P.P.M.','PPM',
                    'UG KG-1','MG KG-1','MG/ML','NG ML-1',"10'-4UMOL/L",'UG ML-1','MM','NMOL/MG',
                    "10'13NM","10'-10M","10'8NM",'10^5/M','MOL/G',"10'-8M",'NM G-1','M-1',"10'-5M",
                    "10'7NM",'/UM',"10'-7M","10'-9M","10'-9L/MOL","10'-10L/MOL",'PM','UM L-1','UG/ML',
                    'NO_UNIT']
    # Detection-only plausibility band for a computed pChEMBL value in Molar Concentration
    # Derived from ChEMBL and BindingDB datasets
    PCHEMBL_PLAUSIBLE_RANGE = (0,11)

    def _filter_database(self, DTC_raw: pd.DataFrame) -> pd.DataFrame:
        valid_standard_types = self.DIRECT_STD_TYPES + self.INVERTED_STD_TYPES
        dimensionless_types = set(self.PSCALE_DIRECT_TYPES) | set(self.LOG_DIRECT_OFFSETS) | {'LOG(10^6/IC50)'}

        unit_ok = (DTC_raw['standard_units'].isin(self.VALID_UNIT_TYPES) |
                   (DTC_raw['standard_units'].isnull() & DTC_raw['standard_type'].isin(dimensionless_types)))

                            # Keep only valid standard types
        DTC_filt = DTC_raw.loc[DTC_raw['standard_type'].isin(valid_standard_types) & 
                            # Keep only valid standard units or dimensionless types
                            unit_ok &
                            # Remove entries with no measurement
                            DTC_raw['standard_value'].notnull() &
                            # Remove entries without target info
                            DTC_raw['target_id'].notnull() &
                            # Remove entries without compound id
                            DTC_raw['compound_id'].notnull()
                            ].drop_duplicates().copy().reset_index(drop=True)  
        return DTC_filt
        

    def _flag_implausible_pchembl(self, DTC_act: pd.DataFrame) -> pd.DataFrame:
        """
        Adds a persistent 'pChEMBL_flag' column (True/False) marking any row whose computed
        pChEMBL_eq/lt/gt falls outside PCHEMBL_PLAUSIBLE_RANGE
        """
        lo, hi = self.PCHEMBL_PLAUSIBLE_RANGE
        flagged = pd.Series(False, index=DTC_act.index)
        for col in ['pChEMBL_eq', 'pChEMBL_lt', 'pChEMBL_gt']:
            if col not in DTC_act.columns:
                continue
            vals = DTC_act[col]
            bad = vals.notnull() & ((vals < lo) | (vals > hi))
            flagged = flagged | bad
            if not (getattr(self, '_verbose', False) and bad.any()):
                continue
            for idx in DTC_act.loc[bad].index:
                print("Warning: implausible standard_units is likely for some rows. Check pChEMBL_flag column.")

        DTC_act['pChEMBL_flag'] = flagged
        return DTC_act

    def _standardize_database(self, DTC_raw: pd.DataFrame, dtc_mutated: bool, pChEMBL_thres: float) -> pd.DataFrame:

        # Only use interactions with reported sources
        DTC_act = DTC_raw[DTC_raw['doc_type'].notnull()].reset_index(drop=True)
        # 'inconclusive' carries no evidence either way and is dropped.
        DTC_act = DTC_act[~DTC_act['activity_comment'].isin(['inconclusive'])].reset_index(drop=True)

        if not dtc_mutated:
            # Remove all mutated data
            DTC_act = DTC_act[(DTC_act['wildtype_or_mutant'] != 'mutated') & 
                        (DTC_act['mutation_info'].isnull())].reset_index(drop=True)

        # Converts all the measured values into nM for the calculation of pChEMBL
        DTC_act = self._standard_converter(DTC_act)
        # Convert 0 and infinite values to nan (these come from implausibly extreme values)
        DTC_act['standardized_val'] = DTC_act['standardized_val'].replace([np.inf, 0, -np.inf], np.nan)
        # Convert strings into numeric values 
        DTC_act['standardized_val'] = DTC_act['standardized_val'].apply(pd.to_numeric, errors='coerce')

        # Get Relation information
        type_inverted = DTC_act['standard_type'].isin(self.INVERTED_STD_TYPES)
        unit_inverted = DTC_act['standard_units'].isin(self.INVERTED_UNIT_TYPES)
        inverted = type_inverted != unit_inverted
        relation = DTC_act['standard_relation'].fillna('=').replace('', '=')

        # Filter rows by presence of affinity values and calculate pChEMBL
        has_value = DTC_act['standardized_val'].notnull() & (DTC_act['standardized_val'] > 0)
        derived = -np.log10(DTC_act['standardized_val'].where(has_value) * 1e-9)

        # Assign relational value rows
        # '>' -> upper bound (weaker), '<' -> lower bound (stronger)
        is_gt = relation.isin(['>', '>=']) & has_value
        is_lt = relation.isin(['<', '<=']) & has_value
        is_eq = (relation == '=') & has_value
        upper_bound = (is_gt & ~inverted) | (is_lt & inverted)
        lower_bound = (is_lt & ~inverted) | (is_gt & inverted)
        DTC_act['pChEMBL_eq'] = derived.where(is_eq)
        DTC_act['pChEMBL_lt'] = derived.where(upper_bound)
        DTC_act['pChEMBL_gt'] = derived.where(lower_bound)

        DTC_act = self._flag_implausible_pchembl(DTC_act)

        # Use relation to determine activity. A 'Not Active' comment is always kept regardless of any derivable bound.
        comment_negative = DTC_act['activity_comment'] == 'Not Active'
        true_negatives = comment_negative | (DTC_act['pChEMBL_lt'].notnull() & (DTC_act['pChEMBL_lt'] <= pChEMBL_thres))
        # A censored positive only counts as evidence if its lower bound clears the threshold.
        censored_positive = DTC_act['pChEMBL_gt'].notnull() & (DTC_act['pChEMBL_gt'] > pChEMBL_thres)
        # An exact value is kept regardless of which side of threshold it falls on
        exact_value = DTC_act['pChEMBL_eq'].notnull()

        DTC_act = DTC_act.loc[exact_value | true_negatives | censored_positive].reset_index(drop=True)

        # Bare 'Not Active' calls with no derivable bound at all get an explicit 0.
        missing_val_negative = (
            (DTC_act['activity_comment'] == 'Not Active') &
            DTC_act['pChEMBL_eq'].isnull() &
            DTC_act['pChEMBL_lt'].isnull()
        )
        DTC_act.loc[missing_val_negative, 'pChEMBL_eq'] = 0

        return DTC_act

    def compounds(self, input_comp: pd.DataFrame, dtc_mutated: bool=False, 
                     pChEMBL_thres: float=3.0, verbose: bool=False) -> tuple[pd.DataFrame, str, pd.DataFrame]:

        """
        Retrieves proteins from DTC database interacting with compound passed as input.

        Steps
        -----
        - Filters input database from unnecessary data \\
        Constraints:
            - Only valid standard types 
            - Only valid standard units
            - Entries with measurement only
            - Entiries with target info only \\
        - Finds matches on the database for the input compound and applies a second filter \\
        Constraints:
            - Only interactions with reported sources
            - No inconclusive activities
            - No mutated data (unless dtc_mutated=True) \\
        - Applies a standard conversion for all values required to compute pChEMBL, then derives
          pChEMBL_eq/pChEMBL_lt/pChEMBL_gt from standard_relation, classified against pChEMBL_thres.
        
        Parameters
        ----------
        input_comp : DataFrame
            Dataframe of input compounds from which interacting proteins are found
        dtc_mutated: bool
            bool to select whether to included mutated targets interactions
        pChEMBL_thres : float
            pChEMBL value used to classify a resolved interaction as positive (above) or
            negative (at or below); censored bounds are only kept when they're tight enough to
            confidently fall on one side of this value.
        verbose : bool
            if True, prints a warning for any computed pChEMBL value outside a wide plausible
            range (see PCHEMBL_PLAUSIBLE_RANGE) - a detection-only safety net for unit-labeling
            issues not covered by today's hardcoded fixes, not an automatic correction.
        
        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values: \\
            entrez, gene_type, hgnc_symbol, description, datasource (DTC), pChEMBL_eq/lt/gt, standard_type
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all DTC info about the input compound
        """

        self._verbose = verbose

        columns = ['entrez','gene_type','hgnc_symbol','description','pChEMBL_eq','pChEMBL_lt','pChEMBL_gt','pChEMBL_flag','standard_type','datasource']
        # Create an empty DataFrame with the specified columns
        DTC_act = pd.DataFrame(columns=columns)
        DTC_raw = pd.DataFrame(columns=columns)

        input_comp = input_comp.dropna(subset=['inchikey']).reset_index(drop=True)
        if len(input_comp) == 0:
            return DTC_act, 'Input compound does not contain inchi key', DTC_raw

        input_comp_id = input_comp['inchikey_fb'][0]
        DTC_raw = self.data_manager.retrieve_raw_data('FirstBlock', input_comp_id)

        # Check if at least one match has been found
        if len(DTC_raw) > 0:

            DTC_filt = self._filter_database(DTC_raw)
            DTC_filt['molecular_weight'] = input_comp['molecular_weight'].iloc[0]

            # Filter database
            DTC_act = self._standardize_database(DTC_filt, dtc_mutated, pChEMBL_thres)
                
            # Note: DTC has no tax_id or species column, can only filter for human via biomart conversion
            if len(DTC_act) > 0:
                # Unify gene identifiers by Harmonizing IDs
                ensembl = self.gene_server

                attributes = ['uniprotswissprot', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                names = ['uniprot','entrez','gene_type','hgnc_symbol','description']
                # Search by uniprot id match to biomart
                input_type = 'uniprotswissprot'
                
                dtc_targets = pd.DataFrame(columns=names)
                    
                input_genes = list(DTC_act['target_id'])
                
                dtc_targets = ensembl.subset_search(input_type, input_genes, attributes, names)

                # For each compound, assign specific biomart column values to the ones from the original DTC database
                for index, row in DTC_act.iterrows():
                    S1 = dtc_targets.loc[dtc_targets['uniprot']==row['target_id']]
                    if len(S1) > 0:
                        DTC_act.loc[index,'entrez'] = S1['entrez'].iloc[0]
                        DTC_act.loc[index,'gene_type'] = S1['gene_type'].iloc[0]
                        DTC_act.loc[index,'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                        DTC_act.loc[index,'description'] = S1['description'].iloc[0]
                        DTC_act.loc[index, 'note'] ='Harmonized gene ID'
                    else:
                        DTC_act.loc[index, 'entrez'] = None
                        DTC_act.loc[index, 'gene_type'] = None
                        DTC_act.loc[index, 'hgnc_symbol'] = None
                        DTC_act.loc[index, 'description'] = None  
                        DTC_act.loc[index, 'note'] ='Failed to harmonize gene ID'
                DTC_act['datasource']='DTC'
                statement='completed'
            else:
                statement='Filter reduced interactions to 0'
        else:
            statement='No interaction data'
    
        return DTC_act, statement, DTC_raw

    def proteins(self, input_protein: pd.DataFrame, dtc_mutated: bool=False, pChEMBL_thres: float=3.0, 
                verbose: bool=False) -> tuple[pd.DataFrame, str, pd.DataFrame]:

        """
        Retrieves compounds from DTC database interacting with proteins passed as input.

        Steps
        -----
        - Retrieves gene data from DTC database to find compounds interacting with gene
        Constraints:
            - Standard type and units not null
            - pChEMBL value > threshold
        - Uses Pubchempy to retrieve the compound info (with additional data from chembl) to return, 
        searching with compound ID or with inchi (obtained from chembl API) if the first search is unsuccessful.
        - Applies a standard conversion for all values required to compute pChembl.
        
        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compound are found
        DTC_data : DataFrame
            Dataframe containing all DTC database info
        dtc_mutated: bool
            bool to select whether to included mutated targets interactions
        pChEMBL_thres : integer
            Threshold for pChEMBL value to be valid
        verbose : bool
            if True, prints a warning for any computed pChEMBL value outside a wide plausible
            range (see PCHEMBL_PLAUSIBLE_RANGE) - a detection-only safety net for unit-labeling
            issues not covered by today's hardcoded fixes, not an automatic correction.

        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values: \\
            inchi, inchikey, smiles, iupac_name, datasource (DTC), pChEMBL_eq/lt/gt, standard_type
        """

        self._verbose = verbose

        columns = ['inchi','inchikey','smiles','connectivity_smiles','iupac_name','molecular_weight','datasource',
                   'pChEMBL_eq','pChEMBL_lt','pChEMBL_gt','pChEMBL_flag','standard_type']

        pubchem_columns = ['inchi','inchikey','smiles','connectivity_smiles','iupac_name','molecular_weight',
                            'datasource','pchembl_value']
        # Create an empty DataFrame with the specified columns
        DTC_c1 = pd.DataFrame(columns=columns) 
        DTC_raw = pd.DataFrame()
        # Drop duplicated and null values of uniprot
        input_protein = input_protein.dropna(subset=['uniprot'])
        input_protein = input_protein.drop_duplicates(subset='uniprot').reset_index(drop=True)  
        # Check if there are any input proteins remaining      
        if len(input_protein) > 0:
            input_protein_id = input_protein['uniprot'][0]

            # Filter DTC by protein of interest using uniprot ID
            DTC_raw = self.data_manager.retrieve_raw_data('target_id', input_protein_id)
            # Check if there are any compounds remaining
            if len(DTC_raw) > 0:

                DTC_filt = self._filter_database(DTC_raw)

                if len(DTC_filt) > 0:

                    pc = PubChemServer()
                    DTC_c = pd.DataFrame()
                    
                    if 'CID' in DTC_filt.columns:
                        # Retrieve compounds using cids
                        compounds = self._pubchem_search_cid(DTC_filt, pubchem_columns, pc)

                        if len(compounds) > 0:
                            # Add additional information to interacting compounds and standardize its values       
                            DTC_c = pd.merge(compounds, DTC_filt, on='CID', how='left')
                    else:
                        # Retrieve compounds using compound id (chembl ids) and inchis if formers miss
                        compounds = self._pubchem_search_chembl(DTC_filt, 'compound_id', pubchem_columns, pc)

                        if len(compounds) > 0:
                            # Add additional information to interacting compounds and standardize its values       
                            DTC_c = pd.merge(compounds, DTC_filt, left_on='id', right_on='compound_id', how='left')
                    
                    # Check if there are any compounds remaining
                    if len(DTC_c) > 0:
                        # Standardize database
                        DTC_act = self._standardize_database(DTC_c, dtc_mutated, pChEMBL_thres)

                        if len(DTC_act) > 0:
                            required_cols = ['inchi','inchikey','smiles','connectivity_smiles','iupac_name','molecular_weight',
                                        'CID', 'compound_id', 'standard_type', 'standard_value', 'standard_units', 
                                        'pChEMBL_eq', 'pChEMBL_lt', 'pChEMBL_gt', 'pChEMBL_flag', 'activity_comment']
                            for col in required_cols:
                                if col not in DTC_act.columns:
                                    DTC_act[col] = None
                            
                            # Extend compounds information
                            DTC_c1 = DTC_act[required_cols].drop_duplicates()
                            
                            # Add additional information
                            DTC_c1['datasource'] = 'DTC'
                            
                            # Ensure final columns match expected format
                            for col in columns:
                                if col not in DTC_c1.columns:
                                    DTC_c1[col] = None
                            
                            DTC_c1 = DTC_c1[columns]
                            
                            statement = 'completed'
                        else:
                            statement = 'Standardization reduced interactions to 0'
                    else:
                        statement = 'Compounds not found using Pubchem and ChEMBL'
                else:
                    statement = 'First filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein does not contain uniprot id'                         
        
        return DTC_c1, statement, DTC_raw


    def _standard_converter(self, DTC_act: pd.DataFrame) -> pd.DataFrame:

        """
        Converts all types from DTC database to those needed for pCHEMBL by storing the respective
        standardized value. Also converts all units to nM as it is required to compute pCHEMBL.

        Parameters
        ----------
        DTC_act : DataFrame
            DTC dataframe with compounds whose types and units need to be standardized

        Returns
        -------
        DataFrame
            DTC dataframe of compounds with standardized types and unit
        """    

        std_types = set(self.DIRECT_STD_TYPES) - set(self.LOG_DIRECT_OFFSETS)
        
        # Converting standard types
        DTC_act['standardized_val'] = None
        D1 = DTC_act.loc[DTC_act['standard_type'].isin(std_types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = DTC_act['standard_value'][h]
        D1 = DTC_act.loc[DTC_act['standard_type'].isin(self.PSCALE_DIRECT_TYPES)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = 10**(-DTC_act['standard_value'][h]) * 1e9
        # LOG_DIRECT_OFFSETS's per-type offset was derived assuming standard_units is blank
        # Class has comment on chosen offsets
        D1 = DTC_act.loc[DTC_act['standard_type'].isin(self.LOG_DIRECT_OFFSETS)]
        for h in D1.index:
            offset = self.LOG_DIRECT_OFFSETS[DTC_act['standard_type'][h]]
            DTC_act.loc[h,'standardized_val'] = 10**(DTC_act['standard_value'][h] - offset) * 1e9
        D1 = DTC_act.loc[DTC_act['standard_type']=='LOG(10^6/IC50)']
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = 10**(6-(DTC_act['standard_value'][h]))
        
        # Converting all standard units to nM
        types = ['NM','NMOL/L']
        D1 = DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]
        types = ['M']
        D1 = DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = DTC_act['standardized_val'][h]*(1e9)
        types = ['/NM']
        D1 = DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = 1/DTC_act['standardized_val'][h]
        types = ['MM']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e9/1e3)
        types=['10^5/M']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=((1/(DTC_act['standardized_val'][h]*(1e5)))*1e9)
        types=['NM G-1']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e3)
        types=['M-1']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=((1/DTC_act['standardized_val'][h])*1e9)
        types=['/UM']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=(1/(DTC_act['standardized_val'][h]))*(1e9/1e6)
        types=['PM']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e9/1e12)
        types=['MOL/G']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e3/1)*(1e9/1)
        types=['NMOL/MG']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e3/1e-3)
        types=['UM L-1','UMOL.KG-1']  
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e9/1e6)
        
        types=['UG.ML-1','UG ML-1','UG/ML']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            if D1['molecular_weight'][h]=='':
                DTC_act.loc[h,'standardized_val']=np.nan
            else:
                DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e3/1)*(1/1e6)*(1/float(D1['molecular_weight'][h]))*(1e9/1)
        types=['MG/ML']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            if D1['molecular_weight'][h]=='':
                DTC_act.loc[h,'standardized_val']=np.nan
            else:
                DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e3/1)*(1/1e3)*(1/float(D1['molecular_weight'][h]))*(1e9/1)
        types=['NG ML-1']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            if D1['molecular_weight'][h]=='':
                DTC_act.loc[h,'standardized_val']=np.nan
            else:
                DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e3/1)*(1/1e9)*(1/float(D1['molecular_weight'][h]))*(1e9/1)
        types=['MG.KG-1','MG KG-1']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            if D1['molecular_weight'][h]=='':
                DTC_act.loc[h,'standardized_val']=np.nan
            else:
                DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1/1e3)*(1/float(D1['molecular_weight'][h]))*(1e9/1)
        types=['UG KG-1']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            if D1['molecular_weight'][h]=='':
                DTC_act.loc[h,'standardized_val']=np.nan
            else:
                DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1/1e6)*(1/float(D1['molecular_weight'][h]))*(1e9/1)
        types=['P.P.M.','PPM']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            if D1['molecular_weight'][h]=='':
                DTC_act.loc[h,'standardized_val']=np.nan
            else:
                DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1/float(D1['molecular_weight'][h]))*(1/1e3)*(1e9/1)
        types=["10'13NM","10'8NM","10'7NM"]
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            E=DTC_act['standard_units'][h][DTC_act['standard_units'][h].find("'")+1:DTC_act['standard_units'][h].find("N")]
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(10**int(E))
        types=["10'-10M","10'-8M","10'-5M","10'-7M","10'-9M"]
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            E=DTC_act['standard_units'][h][DTC_act['standard_units'][h].find("'")+1:DTC_act['standard_units'][h].find("M")]
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(10**int(E))*(1e9)
        types=["10'-4UMOL/L"]
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            E=DTC_act['standard_units'][h][DTC_act['standard_units'][h].find("'")+1:DTC_act['standard_units'][h].find("U")]
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(10**int(E))*(1e9/1e6)
        types=["10'-9L/MOL","10'-10L/MOL"]
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            E=DTC_act['standard_units'][h][DTC_act['standard_units'][h].find("'")+1:DTC_act['standard_units'][h].find("L")]
            DTC_act.loc[h,'standardized_val']=(1/(DTC_act['standardized_val'][h]*(10**int(E))))*(1e9)
        
        return DTC_act