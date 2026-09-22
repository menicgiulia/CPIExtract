'''Template for a pipeline.Load,search,filter,postprocess datas from all of the databases.'''

import pandas as pd
import numpy as np
import os

from abc import ABC
from ..databases import *
from ..sql_connection import connect_to_mysql
from ..servers.BiomartServer import BiomartServer
from ..servers.MyGeneServer import MyGeneServer

class Pipeline(ABC):
    '''Template for a pipeline.Load,search,filter,postprocess datas from all of the databases.'''

    def __init__(self, execution_mode: str='local', dbs: dict|None=None, server_info: dict|None=None, 
                 pubchem_files: dict|None=None, server_select: str='mygene') -> None:
        """
        Parameters
        ----------
        execution_mode : {'local', 'online', 'server'}, default 'local'
            The execution modality of the pipeline:
            - local: retrieves most databases info from local csv files, if available
            - server: retrieves most databases info from MySQL server

        dbs : dict
            Dictionary containing all the database stored into DataFrame objects

        server_info : dict
            Dictionary containing the configuration info to connect to the SQL server
        
        pubchem_files : dict, optional
            Dictionary containing path to PubChem local database for faster processing.
            Expects local in a dubckdb structure: pubchem.duckdb
            If database not found, PubChem will use API calls.

        server_select : {'mygene', 'biomart'}, default 'mygene'
            The service used for gene identifier harmonization.
            - mygene: queries MyGene.info (recommended, more stable)
            - biomart: queries Ensembl BioMart
        """

        self.cnx = None

        if execution_mode not in ['local', 'server']:
            raise ValueError("Invalid execution mode provided. Must be one of the following: \
                             'local', 'server'.")
        
        if execution_mode in ['local']:
            if dbs is None:
                raise ValueError("Databases not specified.")
        
        if execution_mode in ['server']: 
            if server_info is None:
                raise ValueError("Server info not specified.")
            
            dbs = {}
            self.cnx = connect_to_mysql(config=server_info)

            if not self.cnx or not self.cnx.is_connected():
                raise ConnectionError("Couldn't connect to SQL server.")
        
        # Get PubChem database path with automatic detection
        pc_bioact = self._get_pubchem_path(pubchem_files)
        
        if server_select == 'mygene':
            self.gene_server = MyGeneServer()
        elif server_select == 'biomart':
            self.gene_server = BiomartServer()
        else:
            raise ValueError(f"server_select must be 'mygene' or 'biomart', got '{server_select}'")

        self.databases = {
            'pc': PubChem(bioact_file=pc_bioact,gene_server=self.gene_server),
            'chembl': ChEMBL(database=dbs.get('chembl', None), connection=self.cnx,gene_server=self.gene_server),
            'bdb': BindingDB(database=dbs.get('bdb', None), connection=self.cnx,gene_server=self.gene_server),
            'stitch': Stitch(database=dbs.get('stitch', None), connection=self.cnx,gene_server=self.gene_server),
            'ctd': CTD(database=dbs.get('ctd', None), connection=self.cnx,gene_server=self.gene_server),
            'dtc': DTC(database=dbs.get('dtc', None), connection=self.cnx,gene_server=self.gene_server),
            'otp': OTP(database=dbs.get('otp', None),chembl_database=dbs.get('chembl', None),gene_server=self.gene_server),
            'dc': DrugCentral(database=dbs.get('dc', None), connection=self.cnx,gene_server=self.gene_server),
            'db': DB(database=dbs.get('db', None), connection=self.cnx,gene_server=self.gene_server)
        }

        self.database_args = {}

        self.sources = ['PubChem', 'ChEMBL', 'BindingDB', 'Stitch', 'CTD', 'DTC', 'OTP', 'DrugCentral', 'DrugBank']

    TYPE_GROUPS = {
        # K-group: binding/dissociation/affinity constants
        'Ki': 'K', 'KI': 'K', 'Kd': 'K', 'KD': 'K', '1/KI': 'K', 'Km': 'K', 'app Km': 'K',
        'PKI': 'K', 'PKD': 'K', 'PA2': 'K', 'PKB': 'K', 'PED50': 'K',
        'LOGKI': 'K', 'LOGKD': 'K', 'LOG KI': 'K', 'LOG KD': 'K',
        '-LOGKD': 'K', 'LOG1/KD': 'K',
        # C50-group: 50%-response functional potency measures
        'IC50': 'C50', 'EC50': 'C50', 'AC50': 'C50', 'XC50': 'C50', 'GC50': 'C50',
        'GI50': 'C50', 'SC50': 'C50', 'DC50': 'C50', 'CC50': 'C50',
        'ACTIVITYEC50': 'C50', 'AVERAGEIC50': 'C50', 'Potency': 'C50',
        'PIC50': 'C50', 'PIC50(CALC)': 'C50', 'PEC50': 'C50', 'PXC50': 'C50',
        '-LOGIC50': 'C50', 'LOGIC50': 'C50', 'logIC50': 'C50', 'LOGEC50': 'C50',
        'LOG IC50': 'C50', 'LOG EC50': 'C50', 'LOG1/IC50': 'C50', 'LOG(1/IC50)': 'C50',
        'LOG(10^6/IC50)': 'C50',
    }

    def _harmonize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Extracts/renames a per-database raw result (whatever its native column names are) into
        the pipeline's unified schema. Quantitative databases (PubChem/ChEMBL/BindingDB/DTC/
        DrugCentral) output mixed-case 'pChEMBL_eq'/'pChEMBL_lt'/'pChEMBL_gt' directly.
        Qualitative-only databases (Stitch/CTD/OTP/DrugBank - no quantitative measurements at
        all) output a single lowercase 'pchembl_eq' column with no lt/gt distinction (always
        NaN in practice, since these sources have no way to derive a real value or a censored
        bound). The original 'pchembl_value' name is checked last, purely as a backward-compat
        fallback for any not-yet-migrated local data using the pre-rewrite naming.
        """
        out = pd.DataFrame(index=df.index)

        if 'pChEMBL_eq' in df.columns:
            out['pchembl_eq'] = df['pChEMBL_eq']
            out['pchembl_lt'] = df['pChEMBL_lt'] if 'pChEMBL_lt' in df.columns else np.nan
            out['pchembl_gt'] = df['pChEMBL_gt'] if 'pChEMBL_gt' in df.columns else np.nan
        elif 'pchembl_eq' in df.columns:
            out['pchembl_eq'] = df['pchembl_eq']
            out['pchembl_lt'] = df['pchembl_lt'] if 'pchembl_lt' in df.columns else np.nan
            out['pchembl_gt'] = df['pchembl_gt'] if 'pchembl_gt' in df.columns else np.nan
        else:
            out['pchembl_eq'] = df['pchembl_value'] if 'pchembl_value' in df.columns else np.nan
            out['pchembl_lt'] = np.nan
            out['pchembl_gt'] = np.nan

        for col in ['inchi', 'inchikey', 'CID', 'smiles', 'connectivity_smiles', 'iupac_name', 'synonyms',
                    'entrez', 'hgnc_symbol', 'gene_type', 'description', 'standard_type',
                    'datasource']:
            out[col] = df[col] if col in df.columns else np.nan

        if 'pChEMBL_flag' in df.columns:
            out = out.loc[~df['pChEMBL_flag'].fillna(False)]

        return out

    def _top_synonyms(self, synonyms_series: pd.Series, n: int = 10) -> str:
        """
        Collects up to n unique synonyms across all rows for one compound/protein, pipe-joined.
        Handles both representations seen in this codebase: a pipe-joined string (most
        database outputs, e.g. DrugBank) and a native Python list (compound_identifiers()'s
        own output, confirmed against real data) - a value in either form is accepted, rather
        than silently skipping list-valued entries the way an isinstance(val, str) check alone
        would.
        """
        all_syn = []
        for val in synonyms_series.dropna():
            if isinstance(val, str):
                all_syn.extend(s.strip() for s in val.split('|') if s.strip())
            elif isinstance(val, (list, tuple, set)):
                all_syn.extend(str(s).strip() for s in val if str(s).strip())
        # Preserve order of first appearance while deduplicating, then cap at n
        seen = []
        for s in all_syn:
            if s not in seen:
                seen.append(s)
        return ' | '.join(seen[:n])

    def _get_pubchem_path(self, pubchem_files):
        """
        Get PubChem database path with automatic detection.
        """
        if pubchem_files:
            # Check for direct database file specification first
            if 'db_file' in pubchem_files:
                return pubchem_files['db_file']
    
        # Try environment variable
        pubchem_dir = os.environ.get('PUBCHEM_DATA')
    
        # Fallback to default location
        if not pubchem_dir:
            pubchem_dir = os.path.join(os.path.expanduser('~'), 'CPE_data', 'pubchem')
        if not os.path.exists(pubchem_dir):
            return None
    
        # Return path to database
        db_file = os.path.join(pubchem_dir, 'pubchem.duckdb')
        if os.path.exists(db_file):
            return db_file
    
        return None

    def _aggregate_pchembl(self, data: pd.DataFrame, index: int, comp: pd.DataFrame,
                            pchembl_thres: float, pchembl_grouping: str = 'combined',
                            strong_positive_thres: float = 6.0) -> pd.DataFrame:
        """
        Aggregates pchembl_eq across all measurements for one compound-protein pair, and
        classifies the pair into a new 'interaction_class' column: 'strong positive'/'weak
        positive' (quantitative, split at strong_positive_thres)/'negative'/'ambiguous' when a
        real number places the pair on the pChEMBL scale, or 'qualitative positive'/
        'qualitative negative' when nothing does - see the classification logic below for
        exactly which case produces which label.

        pchembl_grouping controls how the average/std are computed:
            'combined' - computes BOTH a single combined average across every standard_type
                         (columns 'pchembl_count'/'ave_pchembl'/'std_pchembl', no suffix) AND
                         a separate average per type-group - K-types (Ki/Kd/...) vs C50-types
                         (IC50/EC50/...), see TYPE_GROUPS - suffixed '_K'/'_C50' - together in
                         the same output row, rather than needing to run this twice to get both
            'unique'   - a separate average per exact standard_type value present in the data;
                         columns dynamically suffixed with the type itself (e.g. '_IC50')

        """
        numeric_eq = pd.to_numeric(comp['pchembl_eq'], errors='coerce')
        has_zero = (numeric_eq == 0).any()
        real_mask = numeric_eq.notnull() & (numeric_eq != 0)
        real_values = numeric_eq[real_mask]

        if pchembl_grouping == 'combined':
            type_group = comp['standard_type'].map(self.TYPE_GROUPS)
            groups = {'': real_values}
            groups.update({f'_{g}': numeric_eq[real_mask & (type_group == g)] for g in ['K', 'C50']})
        elif pchembl_grouping == 'unique':
            groups = {f'_{t}': numeric_eq[real_mask & (comp['standard_type'] == t)]
                      for t in comp['standard_type'].dropna().unique()}
        else:
            raise ValueError(f"pchembl_grouping must be 'combined' or 'unique', got '{pchembl_grouping}'")

        for suffix, vals in groups.items():
            count_col, ave_col, std_col = f'pchembl_count{suffix}', f'ave_pchembl{suffix}', f'std_pchembl{suffix}'

            for c in (ave_col, std_col):
                if c not in data.columns:
                    data[c] = pd.Series([None] * len(data), index=data.index, dtype='object')
                elif data[c].dtype != 'object':
                    data = data.astype({c: 'object'})

            if len(vals) > 0:
                data.loc[index, count_col] = len(vals)
                data.loc[index, ave_col] = vals.mean()
                with np.errstate(invalid='ignore'):
                    data.loc[index, std_col] = np.std(vals) if len(vals) > 1 else 0.0
            else:
                data.loc[index, count_col] = 0
                data.loc[index, ave_col] = 'Sources do not provide activity data'
                data.loc[index, std_col] = 'Sources do not provide activity data'

        def _positive_tier(value: float) -> str:
            return 'strong positive' if value > strong_positive_thres else 'weak positive'

        if len(real_values) > 0:
            ave = real_values.mean()
            classification = _positive_tier(ave) if ave > pchembl_thres else 'negative'
        elif has_zero:
            classification = 'qualitative negative'
        else:
            lt_vals = pd.to_numeric(comp['pchembl_lt'], errors='coerce').dropna()
            gt_vals = pd.to_numeric(comp['pchembl_gt'], errors='coerce').dropna()
            if len(gt_vals) > 0 and gt_vals.mean() > pchembl_thres:
                classification = _positive_tier(gt_vals.mean())
            elif len(lt_vals) > 0 and lt_vals.mean() <= pchembl_thres:
                classification = 'negative'
            elif len(lt_vals) > 0 and lt_vals.mean() > pchembl_thres:
                classification = 'ambiguous'
            else:
                classification = 'qualitative positive'

        data.loc[index, 'interaction_class'] = classification

        return data