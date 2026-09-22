"""
Loads a data folder produced by prepare_cpiextract_data.py into the exact arguments
Pipeline.__init__() expects: a `dbs` dict of pandas DataFrames, keyed the same way
Pipeline.py's own `self.databases` dict is (pc/chembl/bdb/stitch/ctd/dtc/otp/dc/db), and a
`pubchem_files` dict pointing at the local PubChem DuckDB build.

Usage
-----
    from load_cpiextract_data import load_dbs, load_pubchem_files
    from cpiextract.pipelines.Comp2Prot import Comp2Prot

    dbs = load_dbs('~/CPE_data')
    pubchem_files = load_pubchem_files('~/CPE_data')
    C2P = Comp2Prot(execution_mode='local', dbs=dbs, pubchem_files=pubchem_files)

Missing files are skipped by default (skip_missing=True) rather than raising, since
Pipeline.py itself already handles a missing database gracefully (dbs.get(key, None) ->
that database's data_manager falls back to SQLManager(None, ...) internally) - matching
CPIExtract's own existing tolerance for partially-available local data, not stricter
than the pipeline itself requires. Pass skip_missing=False to fail loudly instead, e.g.
in a CI/validation context where you expect every database to be present.
"""

import os

import pandas as pd

# Maps Pipeline.py's own self.databases keys to the exact standardized filenames
# prepare_cpiextract_data.py produces. PubChem is handled separately (load_pubchem_files)
# since it isn't a DataFrame - it's a DuckDB file path.
DB_FILENAMES = {
    'chembl': 'ChEMBL_standardized.csv',
    'bdb': 'BindingDB_standardized.csv',
    'stitch': 'STITCH_standardized.csv',
    'ctd': 'CTD_standardized.csv',
    'dtc': 'DTC_standardized.csv',
    'otp': 'OTP_standardized.csv',
    'dc': 'DrugCentral_standardized.csv',
    'db': 'DrugBank_standardized.csv',
}

# Columns each database class actually filters on with an EXACT '==' or '.isin()'
# comparison (confirmed against the current database files' own retrieve_raw_data()/
# self._local_data[...] calls - see conversation), and the dtype the querying code's own
# filter_value uses. LocalManager (and OTP's own inline filtering) does no type coercion
# at all, so a CSV round-trip silently producing the wrong in-memory dtype (e.g. pandas
# inferring a numeric-looking ID column as int64, or leaving a mixed column as float64
# because of NaNs) makes every query against that column return zero rows - no error,
# just silent, wrong-looking-like-"no data" results. Applied on write (prepare_
# cpiextract_data.py, for a clean on-disk text representation with no accidental ".0"
# suffixes) and on read (here, which is what actually fixes the in-memory dtype
# LocalManager compares against - the more important of the two).
#
# Excluded: PubChem, whose local-data path bypasses LocalManager/these CSVs entirely -
# it queries the DuckDB file directly (prepare_pubchem() already declares that schema's
# column types explicitly in SQL).
LOOKUP_DTYPES: dict[str, dict[str, str]] = {
    'chembl':  {'FirstBlock': 'str', 'target_chembl_id': 'str'},
    'bdb':     {'FirstBlock': 'str', 'UniProt (SwissProt) Primary ID of Target Chain 1': 'str'},
    'ctd':     {'FirstBlock': 'str', 'GeneID': 'Int64'},  # nullable int - GeneID can be missing
    'dtc':     {'FirstBlock': 'str', 'target_id': 'str'},
    'db':      {'FirstBlock': 'str', 'HUGO Gene Nomenclature Committee (HGNC)': 'str'},
    'dc':      {'FirstBlock': 'str', 'ACCESSION': 'str'},
    'stitch':  {'FirstBlock': 'str', 'protein': 'str'},
    # OTP: 'inchikey' for compound-input lookups, 'chemblIds'/'targets' for protein-input
    # lookups - all via self._local_data[...] rather than LocalManager, same risk.
    'otp':     {'inchikey': 'str', 'chemblIds': 'str', 'targets': 'str'},
}


def load_dbs(data_path: str, skip_missing: bool = True) -> dict[str, pd.DataFrame]:
    """
    Reads every standardized CSV present in data_path into a dict of DataFrames keyed
    to match Pipeline.py's own self.databases dict - suitable to pass directly as
    Pipeline(..., dbs=load_dbs(data_path)).

    Applies LOOKUP_DTYPES on read, so the columns each database class actually filters
    on come back as the dtype its own query code expects, regardless of what pandas'
    automatic type inference on the raw CSV text would otherwise produce.
    """
    data_path = os.path.expanduser(data_path)
    dbs = {}
    missing = []
    for key, filename in DB_FILENAMES.items():
        path = os.path.join(data_path, filename)
        if os.path.exists(path):
            dtype = {col: dt for col, dt in LOOKUP_DTYPES.get(key, {}).items() if dt != 'Int64'}
            df = pd.read_csv(path, low_memory=False, dtype=dtype)
            # Int64 (nullable) columns need to go through pd.to_numeric first - read_csv's
            # own dtype= parameter can't coerce arbitrary text directly to a nullable type
            # the way it can for plain str/int64.
            for col, dt in LOOKUP_DTYPES.get(key, {}).items():
                if dt == 'Int64' and col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')
            dbs[key] = df
        else:
            missing.append(filename)

    if missing:
        message = (f"Missing {len(missing)} of {len(DB_FILENAMES)} expected files in "
                    f"{data_path}: {', '.join(missing)}")
        if skip_missing:
            print(f"Warning: {message} - continuing without them.")
        else:
            raise FileNotFoundError(f"{message}. Run prepare_cpiextract_data.py first, "
                                     f"or pass skip_missing=True to proceed without them.")

    return dbs


def load_pubchem_files(data_path: str) -> dict[str, str] | None:
    """
    Returns the pubchem_files dict Pipeline.__init__() expects (a dict with a 'db_file'
    key pointing at the local PubChem DuckDB build), or None if it hasn't been prepared
    yet - in which case PubChem falls back to live API calls, matching Pipeline.py's own
    documented behavior when no local DuckDB file is found.
    """
    data_path = os.path.expanduser(data_path)
    db_file = os.path.join(data_path, 'pubchem', 'pubchem.duckdb')
    if os.path.exists(db_file):
        return {'db_file': db_file}
    print(f"Note: no pubchem.duckdb found at {db_file} - PubChem will use live API calls "
          f"instead of the local database.")
    return None


def check_status(data_path: str) -> None:
    """Prints which of the 9 databases are ready to load and which are missing - a quick
    sanity check to run after prepare_cpiextract_data.py, or before a long pipeline run."""
    data_path = os.path.expanduser(data_path)
    print(f"Checking {data_path}...")

    pubchem_ready = os.path.exists(os.path.join(data_path, 'pubchem', 'pubchem.duckdb'))
    print(f"  {'[OK]  ' if pubchem_ready else '[MISS]'} pubchem")

    for key, filename in DB_FILENAMES.items():
        ready = os.path.exists(os.path.join(data_path, filename))
        print(f"  {'[OK]  ' if ready else '[MISS]'} {key} ({filename})")
