"""
Downloads and standardizes source data for all 9 CPIExtract databases into a single,
flat data folder that the pipeline can load directly (see load_cpiextract_data.py).

Automated end-to-end (no manual download needed):
    - PubChem     (FTP downloads + local DuckDB build)
    - ChEMBL      (FTP download of the SQLite dump + SQL extraction)
    - OTP         (web-scraped Open Targets Platform parquet files)

Require a manually-downloaded raw file first (licensing/access restrictions, or no
stable scriptable URL) - this script checks for the file, standardizes it if present,
and prints exactly what to download and where to place it if not:
    - BindingDB
    - CTD
    - DrugBank    (requires an academic license - see https://go.drugbank.com/releases/latest)
    - STITCH
    - DTC
    - DrugCentral

Usage
-----
    python -m cpiextract.utils.prepare_cpiextract_data --data-path ~/CPE_data
    python -m cpiextract.utils.prepare_cpiextract_data --data-path ~/CPE_data --skip drugbank stitch
    python -m cpiextract.utils.prepare_cpiextract_data --data-path ~/CPE_data --only pubchem chembl

Re-running is safe: raw downloads are skipped if already present (use --force-download
to re-fetch), and standardization always re-runs against whatever raw data is on disk,
since it's cheap relative to the downloads themselves.

Output (flat folder at --data-path):
    pubchem/pubchem.duckdb              - PubChem's local bioactivity database
    ChEMBL_standardized.csv
    BindingDB_standardized.csv
    CTD_standardized.csv
    DrugBank_standardized.csv
    DrugCentral_standardized.csv
    DTC_standardized.csv
    STITCH_standardized.csv
    OTP_standardized.csv
    pubchem-source-SIDs.csv             - PubChem Source+ExternalID -> CID reference
    pubchem-cid-inchikey.csv            - PubChem CID -> InChIKey reference
"""

import argparse
import glob
import gzip
import os
import re
import shutil
import sqlite3
import tarfile
import xml.etree.ElementTree as ET
import zipfile

import duckdb
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm

try:
    # Normal case: run as part of the package (python -m cpiextract.utils.prepare_cpiextract_data)
    from .load_cpiextract_data import LOOKUP_DTYPES
except ImportError:
    # Fallback: this file executed directly by path (python path/to/prepare_cpiextract_data.py),
    # which has no package context for a relative import - load_cpiextract_data.py is a sibling
    # in the same folder in that case too, just imported as a plain top-level module instead.
    from load_cpiextract_data import LOOKUP_DTYPES

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def download_with_progress(url: str, dest_path: str, force: bool = False,
                            desc: str | None = None) -> str:
    """Streams url to dest_path with a progress bar. Skips if dest_path already exists
    and force=False, so re-running the script doesn't re-download multi-GB files.

    desc overrides the progress bar's label - pass this whenever dest_path is a fixed,
    generic local filename that doesn't reveal what's actually being fetched (e.g.
    ChEMBL is always saved to a fixed 'chembl.tar.gz' locally regardless of the
    upstream version, specifically so the extraction step downstream can find it
    predictably - without an explicit desc, the progress bar would show that generic
    name instead of which version is actually downloading)."""
    if os.path.exists(dest_path) and not force:
        label = desc if desc is not None else os.path.basename(dest_path)
        print(f"  [skip] {label} already present")
        return dest_path

    response = requests.head(url, allow_redirects=True)
    total_size = int(response.headers.get('content-length', 0))
    response = requests.get(url, stream=True)
    response.raise_for_status()

    label = desc if desc is not None else os.path.basename(dest_path)
    progress_bar = tqdm(total=total_size, unit='iB', unit_scale=True,
                         desc=f"Downloading {label}")
    tmp_path = dest_path + '.part'
    with open(tmp_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                progress_bar.update(len(chunk))
                f.write(chunk)
    progress_bar.close()
    os.replace(tmp_path, dest_path)
    return dest_path


def enforce_lookup_dtypes(df: pd.DataFrame, key: str) -> pd.DataFrame:
    """
    Casts each database's known exact-match lookup columns (see LOOKUP_DTYPES in
    load_cpiextract_data.py) to a clean, consistent representation before writing to
    CSV - mainly to avoid the on-disk text itself being misleading (e.g. a string ID
    that happens to look numeric getting written as "123" instead of "00123", losing
    the leading zeros before load_cpiextract_data.py's read-side dtype= enforcement
    even gets a chance to run). The read-side enforcement in load_dbs() is what
    actually guarantees the correct in-memory dtype; this is a belt-and-suspenders
    write-side pass so the CSV itself doesn't already have the information silently
    destroyed.
    """
    df = df.copy()
    for col, dt in LOOKUP_DTYPES.get(key, {}).items():
        if col not in df.columns:
            continue
        if dt == 'Int64':
            df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')
        else:
            # Only stringify non-null values - leaves genuine NaNs as NaN rather than
            # the literal string "nan", which would otherwise itself become a
            # silently-wrong match target downstream.
            df[col] = df[col].apply(lambda v: str(v) if pd.notna(v) else v)
    return df


def add_firstblock_and_link_to_pubchem(df: pd.DataFrame, external_id_col: str,
                                        pc_source_name: str, pc_sid: pd.DataFrame,
                                        pc_cid: pd.DataFrame) -> pd.DataFrame:
    """
    The one standardization step every non-PubChem database goes through: link this
    database's own external ID to a PubChem CID (via pc_sid, PubChem's own
    Source+ExternalID->CID mapping), then link that CID to PubChem's canonical InChIKey
    (via pc_cid), drop rows with no match on either step, and compute FirstBlock (the
    first 14-char block of the InChIKey) - exactly matching db_standardization.ipynb's
    per-database pattern (verified identical across BindingDB/ChEMBL/CTD/DrugBank/
    DrugCentral/DTC/OTP; STITCH differs only in how CID itself is derived, handled
    separately in prepare_stitch()).
    """
    pc_this_source = pc_sid[pc_sid['Source'] == pc_source_name].rename(
        columns={'External_ID': external_id_col})[[external_id_col, 'CID']]

    df = df.merge(pc_this_source, on=external_id_col, how='left')
    df['CID'] = df['CID'].fillna(-1).astype(int)
    df = df.merge(pc_cid[['CID', 'InChIKey']], on='CID', how='left')
    df = df.rename(columns={'InChIKey': 'inchikey'})

    total = len(df)
    df = df[df['CID'] != -1].reset_index(drop=True)
    df = df.dropna(subset=['inchikey']).reset_index(drop=True)
    print(f"    {total} interactions -> {len(df)} with a PubChem CID and inchikey match")

    df['FirstBlock'] = df['inchikey'].str.split('-').str[0]
    return df


# ---------------------------------------------------------------------------
# PubChem reference tables (built once, reused by every other database below)
# ---------------------------------------------------------------------------

def _iter_gzip_lines_with_progress(gz_path: str, desc: str):
    """
    Yields lines from a gzip-compressed file with a real, byte-based progress bar -
    tracking bytes read from the underlying COMPRESSED file (via tqdm.wrapattr) against
    its known total size. Files like SID-Map.gz and CID-InChI-Key.gz map every PubChem
    substance/compound record and can be several GB compressed; parsing them line-by-line
    in pure Python with zero progress feedback previously left no way to tell the
    difference between "still working" and "hung" (see conversation - this was a
    regression from the original notebook, which at least had a periodic print every 5M
    lines; this replaces that with an actual percentage-based progress bar instead).
    """
    total_size = os.path.getsize(gz_path)
    with open(gz_path, 'rb') as raw_f:
        with tqdm.wrapattr(raw_f, 'read', total=total_size, desc=desc,
                            unit='B', unit_scale=True) as wrapped_f:
            with gzip.open(wrapped_f, 'rt') as f:
                for line in f:
                    yield line


def build_pubchem_reference_tables(data_path: str, force_download: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Builds pc_sid (Source+External_ID -> CID) and pc_cid (CID -> InChIKey), PubChem's
    own cross-database identity mappings that every other database's standardization
    step below merges against. Caches both as CSVs so a re-run doesn't rebuild them.

    Reuses CID-InChI-Key.gz if prepare_pubchem() already downloaded it for the DuckDB
    build below, rather than fetching the same multi-GB file twice (the original
    notebooks downloaded it separately in each notebook).
    """
    print("Building PubChem reference tables (Source/CID/InChIKey mappings)...")
    pubchem_dir = os.path.join(data_path, 'pubchem')
    os.makedirs(pubchem_dir, exist_ok=True)

    sid_csv = os.path.join(data_path, 'pubchem-source-SIDs.csv')
    cid_csv = os.path.join(data_path, 'pubchem-cid-inchikey.csv')

    # ---- pc_sid: Source, External_ID, CID (only sources CPIExtract actually uses) ----
    if os.path.exists(sid_csv) and not force_download:
        print(f"  [skip] {os.path.basename(sid_csv)} already present")
        pc_sid = pd.read_csv(sid_csv, low_memory=False, index_col=0)
    else:
        sid_map_gz = os.path.join(pubchem_dir, 'SID-Map.gz')
        download_with_progress("https://ftp.ncbi.nlm.nih.gov/pubchem/Substance/Extras/SID-Map.gz",
                                sid_map_gz, force=force_download)

        src_list = ['BindingDB', 'ChEMBL', 'Comparative Toxicogenomics Database (CTD)',
                    'DrugBank', 'DrugCentral']
        data_4col = []
        for line in _iter_gzip_lines_with_progress(sid_map_gz, "Parsing SID-Map.gz"):
            parts = line.strip().split('\t')
            if len(parts) == 4:
                data_4col.append(parts)
        pc_sid = pd.DataFrame(data_4col, columns=['SID', 'Source', 'External_ID', 'CID'])
        pc_sid['SID'] = pc_sid['SID'].astype(int)
        pc_sid['CID'] = pc_sid['CID'].astype(int)
        pc_sid = pc_sid[pc_sid['Source'].isin(src_list)]
        pc_sid.to_csv(sid_csv)
        print(f"  Built {os.path.basename(sid_csv)}: {len(pc_sid):,} rows")

    # ---- pc_cid: CID -> InChIKey ----
    if os.path.exists(cid_csv) and not force_download:
        print(f"  [skip] {os.path.basename(cid_csv)} already present")
        pc_cid = pd.read_csv(cid_csv, low_memory=False, index_col=0)
    else:
        cid_inchikey_gz = os.path.join(pubchem_dir, 'CID-InChI-Key.gz')
        # Reuse the file if prepare_pubchem() already fetched it for the DuckDB build
        download_with_progress("https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz",
                                cid_inchikey_gz, force=force_download)

        data_2col = []
        for line in _iter_gzip_lines_with_progress(cid_inchikey_gz, "Parsing CID-InChI-Key.gz"):
            parts = line.strip().split('\t')
            data_2col.append([parts[0], parts[2]])
        pc_cid = pd.DataFrame(data_2col, columns=['CID', 'InChIKey'])
        pc_cid['CID'] = pc_cid['CID'].astype(int)
        pc_cid.to_csv(cid_csv)
        print(f"  Built {os.path.basename(cid_csv)}: {len(pc_cid):,} rows")

    return pc_sid, pc_cid


# ---------------------------------------------------------------------------
# Fully-automated databases
# ---------------------------------------------------------------------------

def _find_pubchem_target_file_url(filename_prefix: str) -> str:
    """
    PubChem's Target/ directory listing changes its exact filenames over time -
    confirmed live (see conversation) that 'gene2info.gz' (what the original notebook
    used) has since become 'gene2info.tsv.gz'. Rather than hardcoding a specific
    filename that can go stale the same way again, this scrapes the directory listing
    for whatever file currently starts with the given prefix (e.g. 'gene2info') -
    same technique already used for ChEMBL's release filename and OTP's version
    directory.
    """
    listing_url = "https://ftp.ncbi.nlm.nih.gov/pubchem/Target/"
    response = requests.get(listing_url)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, 'html.parser')

    candidates = [a['href'] for a in soup.find_all('a') if a['href'].startswith(filename_prefix)]
    if not candidates:
        raise RuntimeError(
            f"Could not find a file starting with '{filename_prefix}' at {listing_url} - "
            f"check that URL directly to confirm the current filename."
        )
    if len(candidates) > 1:
        print(f"  Note: multiple files match prefix '{filename_prefix}' at {listing_url} "
              f"({candidates}) - using the first one alphabetically: {sorted(candidates)[0]}")
    filename = sorted(candidates)[0]
    return listing_url + filename


def prepare_pubchem(data_path: str, force_download: bool = False) -> None:
    print("\n=== PubChem ===")
    pubchem_dir = os.path.join(data_path, 'pubchem')
    os.makedirs(pubchem_dir, exist_ok=True)

    bioact_file = os.path.join(pubchem_dir, 'pc_bioactivities.tsv.gz')
    gene2info_gz = os.path.join(pubchem_dir, 'gene2info.tsv.gz')
    gene2info_csv = os.path.join(pubchem_dir, 'gene2info.csv')
    cid_inchikey_file = os.path.join(pubchem_dir, 'CID-InChI-Key.gz')
    db_file = os.path.join(pubchem_dir, 'pubchem.duckdb')

    if os.path.exists(db_file) and not force_download:
        print(f"  [skip] pubchem.duckdb already present")
        return

    download_with_progress("https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/bioactivities.tsv.gz",
                            bioact_file, force=force_download)
    gene2info_url = _find_pubchem_target_file_url('gene2info')
    download_with_progress(gene2info_url, gene2info_gz, force=force_download,
                            desc=gene2info_url.rsplit('/', 1)[-1])
    download_with_progress("https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz",
                            cid_inchikey_file, force=force_download)

    gene2info = pd.read_csv(gene2info_gz, sep='\t', compression='gzip')
    gene2info.to_csv(gene2info_csv)

    print("  Building pubchem.duckdb (bioactivities + CID-InChIKey + gene_info)...")
    con = duckdb.connect(db_file)
    con.execute(f"""
        CREATE TABLE bioactivities AS
        SELECT * FROM read_csv('{bioact_file}', delim='\t', compression='gzip', header=true)
    """)
    con.execute(f"""
        CREATE TABLE cid_inchikey AS
        SELECT
            column0 AS CID,
            column1 AS InChI,
            column2 AS InChIKey,
            substr(column2, 1, 14) AS firstblock
        FROM read_csv('{cid_inchikey_file}', delim='\t', compression='gzip',
                       header=false, auto_detect=false,
                       columns={{'column0': 'INTEGER', 'column1': 'VARCHAR', 'column2': 'VARCHAR'}})
    """)
    con.execute(f"""
        CREATE TABLE gene_info AS SELECT * FROM read_csv('{gene2info_csv}', header=true)
    """)
    con.execute("CREATE INDEX idx_bioact_cid ON bioactivities(CID)")
    con.execute('CREATE INDEX idx_bioact_gene ON bioactivities("Gene ID")')
    con.execute("CREATE INDEX idx_cid_inchikey ON cid_inchikey(CID)")
    con.execute("CREATE INDEX idx_firstblock ON cid_inchikey(firstblock)")
    con.execute("CREATE INDEX idx_gene_id ON gene_info(GeneID)")
    con.execute("CREATE INDEX idx_gene_tax ON gene_info(TaxonomyID)")

    bioact_count = con.execute("SELECT COUNT(*) FROM bioactivities").fetchone()[0]
    cid_count = con.execute("SELECT COUNT(*) FROM cid_inchikey").fetchone()[0]
    gene_count = con.execute("SELECT COUNT(*) FROM gene_info").fetchone()[0]
    con.close()

    print(f"  Done: {bioact_count:,} bioactivities, {cid_count:,} CID-InChIKey rows, "
          f"{gene_count:,} genes -> {db_file}")


def _find_latest_chembl_sqlite_url() -> str:
    """
    ChEMBL's 'latest' directory always points at the current release, but the SQLite
    dump's filename itself embeds the version number (e.g. chembl_37_sqlite.tar.gz) -
    hardcoding that number, as the original notebook did, silently breaks the moment a
    new ChEMBL version ships (the file that version number points to either stops
    existing or, worse, could resolve to something stale). Scraping the directory
    listing for whatever chembl_NN_sqlite.tar.gz is actually present avoids hardcoding
    the version at all - same technique already used for OTP's parquet file listing.
    """
    listing_url = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
    response = requests.get(listing_url)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, 'html.parser')

    candidates = [a['href'] for a in soup.find_all('a')
                  if re.match(r'^chembl_\d+_sqlite\.tar\.gz$', a['href'])]
    if not candidates:
        raise RuntimeError(
            f"Could not find a chembl_NN_sqlite.tar.gz file at {listing_url} - "
            f"the naming convention or directory layout may have changed."
        )
    # Sort by the embedded version number (not lexicographically - "9" would
    # otherwise sort after "10" as plain text) and take the highest.
    candidates.sort(key=lambda name: int(re.search(r'\d+', name).group()))
    latest = candidates[-1]
    print(f"  Latest ChEMBL release found: {latest}")
    return listing_url + latest


def _find_chembl_release_sqlite_url(release: str) -> str:
    """
    Finds the SQLite dump for a SPECIFIC ChEMBL release number (e.g. release='33'),
    under ChEMBLdb/releases/ rather than ChEMBLdb/latest/ - same directory-scraping
    approach as _find_latest_chembl_sqlite_url(), just pointed at a specific version's
    subdirectory.

    Confirmed live (see conversation) that early releases use a zero-padded 2-digit
    directory/filename convention (chembl_01, not chembl_1) while later two-digit-plus
    releases (chembl_33, etc.) don't need padding at all. Rather than requiring the
    caller to know which convention a given release number uses, this tries the value
    as given first, then a zero-padded 2-digit form as a fallback for single-digit
    numeric releases.
    """
    candidates_to_try = [release]
    if release.isdigit() and len(release) < 2:
        candidates_to_try.append(release.zfill(2))

    errors = []
    for candidate in candidates_to_try:
        listing_url = f"https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{candidate}/"
        try:
            response = requests.get(listing_url)
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            errors.append(f"{listing_url} -> {e}")
            continue

        soup = BeautifulSoup(response.text, 'html.parser')
        expected_name = f"chembl_{candidate}_sqlite.tar.gz"
        matches = [a['href'] for a in soup.find_all('a') if a['href'] == expected_name]
        if matches:
            print(f"  Using requested ChEMBL release: {matches[0]}")
            return listing_url + matches[0]
        errors.append(f"{listing_url} -> directory exists but no {expected_name} found in it")

    raise RuntimeError(
        f"Could not find a ChEMBL release matching '{release}'. Tried:\n  " +
        "\n  ".join(errors) +
        f"\nCheck https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/ directly "
        f"to confirm the release number and directory naming."
    )



def prepare_chembl(data_path: str, force_download: bool = False,
                    chembl_release: str | None = None) -> None:
    """
    chembl_release=None (default) auto-detects and uses the latest ChEMBL release.
    Pass a specific release number (e.g. chembl_release='33') to fetch that release
    instead, from ChEMBLdb/releases/ rather than ChEMBLdb/latest/.
    """
    print("\n=== ChEMBL ===")
    out_csv = os.path.join(data_path, 'ChEMBL_standardized.csv')
    if os.path.exists(out_csv) and not force_download:
        print(f"  [skip] {os.path.basename(out_csv)} already present")
        return

    tar_path = os.path.join(data_path, 'chembl.tar.gz')
    if chembl_release is None:
        chembl_url = _find_latest_chembl_sqlite_url()
    else:
        chembl_url = _find_chembl_release_sqlite_url(chembl_release)
    download_with_progress(chembl_url, tar_path, force=force_download,
                            desc=chembl_url.rsplit('/', 1)[-1])

    print("  Extracting...")
    with tarfile.open(tar_path, 'r:gz') as tar:
        tar.extractall(data_path)

    db_matches = glob.glob(os.path.join(data_path, 'chembl_*', 'chembl_*_sqlite', 'chembl_*.db'))
    if not db_matches:
        raise FileNotFoundError("Could not locate the extracted ChEMBL SQLite file - "
                                 "the archive layout may have changed from what this script expects.")
    conn = sqlite3.connect(db_matches[0])

    query = """
    SELECT DISTINCT
        md.chembl_id AS molecule_chembl_id,
        md.pref_name AS molecule_name,
        cs.standard_inchi_key AS inchikey,
        td.chembl_id AS target_chembl_id,
        td.pref_name AS target_name,
        td.target_type,
        td.organism,
        cs_comp.accession AS uniprot,
        act.standard_relation,
        act.standard_type,
        act.standard_value,
        act.standard_units,
        act.pchembl_value,
        act.activity_comment,
        act.action_type,
        act.data_validity_comment,
        src.src_id
    FROM molecule_dictionary md
    JOIN compound_structures cs ON md.molregno = cs.molregno
    JOIN activities act ON md.molregno = act.molregno
    JOIN assays ass ON act.assay_id = ass.assay_id
    JOIN target_dictionary td ON ass.tid = td.tid
    LEFT JOIN target_components tc ON td.tid = tc.tid
    LEFT JOIN component_sequences cs_comp ON tc.component_id = cs_comp.component_id
    LEFT JOIN source src ON ass.src_id = src.src_id
    WHERE act.standard_value IS NOT NULL
    AND cs_comp.component_type = 'PROTEIN'
    AND cs_comp.tax_id = 9606
    AND td.target_type = 'SINGLE PROTEIN'
    AND td.organism = 'Homo sapiens'
    AND cs_comp.accession IS NOT NULL
    AND ass.variant_id IS NULL
    """
    print("  Running extraction query...")
    chembl = pd.read_sql_query(query, conn)
    conn.close()
    print(f"  Retrieved {len(chembl)} compound-protein interactions")

    print("  Standardizing (linking to PubChem CID/InChIKey)...")
    pc_sid, pc_cid = build_pubchem_reference_tables(data_path)
    chembl = chembl.rename(columns={'inchikey': 'db_inchikey'})
    chembl = add_firstblock_and_link_to_pubchem(chembl, 'molecule_chembl_id', 'ChEMBL', pc_sid, pc_cid)

    chembl = enforce_lookup_dtypes(chembl, 'chembl')
    chembl.to_csv(out_csv, index=False)
    print(f"  Saved {out_csv}")


def _find_latest_otp_version() -> str:
    """
    Open Targets Platform releases live in version-numbered subdirectories (e.g. 25.12,
    25.09) under a stable parent path - scraping that parent directory's listing for the
    highest YY.MM-formatted subdirectory avoids hardcoding a version that goes stale the
    moment a new release ships.
    """
    parent_url = 'http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/'
    response = requests.get(parent_url)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, 'html.parser')

    version_pattern = re.compile(r'^(\d{2})\.(\d{2})/?$')
    versions = []
    for a in soup.find_all('a'):
        match = version_pattern.match(a['href'])
        if match:
            versions.append((int(match.group(1)), int(match.group(2)), match.group(0).rstrip('/')))
    if not versions:
        raise RuntimeError(
            f"Could not find any YY.MM-formatted version directory at {parent_url} - "
            f"the naming convention or directory layout may have changed."
        )
    versions.sort()
    latest = versions[-1][2]
    print(f"  Latest Open Targets Platform version found: {latest}")
    return latest


def prepare_otp(data_path: str, force_download: bool = False,
                 otp_version: str | None = None) -> None:
    print("\n=== Open Targets Platform (OTP) ===")
    out_csv = os.path.join(data_path, 'OTP_standardized.csv')
    if os.path.exists(out_csv) and not force_download:
        print(f"  [skip] {os.path.basename(out_csv)} already present")
        return

    otp_dir = os.path.join(data_path, 'otp_raw')
    os.makedirs(otp_dir, exist_ok=True)

    if otp_version is None:
        otp_version = _find_latest_otp_version()
    else:
        print(f"  Using explicitly-requested OTP version: {otp_version}")

    url = f'http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{otp_version}/output/drug_mechanism_of_action/'
    print(f"  Listing parquet files at {url}")
    response = requests.get(url)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, 'html.parser')
    parquet_links = [a['href'] for a in soup.find_all('a') if a['href'].endswith('.parquet')]
    print(f"  Found {len(parquet_links)} parquet files")

    for filename in tqdm(parquet_links, desc="Downloading OTP parquet files"):
        download_with_progress(url + filename, os.path.join(otp_dir, filename), force=force_download)

    parquet_files = glob.glob(os.path.join(otp_dir, '*.parquet'))
    OTP_data = pd.concat([pd.read_parquet(f) for f in parquet_files], ignore_index=True)
    print(f"  Retrieved {len(OTP_data)} rows")

    def parse_list_col(x):
        if isinstance(x, str):
            cleaned = x.strip("[]'\"").replace("'", "").replace('"', '')
            return [item.strip() for item in cleaned.split() if item.strip()]
        return x

    OTP_data['chemblIds'] = OTP_data['chemblIds'].apply(parse_list_col)
    OTP_data = OTP_data.explode('chemblIds').reset_index(drop=True)
    OTP_data['targets'] = OTP_data['targets'].apply(parse_list_col)
    OTP_data = OTP_data.explode('targets').reset_index(drop=True)
    OTP_data = OTP_data.dropna(subset=['targets', 'chemblIds']).reset_index(drop=True)
    print(f"  After exploding chemblIds/targets: {len(OTP_data)} rows")

    print("  Standardizing (linking to PubChem CID/InChIKey)...")
    pc_sid, pc_cid = build_pubchem_reference_tables(data_path)
    OTP_data = add_firstblock_and_link_to_pubchem(OTP_data, 'chemblIds', 'ChEMBL', pc_sid, pc_cid)

    OTP_data = enforce_lookup_dtypes(OTP_data, 'otp')
    OTP_data.to_csv(out_csv, index=False)
    print(f"  Saved {out_csv}")


# ---------------------------------------------------------------------------
# Manual-download-dependent databases
# ---------------------------------------------------------------------------

# One entry per manual source: a glob pattern for the expected raw file (some
# filenames embed a date/version, e.g. BindingDB's monthly release) and instructions
# printed verbatim if no matching file is found in data_path.
MANUAL_SOURCES = {
    'bindingdb': {
        'glob': 'BindingDB_All_*_tsv.zip',
        'instructions': (
            "BindingDB: download the 'BindingDB_All' TSV archive from\n"
            "    https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp\n"
            "and place the .zip (e.g. BindingDB_All_202512_tsv.zip) directly in "
            "--data-path."
        ),
    },
    'ctd': {
        'glob': 'CTD_chem_gene_ixns.csv.gz',
        'instructions': (
            "CTD: download 'CTD_chem_gene_ixns.csv.gz' from\n"
            "    https://ctdbase.org/downloads/#cg\n"
            "and place it directly in --data-path."
        ),
    },
    'drugbank': {
        'glob': 'drugbank_all_full_database.xml.zip',
        'instructions': (
            "DrugBank: requires an academic license. Download "
            "'drugbank_all_full_database.xml.zip' from\n"
            "    https://go.drugbank.com/releases/latest\n"
            "and place it directly in --data-path."
        ),
    },
    'stitch': {
        'glob': '9606.protein_chemical.links.detailed.v*.tsv.gz',
        'instructions': (
            "STITCH: download the human-only ('9606') detailed protein-chemical "
            "links file from\n"
            "    http://stitch.embl.de/cgi/download.pl\n"
            "and place it directly in --data-path."
        ),
    },
    'dtc': {
        'glob': 'DTC.csv',
        'instructions': (
            "DTC (DrugTargetCommons): download the full activities dataset from\n"
            "    https://drugtargetcommons.fimm.fi/\n"
            "and place it as DTC.csv directly in --data-path."
        ),
    },
    'drugcentral': {
        'glob': 'drug.target.interaction.tsv.gz',
        'instructions': (
            "DrugCentral: download 'drug.target.interaction.tsv.gz' from\n"
            "    https://drugcentral.org/download\n"
            "and place it directly in --data-path."
        ),
    },
}


def _find_manual_file(data_path: str, key: str) -> str | None:
    matches = glob.glob(os.path.join(data_path, MANUAL_SOURCES[key]['glob']))
    return matches[0] if matches else None


def _print_missing_instructions(key: str) -> None:
    print(f"  [MISSING] {MANUAL_SOURCES[key]['instructions']}")
    print(f"  Skipping {key} - re-run this script after placing the file.")


def prepare_bindingdb(data_path: str) -> None:
    print("\n=== BindingDB ===")
    out_csv = os.path.join(data_path, 'BindingDB_standardized.csv')
    raw_zip = _find_manual_file(data_path, 'bindingdb')
    if raw_zip is None:
        _print_missing_instructions('bindingdb')
        return

    print(f"  Extracting {os.path.basename(raw_zip)}...")
    with zipfile.ZipFile(raw_zip, 'r') as zip_ref:
        zip_ref.extractall(data_path)

    tsv_path = os.path.join(data_path, 'BindingDB_All.tsv')
    print("  Loading (this is a large file, can take a while)...")
    BDB_data = pd.read_csv(tsv_path, sep='\t', low_memory=False, usecols=[
        'Ligand SMILES', 'Ligand InChI', 'BindingDB MonomerID', 'Ligand InChI Key',
        'BindingDB Ligand Name', 'Target Name',
        'Target Source Organism According to Curator or DataSource',
        'Ki (nM)', 'IC50 (nM)', 'Kd (nM)', 'EC50 (nM)', 'pH', 'Temp (C)',
        'Curation/DataSource',
        'Number of Protein Chains in Target (>1 implies a multichain complex)',
        'UniProt (SwissProt) Entry Name of Target Chain 1',
        'UniProt (SwissProt) Primary ID of Target Chain 1',
    ])
    BDB_data = BDB_data.rename(columns={'Ligand InChI Key': 'inchikey', 'BindingDB MonomerID': 'BindingDB MonomerID'})
    print(f"  Retrieved {len(BDB_data)} compound-protein interactions")

    print("  Standardizing (linking to PubChem CID/InChIKey)...")
    pc_sid, pc_cid = build_pubchem_reference_tables(data_path)
    BDB_data = BDB_data.rename(columns={'inchikey': 'old_inchikey'})
    BDB_data = add_firstblock_and_link_to_pubchem(BDB_data, 'BindingDB MonomerID', 'BindingDB', pc_sid, pc_cid)

    BDB_data = enforce_lookup_dtypes(BDB_data, 'bdb')
    BDB_data.to_csv(out_csv, index=False)
    print(f"  Saved {out_csv}")


def prepare_ctd(data_path: str) -> None:
    print("\n=== CTD ===")
    out_csv = os.path.join(data_path, 'CTD_standardized.csv')
    raw_gz = _find_manual_file(data_path, 'ctd')
    if raw_gz is None:
        _print_missing_instructions('ctd')
        return

    csv_path = os.path.join(data_path, 'CTD_chem_gene_ixns.csv')
    with gzip.open(raw_gz, 'rb') as f_in, open(csv_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    # CTD's CSV export has a fixed 28-line preamble (rows 0-26 and 28 are metadata,
    # not data - matches db_download_2025.ipynb exactly).
    rows_to_skip = list(range(27)) + [28]
    CTD_data = pd.read_csv(csv_path, skiprows=rows_to_skip, header=0)
    CTD_data = CTD_data.rename(columns={'# ChemicalName': 'ChemicalName'})
    print(f"  Retrieved {len(CTD_data)} compound-protein interactions")

    print("  Standardizing (linking to PubChem CID/InChIKey)...")
    pc_sid, pc_cid = build_pubchem_reference_tables(data_path)
    CTD_data = add_firstblock_and_link_to_pubchem(
        CTD_data, 'ChemicalID', 'Comparative Toxicogenomics Database (CTD)', pc_sid, pc_cid)

    CTD_data = enforce_lookup_dtypes(CTD_data, 'ctd')
    CTD_data.to_csv(out_csv, index=False)
    print(f"  Saved {out_csv}")


def _extract_drugbank_protein_info(protein_element, ns) -> dict:
    protein_data = {}
    p_id = protein_element.find('db:id', ns)
    p_name = protein_element.find('db:name', ns)
    p_organism = protein_element.find('db:organism', ns)
    protein_data['protein_id'] = p_id.text if p_id is not None else ''
    protein_data['protein_name'] = p_name.text if p_name is not None else ''
    protein_data['organism'] = p_organism.text if p_organism is not None else ''

    external_ids = {}
    polypeptide = protein_element.find('db:polypeptide', ns)
    if polypeptide is not None:
        for ext_id in polypeptide.findall('.//db:external-identifiers/db:external-identifier', ns):
            resource = ext_id.find('db:resource', ns)
            identifier = ext_id.find('db:identifier', ns)
            if resource is not None and identifier is not None:
                external_ids[resource.text] = identifier.text
    protein_data['external_ids'] = external_ids
    return protein_data


def _parse_drugbank_xml(root) -> list[dict]:
    ns = {'db': 'http://www.drugbank.ca'}
    all_rows = []
    for drug in root.findall('db:drug', ns):
        primary_id = drug.find("db:drugbank-id[@primary='true']", ns)
        drugbank_id = primary_id.text if primary_id is not None else ''
        name_elem = drug.find('db:name', ns)
        drug_name = name_elem.text if name_elem is not None else ''

        synonyms = drug.findall('.//db:synonyms/db:synonym', ns)
        synonyms_str = ' | '.join(syn.text for syn in synonyms if syn.text)

        pubchem_cid = ''
        for ext_id in drug.findall('.//db:external-identifiers/db:external-identifier', ns):
            resource = ext_id.find('db:resource', ns)
            identifier = ext_id.find('db:identifier', ns)
            if resource is not None and identifier is not None and 'PubChem' in resource.text and 'Compound' in resource.text:
                pubchem_cid = identifier.text
                break

        inchikey = ''
        for prop in drug.findall('.//db:calculated-properties/db:property', ns):
            kind = prop.find('db:kind', ns)
            value = prop.find('db:value', ns)
            if kind is not None and value is not None and kind.text == 'InChIKey':
                inchikey = value.text
                break

        for protein_type, xpath in [('target', './/db:targets/db:target'),
                                     ('enzyme', './/db:enzymes/db:enzyme'),
                                     ('carrier', './/db:carriers/db:carrier'),
                                     ('transporter', './/db:transporters/db:transporter')]:
            for element in drug.findall(xpath, ns):
                protein_info = _extract_drugbank_protein_info(element, ns)
                row = {
                    'drugbank_id': drugbank_id, 'name': drug_name, 'synonyms': synonyms_str,
                    'pubchem_cid': pubchem_cid, 'inchikey': inchikey,
                    'protein': protein_info['protein_id'], 'protein_type': protein_type,
                    'protein_name': protein_info['protein_name'], 'organism': protein_info['organism'],
                }
                row.update(protein_info['external_ids'])
                all_rows.append(row)
    return all_rows


def prepare_drugbank(data_path: str) -> None:
    print("\n=== DrugBank ===")
    out_csv = os.path.join(data_path, 'DrugBank_standardized.csv')
    raw_zip = _find_manual_file(data_path, 'drugbank')
    if raw_zip is None:
        _print_missing_instructions('drugbank')
        return

    print(f"  Extracting {os.path.basename(raw_zip)}...")
    with zipfile.ZipFile(raw_zip, 'r') as zip_ref:
        zip_ref.extractall(data_path)

    xml_candidates = glob.glob(os.path.join(data_path, '*full database*.xml'))
    if not xml_candidates:
        raise FileNotFoundError("Could not locate the extracted DrugBank XML file.")
    print("  Parsing XML (this can take several minutes for the full database)...")
    tree = ET.parse(xml_candidates[0])
    root = tree.getroot()
    rows = _parse_drugbank_xml(root)

    DB_data = pd.DataFrame(rows)
    base_columns = ['drugbank_id', 'name', 'synonyms', 'pubchem_cid', 'inchikey',
                     'protein', 'protein_type', 'protein_name', 'organism']
    if len(DB_data) > 0:
        external_id_columns = sorted(c for c in DB_data.columns if c not in base_columns)
        DB_data = DB_data[base_columns + external_id_columns]
    print(f"  Extracted {len(DB_data)} drug-protein relationships")

    print("  Standardizing (linking to PubChem CID/InChIKey)...")
    pc_sid, pc_cid = build_pubchem_reference_tables(data_path)
    DB_data = DB_data.rename(columns={'inchikey': 'db_inchikey'})
    DB_data = add_firstblock_and_link_to_pubchem(DB_data, 'drugbank_id', 'DrugBank', pc_sid, pc_cid)

    DB_data = enforce_lookup_dtypes(DB_data, 'db')
    DB_data.to_csv(out_csv, index=False)
    print(f"  Saved {out_csv}")


def prepare_stitch(data_path: str) -> None:
    print("\n=== STITCH ===")
    out_csv = os.path.join(data_path, 'STITCH_standardized.csv')
    raw_gz = _find_manual_file(data_path, 'stitch')
    if raw_gz is None:
        _print_missing_instructions('stitch')
        return

    tsv_path = raw_gz[:-3]  # strip .gz
    with gzip.open(raw_gz, 'rb') as f_in, open(tsv_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    stch_data = pd.read_csv(tsv_path, sep='\t')
    print(f"  Retrieved {len(stch_data)} compound-protein interactions")

    def extract_cid(chemical_id):
        match = re.search(r'CID[a-z]*(\d+)', chemical_id)
        return int(match.group(1)) if match else None

    stch_data['CID'] = stch_data['chemical'].apply(extract_cid)

    print("  Standardizing (linking to PubChem InChIKey)...")
    _, pc_cid = build_pubchem_reference_tables(data_path)
    stch_data = stch_data.merge(pc_cid[['CID', 'InChIKey']], on='CID', how='left')
    stch_data = stch_data.rename(columns={'InChIKey': 'inchikey'})
    total = len(stch_data)
    stch_data = stch_data[stch_data['CID'].notna()].reset_index(drop=True)
    stch_data = stch_data.dropna(subset=['inchikey']).reset_index(drop=True)
    print(f"    {total} interactions -> {len(stch_data)} with a PubChem CID and inchikey match")
    stch_data['FirstBlock'] = stch_data['inchikey'].str.split('-').str[0]

    stch_data = enforce_lookup_dtypes(stch_data, 'stitch')
    stch_data.to_csv(out_csv, index=False)
    print(f"  Saved {out_csv}")


def prepare_dtc(data_path: str) -> None:
    print("\n=== DTC (DrugTargetCommons) ===")
    out_csv = os.path.join(data_path, 'DTC_standardized.csv')
    raw_csv = _find_manual_file(data_path, 'dtc')
    if raw_csv is None:
        _print_missing_instructions('dtc')
        return

    DTC_data = pd.read_csv(raw_csv, sep=',')
    DTC_data = DTC_data.dropna(subset=['compound_id'])
    print(f"  Retrieved {len(DTC_data)} compound-protein interactions")

    print("  Standardizing (linking to PubChem CID/InChIKey via ChEMBL compound_id)...")
    pc_sid, pc_cid = build_pubchem_reference_tables(data_path)
    if 'CID' in DTC_data.columns:
        DTC_data = DTC_data.drop(columns=['CID'])
    DTC_data = add_firstblock_and_link_to_pubchem(DTC_data, 'compound_id', 'ChEMBL', pc_sid, pc_cid)

    DTC_data = enforce_lookup_dtypes(DTC_data, 'dtc')
    DTC_data.to_csv(out_csv, index=False)
    print(f"  Saved {out_csv}")


def prepare_drugcentral(data_path: str) -> None:
    print("\n=== DrugCentral ===")
    out_csv = os.path.join(data_path, 'DrugCentral_standardized.csv')
    raw_gz = _find_manual_file(data_path, 'drugcentral')
    if raw_gz is None:
        _print_missing_instructions('drugcentral')
        return

    DC_data = pd.read_csv(raw_gz, sep='\t', compression='gzip', low_memory=False)
    DC_data = DC_data[['DRUG_NAME', 'STRUCT_ID', 'TARGET_NAME', 'TARGET_CLASS', 'ACCESSION',
                        'GENE', 'SWISSPROT', 'ACT_VALUE', 'ACT_UNIT', 'ACT_TYPE', 'ACT_COMMENT',
                        'ACT_SOURCE', 'RELATION', 'ACTION_TYPE', 'ORGANISM']]
    print(f"  Retrieved {len(DC_data)} compound-protein interactions")

    print("  Standardizing (linking to PubChem CID/InChIKey)...")
    pc_sid, pc_cid = build_pubchem_reference_tables(data_path)
    DC_data['STRUCT_ID'] = DC_data['STRUCT_ID'].astype(str)
    DC_data = add_firstblock_and_link_to_pubchem(DC_data, 'STRUCT_ID', 'DrugCentral', pc_sid, pc_cid)

    DC_data = enforce_lookup_dtypes(DC_data, 'dc')
    DC_data.to_csv(out_csv, index=False)
    print(f"  Saved {out_csv}")


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

ALL_DATABASES = ['pubchem', 'chembl', 'otp', 'bindingdb', 'ctd', 'drugbank',
                  'stitch', 'dtc', 'drugcentral']

AUTOMATED = {
    'pubchem': prepare_pubchem,
    'chembl': prepare_chembl,
    'otp': prepare_otp,
}

MANUAL = {
    'bindingdb': prepare_bindingdb,
    'ctd': prepare_ctd,
    'drugbank': prepare_drugbank,
    'stitch': prepare_stitch,
    'dtc': prepare_dtc,
    'drugcentral': prepare_drugcentral,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data-path', default=os.path.join(os.path.expanduser('~'), 'CPE_data'),
                         help="Output folder (default: ~/CPE_data)")
    parser.add_argument('--only', nargs='*', choices=ALL_DATABASES, default=None,
                         help="Only prepare these databases")
    parser.add_argument('--skip', nargs='*', choices=ALL_DATABASES, default=[],
                         help="Skip these databases")
    parser.add_argument('--force-download', action='store_true',
                         help="Re-download and re-build even if output files already exist")
    parser.add_argument('--otp-version', default=None,
                         help="Open Targets Platform release version to fetch "
                              "(default: auto-detect the latest available)")
    parser.add_argument('--chembl-release', default=None,
                         help="Specific ChEMBL release number to fetch, e.g. '33' "
                              "(default: auto-detect the latest available)")
    args = parser.parse_args()

    os.makedirs(args.data_path, exist_ok=True)
    targets = args.only if args.only is not None else ALL_DATABASES
    targets = [t for t in targets if t not in args.skip]

    print(f"Preparing CPIExtract data in: {args.data_path}")
    print(f"Databases: {', '.join(targets)}\n")

    for key in targets:
        if key in AUTOMATED:
            if key == 'otp':
                AUTOMATED[key](args.data_path, force_download=args.force_download,
                                otp_version=args.otp_version)
            elif key == 'chembl':
                AUTOMATED[key](args.data_path, force_download=args.force_download,
                                chembl_release=args.chembl_release)
            else:
                AUTOMATED[key](args.data_path, force_download=args.force_download)
        elif key in MANUAL:
            MANUAL[key](args.data_path)

    print("\nDone. Load this folder into the pipeline with load_cpiextract_data.py.")


if __name__ == '__main__':
    main()