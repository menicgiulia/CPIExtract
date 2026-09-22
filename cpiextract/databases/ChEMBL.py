'''Loading,searching,filtering and preprocessing data from ChEMBL.'''

import pandas as pd
import numpy as np
import re

from ..utils.typing import Connection
from ..servers.PubChemServer import PubChemServer
from ..servers.BiomartServer import BiomartServer
from ..servers.MyGeneServer import MyGeneServer
from .Database import Database
from ..data_manager import *

class ChEMBL(Database):
    '''Loading,searching,filtering and preprocessing data from ChEMBL.'''

    # These match the native chembl_webresource_client field names, so the local/SQL path (the
    # SQL query supplied to LocalManager/SQLManager) must alias to these exact names.
    ACTIVITY_FIELDS = ['molecule_chembl_id', 'inchikey', 'activity_comment', 'data_validity_comment',
                        'organism', 'src_id', 'pchembl_value', 'target_chembl_id',
                        'standard_type', 'standard_relation', 'standard_value', 'standard_units']


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
        elif connection is not None:
            self.data_manager = SQLManager(connection, 'CHEMBL')
        else:
            raise ValueError(
                "ChEMBL requires either a local 'database' DataFrame (a bulk ChEMBL dump) or a "
                "SQL 'connection' - live API-only mode is no longer supported."
            )

    # Recognized molar concentration units within ChEMBL and their multiplier to nM.
    UNIT_TO_NM = {'M': 1e9, 'mM': 1e6, 'uM': 1e3, 'µM': 1e3, 'microM': 1e3, 
                  'nM': 1, 'pM': 1e-3, 'fM': 1e-6}

    # Mass/volume concentration units (needs a molecular_weight to finish the conversion to molar).
    MASS_VOLUME_TO_G_PER_L = {'ug.mL-1': 1e-3, 'ug ml-1': 1e-3, 'ug/mL': 1e-3, 'ug/ml': 1e-3,
        'ug l-1': 1e-6, 'ug/L': 1e-6, 'ug/l': 1e-6,'mg/L': 1e-3, 'mg/l': 1e-3,
        'mg/dl': 1e-2, 'mg dl-1': 1e-2,'mg/ml': 1, 'mg.mL-1': 1,
        'ng/ml': 1e-6, 'ng ml-1': 1e-6, 'ng/mL': 1e-6,'pg/ml': 1e-9, 'pg/mL': 1e-9}

    # Association-constant units (Ka, roughly 1/Kd) 
    ASSOC_UNIT_TO_INV_M = {'L/mol': 1, 'l M-1': 1, 'L mol-1': 1, 'l mol-1': 1, 'M-1': 1, 
                           '1/M': 1, '/M': 1}

    # Units that are ALREADY on a pChEMBL-equivalent (-log10) scale
    P_SCALE_UNITS = {'pKi', 'pKd', 'pIC50', 'pEC50', 'pA2', 'pKB'}

    def _standard_value_to_nM(self, standard_units: str) -> float:
        """
        Returns the multiplier that converts a standard_value reported in
        `standard_units` into nM, or NaN if the unit isn't a convertible molar unit.
        """
        if pd.isnull(standard_units):
            return np.nan

        units = standard_units.strip()

        if units in self.UNIT_TO_NM:
            return self.UNIT_TO_NM[units]

        # older ASCII-apostrophe form '10'5nM' included
        match = re.match(r"^10[\^'](-?\d+)\s*(M|mM|uM|µM|microM|nM|pM|fM)$", units)
        if match:
            exponent, base_unit = match.groups()
            return (10 ** int(exponent)) * self.UNIT_TO_NM[base_unit]

        return np.nan

    def _assoc_to_invM(self, standard_units: str) -> float:
        """
        Returns the multiplier that converts a standard_value reported in `standard_units`
        into M^-1 (Ka), or NaN if the unit isn't a recognized plain association constant.
        """
        if pd.isnull(standard_units):
            return np.nan

        units = standard_units.strip()

        if units in self.ASSOC_UNIT_TO_INV_M:
            return self.ASSOC_UNIT_TO_INV_M[units]

        match = re.match(r"^10[\^'](-?\d+)\s*(L/mol|M-1|1/M|/M)$", units)
        if match:
            exponent, base_unit = match.groups()
            return (10 ** int(exponent)) * self.ASSOC_UNIT_TO_INV_M[base_unit]

        return np.nan

    def _filter_database_pre(self, chembl_raw: pd.DataFrame, pChEMBL_thres: float) -> pd.DataFrame:
        # Records with no evidence either way -> excluded entirely
        invalid_activities = ['Not Determined', 'Not Evaluated', 'inconclusive', 'undetermined', 'No data']
        # Records that are True Negatives (tested, no interaction found)
        negative_activities = ['Not Active', 'inactive']

        # ChEMBL only computes pchembl_value when data_validity_comment is NULL or exactly 'Manually validated'
        valid_validity = (chembl_raw['data_validity_comment'].isnull() |
            (chembl_raw['data_validity_comment'] == 'Manually validated'))

        chembl_act = chembl_raw.loc[
                    # Remove all non-human
                    (chembl_raw['organism'] == 'Homo sapiens') &  
                    # Remove various activity comments
                    (~chembl_raw['activity_comment'].isin(invalid_activities)) &  
                    # Only use interactions extracted from scientific literature
                    (chembl_raw['src_id'] == 1) &  
                    # Only use interactions with no (or manually-cleared) data validity comments
                    valid_validity
                ].drop_duplicates(subset=['activity_comment', 'molecule_chembl_id', 'pchembl_value', 'standard_value']).copy()

        # Clean/coerce numeric fields
        chembl_act['pchembl_value'] = chembl_act['pchembl_value'].apply(lambda x: re.sub(r'[^0-9.]', '', x) if type(x) is str else x)
        chembl_act['pchembl_value'] = chembl_act['pchembl_value'].apply(pd.to_numeric, errors='coerce')
        chembl_act['standard_value'] = pd.to_numeric(chembl_act['standard_value'], errors='coerce')

        # Comment-based true negatives
        comment_negative = chembl_act['activity_comment'].isin(negative_activities)

        # Censored (inequality) evidence: kept OUT of pchembl_value/ave_pChEMBL/std_pChEMBL
        no_pchembl = chembl_act['pchembl_value'].isnull()

        nM_multiplier = chembl_act['standard_units'].apply(self._standard_value_to_nM)
        concentration_nM = chembl_act['standard_value'] * nM_multiplier
        concentration_has_value = concentration_nM.notnull() & (concentration_nM > 0)
        concentration_derived = -np.log10(concentration_nM.where(concentration_has_value) * 1e-9)

        assoc_multiplier = chembl_act['standard_units'].apply(self._assoc_to_invM)
        Ka_M = chembl_act['standard_value'] * assoc_multiplier
        assoc_has_value = Ka_M.notnull() & (Ka_M > 0)
        assoc_derived = np.log10(Ka_M.where(assoc_has_value))

        pscale_has_value = chembl_act['standard_units'].isin(self.P_SCALE_UNITS) & chembl_act['standard_value'].notnull()
        pscale_derived = chembl_act['standard_value'].where(pscale_has_value)

        # 'same_' = same polarity as pChEMBL itself
        same_has_value = assoc_has_value | pscale_has_value
        same_derived = assoc_derived.where(assoc_has_value, pscale_derived)

        has_value = concentration_has_value | same_has_value
        derived_value = concentration_derived.where(concentration_has_value, same_derived)
        # True where this row's derived value came from a concentration-type
        # source; False for association/p-scale sources.
        inverted = concentration_has_value

        is_gt = chembl_act['standard_relation'].isin(['>', '>='])
        is_lt = chembl_act['standard_relation'].isin(['<', '<='])
        is_eq = chembl_act['standard_relation'] == '='

        # relation '>' on an inverted-polarity source, or '<' on a same-polarity source ->
        # true pChEMBL < derived value (upper bound)
        upper_bound = no_pchembl & has_value & ((inverted & is_gt) | (~inverted & is_lt))
        # relation '<' on an inverted-polarity source, or '>' on a same-polarity source ->
        # true pChEMBL > derived value (lower bound)
        lower_bound = no_pchembl & has_value & ((inverted & is_lt) | (~inverted & is_gt))
        # relation '=' with a derivable value but no ChEMBL-provided pchembl_value -> exact fill-in
        fill_eq = no_pchembl & has_value & is_eq

        chembl_act['pChEMBL_lt'] = np.nan
        chembl_act['pChEMBL_gt'] = np.nan
        chembl_act.loc[upper_bound, 'pChEMBL_lt'] = derived_value.loc[upper_bound]
        chembl_act.loc[lower_bound, 'pChEMBL_gt'] = derived_value.loc[lower_bound]
        chembl_act.loc[fill_eq, 'pchembl_value'] = derived_value.loc[fill_eq]

        # A '>'-censored row is negative-leaning evidence regardless of threshold (it's
        # already known to fall below the active range). A '<'-censored row is only kept
        # as positive evidence if its lower bound already clears the threshold.
        true_negatives = comment_negative | (upper_bound & (chembl_act['pChEMBL_lt'] <= pChEMBL_thres))
        censored_positive = lower_bound & (chembl_act['pChEMBL_gt'] > pChEMBL_thres)

        # Rows still unresolved (no exact pchembl_value, no molar/assoc/p-scale-derived value, not a comment-negative)
        is_mass_volume_unit = chembl_act['standard_units'].isin(self.MASS_VOLUME_TO_G_PER_L)
        pending_mw = chembl_act['pchembl_value'].isnull() & is_mass_volume_unit

        chembl_act['_pending_mw'] = pending_mw

        chembl_act = chembl_act.loc[chembl_act['pchembl_value'].notnull() | 
                                    true_negatives | censored_positive | pending_mw].copy()

        return chembl_act.reset_index(drop=True)

    def _filter_database_post(self, chembl_act: pd.DataFrame, pChEMBL_thres: float, molecular_weight: float|pd.Series|None=None) -> pd.DataFrame:
        """
        Second phase: resolves any row flagged '_pending_mw' by _filter_database_pre (mass/volume
        concentration units with no other resolvable evidence) now that a molecular_weight is
        available, applies the final comment-negative fallback (0 when nothing could be derived).
        """
        chembl_act = chembl_act.copy()
        if '_pending_mw' not in chembl_act.columns:
            chembl_act['_pending_mw'] = False
        pending = chembl_act['_pending_mw'].fillna(False)

        if molecular_weight is not None and pending.any():
            mass_g_per_L = chembl_act['standard_units'].map(self.MASS_VOLUME_TO_G_PER_L)
            mass_nM = chembl_act['standard_value'] * mass_g_per_L / molecular_weight * 1e9
            mass_has_value = pending & mass_nM.notnull() & (mass_nM > 0)
            mass_derived = -np.log10(mass_nM.where(mass_has_value) * 1e-9)

            is_gt = chembl_act['standard_relation'].isin(['>', '>='])
            is_lt = chembl_act['standard_relation'].isin(['<', '<='])
            is_eq = chembl_act['standard_relation'] == '='

            upper_bound = pending & mass_has_value & is_gt
            lower_bound = pending & mass_has_value & is_lt
            fill_eq = pending & mass_has_value & is_eq

            chembl_act.loc[upper_bound, 'pChEMBL_lt'] = mass_derived.loc[upper_bound]
            chembl_act.loc[lower_bound, 'pChEMBL_gt'] = mass_derived.loc[lower_bound]
            chembl_act.loc[fill_eq, 'pchembl_value'] = mass_derived.loc[fill_eq]

        negative_activities = ['Not Active', 'inactive']
        comment_negative = chembl_act['activity_comment'].isin(negative_activities)

        # Decide keep/drop only for rows that were left pending
        keep_pending = (chembl_act['pchembl_value'].notnull() |
            comment_negative |
            (chembl_act['pChEMBL_lt'].notnull() & (chembl_act['pChEMBL_lt'] <= pChEMBL_thres)) |
            (chembl_act['pChEMBL_gt'].notnull() & (chembl_act['pChEMBL_gt'] > pChEMBL_thres))
        )
        chembl_act = chembl_act.loc[~(pending.reindex(chembl_act.index, fill_value=False) & ~keep_pending)].copy()


        # Bare true negatives with no numeric evidence at all get an explicit 0 in pchembl_value.
        missing_val_negative = (
            chembl_act['activity_comment'].isin(negative_activities) &
            chembl_act['pchembl_value'].isnull() &
            chembl_act['pChEMBL_lt'].isnull()
        )
        chembl_act.loc[missing_val_negative, 'pchembl_value'] = 0

        return chembl_act.drop(columns=['_pending_mw']).reset_index(drop=True)

    def _filter_database(self, chembl_raw: pd.DataFrame, pChEMBL_thres: float, molecular_weight: float|pd.Series|None=None) -> pd.DataFrame:
        """Convenience wrapper for the single-compound-input path, where molecular_weight (if
        any) is already known before filtering."""
        chembl_act = self._filter_database_pre(chembl_raw, pChEMBL_thres)
        if len(chembl_act) == 0:
            return chembl_act.drop(columns=['_pending_mw'], errors='ignore')
        return self._filter_database_post(chembl_act, pChEMBL_thres, molecular_weight=molecular_weight)

    def _restructure_output(self, chembl_act: pd.DataFrame, compound_id_cols: list, protein_id_cols: list) -> pd.DataFrame:
        """
        Restructures filtered, row-per-activity-record ChEMBL data into the shared
        row-per-passed-interaction shape used across CPIExtract databases:

        compound_identifier cols | protein_identifier cols | standard_units
        | pChEMBL_eq | pChEMBL_lt | pChEMBL_gt | datasource

        pChEMBL_eq holds exact '=' measurements (or 0 for a bare comment-only negative
        with no derivable relation). pChEMBL_lt/pChEMBL_gt hold censored bounds. 
        At most one of the three is populated on any given row.

        Parameters
        ----------
        chembl_act : DataFrame
            Filtered, row-per-activity-record ChEMBL data (output of _filter_database,
            enriched with whatever identifier columns are available).
        compound_id_cols : list
            Columns that uniquely identify a compound (e.g. ['inchikey'] or ['CID']).
            Only columns actually present in chembl_act are used.
        protein_id_cols : list
            Columns that uniquely identify a protein (e.g. ['uniprot']).
            Only columns actually present in chembl_act are used.

        Returns
        -------
        DataFrame
            One row per activity record that passed filtering.
        """

        compound_id_cols = [c for c in compound_id_cols if c in chembl_act.columns]
        protein_id_cols = [c for c in protein_id_cols if c in chembl_act.columns]

        if len(chembl_act) == 0:
            return pd.DataFrame()

        out = chembl_act.rename(columns={'pchembl_value': 'pChEMBL_eq'}).copy()
        out['datasource'] = 'ChEMBL'

        keep = compound_id_cols + protein_id_cols + ['standard_type', 'pChEMBL_eq', 'pChEMBL_lt', 'pChEMBL_gt', 'datasource']
        keep = [c for c in keep if c in out.columns]
        return out[keep].reset_index(drop=True)

    def compounds(self, input_comp: pd.DataFrame, pChEMBL_thres: float=3.0) -> tuple[pd.DataFrame, str, pd.DataFrame]:

        """
        Retrieves proteins from ChEMBL interacting with compound passed as input.

        Steps
        -----
        - Finds chembl id of all synonyms for the compound passed as input \\
        Constraints:
            - Only Homo Sapiens interactions
            - Remove various activity comments (Not Determined, Not Evaluated, inconclusive, undetermined, No data)
            - Only interactions extracted from scientific literature
            - Only interactions with no data validity comments
            - Only interactions with pchembl value
        - Identifies target proteins using chembl target ids
        - Uses Biomart or mygene to obtain proteins info to return,
        searching with chembl id.
        
        Parameters
        ----------
        input_comp : DataFrame
            Dataframe of input compounds from which interacting proteins are found
        chembl_ids : list
            Empty list that will contain all ChEMBL ids found for the input compound
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid

        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values: \\
            entrez, gene_type, hgnc_symbol, description, datasource (chembl), pchembl_value
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all ChEMBL info about the input compound
        """
        
        columns = ['entrez','gene_type','hgnc_symbol','description','datasource', 'pchembl_value']
        # Create an empty DataFrame with the specified columns
        chembl_act = pd.DataFrame(columns=columns)
        chembl_raw = pd.DataFrame()

        input_comp_id = input_comp['inchikey_fb'][0]
        chembl_raw = self.data_manager.retrieve_raw_data('FirstBlock', input_comp_id)
    
        try:
            chembl_raw.reset_index()
        except:
            chembl_raw=pd.DataFrame()
            
        if len(chembl_raw) > 0:
            # Molecular weight if available converts mass/volume-based concentration units (mg/dl, ug/mL, ...) into a proper molar concentration.
            mol_weight = pd.to_numeric(input_comp['molecular_weight'], errors='coerce')[0] if 'molecular_weight' in input_comp.columns else None

            # Filter ChEMBL for high quality interactions. This follows: Bosc, N. et al. J. cheminformatics 11, 1–16 (2019)
            chembl_act = self._filter_database(chembl_raw, pChEMBL_thres, molecular_weight=mol_weight)

            if len(chembl_act) > 0:
                # Keep only relevant columns
                keep_cols = ['molecule_chembl_id', 'inchikey', 'CID', 'activity_comment', 'data_validity_comment', 'organism',
                             'src_id', 'pchembl_value', 'target_chembl_id', 'target_type', 'uniprot',
                             'standard_type', 'standard_relation', 'standard_value', 'standard_units',
                             'pChEMBL_lt', 'pChEMBL_gt']
                chembl_act = chembl_act[[c for c in keep_cols if c in chembl_act.columns]]
                # Only use int. w/proteins, removes cell lines and RNA
                chembl_act = chembl_act.loc[chembl_act['target_type']=='SINGLE PROTEIN'].reset_index(drop=True)
                    
                if len(chembl_act) > 0:
                    # Unify gene identifiers by Harmonizing IDs
                    ensembl = self.gene_server

                    # ChEMBL uses chembl ids
                    input_type='uniprotswissprot' 
                    attributes = ['uniprotswissprot', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                    names = ['uniprot','entrez','gene_type','hgnc_symbol','description']
                    chembl_targets=pd.DataFrame(columns=names)
                    
                    input_genes=list(chembl_act['uniprot'])
                        
                    chembl_targets=ensembl.subset_search(input_type, input_genes, attributes, names)

                    # For each compound, assign specific biomart column values to the ones from the original chembl database
                    for index, row in chembl_act.iterrows():
                        S1 = chembl_targets.loc[chembl_targets['uniprot']==row['uniprot']]
                        if len(S1) > 0:
                            chembl_act.loc[index, 'entrez'] = S1['entrez'].iloc[0]
                            chembl_act.loc[index, 'gene_type'] = S1['gene_type'].iloc[0]
                            chembl_act.loc[index, 'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                            chembl_act.loc[index, 'description'] = S1['description'].iloc[0]
                            chembl_act.loc[index, 'note'] ='Harmonized gene ID'
                        else:
                            chembl_act.loc[index, 'entrez'] = None
                            chembl_act.loc[index, 'gene_type'] = None
                            chembl_act.loc[index, 'hgnc_symbol'] = None
                            chembl_act.loc[index, 'description'] = None  
                            chembl_act.loc[index, 'note'] ='Failed to harmonize gene ID'
                    chembl_act['datasource'] = 'ChEMBL'

                    # Restructure to the shared cross-database output shape
                    compound_id_cols = ['molecule_chembl_id', 'inchikey', 'CID']
                    protein_id_cols = ['target_chembl_id', 'uniprot', 'entrez', 'hgnc_symbol', 'gene_type', 'description']
                    chembl_act = self._restructure_output(chembl_act, compound_id_cols, protein_id_cols)

                    statement = 'completed'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'Filter reduced interactions to 0'
        else:
            statement = 'No interaction data'

        return chembl_act, statement, chembl_raw
    
    def proteins(self, input_protein: pd.DataFrame, pChEMBL_thres: float=3.0) -> tuple[pd.DataFrame, str, pd.DataFrame]:

        """
        Retrieves compounds from ChEMBL database interacting with proteins passed as input.

        Steps
        -----
        - Retrieve compounds interacting with input protein \\
        Constraints:
            - Only Homo Sapiens interactions
            - Remove various activity comments (Not Determined, Not Evaluated, inconclusive, undetermined, No data)
            - Only interactions extracted from scientific literature
            - Only interactions with no data validity comments
            - Only interactions with pchembl value
        - Uses Pubchempy to retrieve the compound info (with additional data from chembl) to return, 
        searching with chEMBL ID or with inchi if the first search is unsuccessful.
        
        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compound are found
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid

        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values: \\
            inchi, inchikey, smiles, iupac_name, datasource (BindingDB), pchembl_value, notes (activity type)
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all BindingDB info about the protein
        """

        columns = ['inchi','inchikey','smiles','connectivity_smiles','iupac_name','molecular_weight','datasource','pchembl_value']
        # Create an empty DataFrame with the specified columns
        chembl_c1 = pd.DataFrame(columns=columns) 
        chembl_raw = pd.DataFrame()
        # Drop duplicated and null value of chembl_id, to make sure the first line of the dataset has the chembl ID
        input_protein = input_protein.dropna(subset=['chembl'])
        input_protein = input_protein.drop_duplicates(subset='chembl').reset_index(drop=True)        
        # Check if there are any input proteins remaining
        if len(input_protein) > 0:
            input_protein_id = input_protein['chembl'].iloc[0]

            chembl_raw = self.data_manager.retrieve_raw_data('target_chembl_id', input_protein_id)
            # Check if therer are any interacting compound
            if len(chembl_raw) > 0:
                # Filter ChEMBL for high quality interactions. Phase 1 only
                chembl_act = self._filter_database_pre(chembl_raw, pChEMBL_thres)

                if len(chembl_act) > 0:

                    pc = PubChemServer()

                    extra_cols = [c for c in ['standard_type', 'standard_relation', 'standard_value', 'standard_units',
                                               'pChEMBL_lt', 'pChEMBL_gt', '_pending_mw', 'activity_comment'] if c in chembl_act.columns]

                    if 'CID' in chembl_act.columns:
                        # Retrieve compounds using cids
                        compounds = self._pubchem_search_cid(chembl_act, columns, pc)
                        
                        if len(compounds) > 0:
                            chembl_info = chembl_act[['CID', 'molecule_chembl_id', 'pchembl_value', 'target_chembl_id', 'uniprot'] + extra_cols].\
                                            rename(columns={'molecule_chembl_id': 'id'}).drop_duplicates()
                            chembl_c1 = pd.merge(compounds, chembl_info, on='CID', how='left')
                        
                    else:
                        # Retrieve compounds using chembl ids and inchis if formers miss
                        compounds = self._pubchem_search_chembl(chembl_act, 'Molecule ChEMBL ID', columns, pc)

                        if len(compounds) > 0:
                            # Update pchembl value with that of filtered compounds
                            chembl_info = chembl_act[['molecule_chembl_id', 'pchembl_value', 'target_chembl_id', 'uniprot'] + extra_cols].\
                                            rename(columns={'molecule_chembl_id': 'id'}).drop_duplicates()
                            chembl_c1 = pd.merge(compounds, chembl_info, on='id', how='left')

                    # Check if at least one compound has been found in total
                    if len(chembl_c1) > 0:
                        # Phase 2: now that each candidate compound's molecular_weight is known, resolve any row still pending
                        mol_weight = pd.to_numeric(chembl_c1['molecular_weight'], errors='coerce') if 'molecular_weight' in chembl_c1.columns else None
                        chembl_c1 = self._filter_database_post(chembl_c1, pChEMBL_thres, molecular_weight=mol_weight)

                    if len(chembl_c1) > 0:
                        # Store into dataframe additional data
                        chembl_c1.loc[:, 'datasource'] = 'ChEMBL'
                        chembl_c1.loc[:, 'notes'] = np.nan

                        # Restructure to the shared cross-database output shape
                        compound_id_cols = ['inchikey', 'CID', 'id']
                        protein_id_cols = ['target_chembl_id', 'uniprot']
                        chembl_c1 = self._restructure_output(chembl_c1, compound_id_cols, protein_id_cols)

                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using Pubchem and ChEMBL'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein does not contain chembl id'

        return chembl_c1, statement, chembl_raw