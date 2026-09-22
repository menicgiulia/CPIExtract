'''Loading,searching,filtering and preprocessing data from BindingDB.'''

import numpy as np
import pandas as pd
import re
import time

from ..utils.typing import Connection
from ..servers.PubChemServer import PubChemServer
from ..servers.BiomartServer import BiomartServer
from ..servers.MyGeneServer import MyGeneServer
from .Database import Database
from ..data_manager import *

class BindingDB(Database):
    '''Loading,searching,filtering and preprocessing data from BindingDB.'''

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
        
       # if not connection and not database:
        #     raise ValueError('Either SQL connection or database should be not None')
        if database is not None:
            self.data_manager = LocalManager(database)
        else:
            self.data_manager = SQLManager(connection, 'BDB')

    # 'Target Name' follows a consistent, parseable bracket convention: 'Protein Name' (wild-type, no brackets),
    # 'Protein Name [123-456]' (a truncated/domain construct, no mutation), or
    # 'Protein Name [123-456,Q491K,I497V,...]' (explicit point mutations after the residue range).
    _MUTATION_TOKEN_RE = re.compile(r'[A-Z]\d+[A-Z]')

    @classmethod
    def _has_mutation(cls, target_name) -> bool:
        if pd.isnull(target_name):
            return False
        for bracket_content in re.findall(r'\[([^\]]*)\]', str(target_name)):
            if cls._MUTATION_TOKEN_RE.search(bracket_content):
                return True
        return False

    def _filter_database(self, bdb_raw: pd.DataFrame) -> pd.DataFrame:
        """
        Filters the bdb database with the following constraints:
            - Only Homo Sapiens interactions
            - Temp(°C) < 40 or null
            - 5 < pH < 9 or null
            - Single-chain targets only (excludes multi-chain complexes)
            - Excludes targets whose name indicates a point mutation
        
        Parameters
        ----------
        bdb_raw : DataFrame
            bdb database

        Returns
        -------
        DataFrame
            Filtered database
        """
        # Remove all non-human interactions
        bdb_filt = bdb_raw.loc[bdb_raw['Target Source Organism According to Curator or DataSource']=='Homo sapiens'].copy()
        # Remove >, <, and C characters
        bdb_filt['Temp (C)'] = bdb_filt['Temp (C)'].str.replace(r'[^0-9\.]','',regex=True)
        # Convert strings to numeric values 
        bdb_filt['Temp (C)'] = bdb_filt['Temp (C)'].apply(pd.to_numeric, errors='coerce') 
        # Filter temperature (also accept null values)
        bdb_filt = bdb_filt.loc[(bdb_filt['Temp (C)'].isnull()) | (bdb_filt['Temp (C)'] < 40) &
                            # Filter pH (also accept null values)
                            ((bdb_filt['pH'].isnull()) | (bdb_filt['pH'] < 9) & (bdb_filt['pH'] > 5))]\
                                .drop_duplicates(ignore_index=True)

        # Exclude multi-chain protein complexes
        n_chains = pd.to_numeric(bdb_filt['Number of Protein Chains in Target (>1 implies a multichain complex)'], errors='coerce')
        bdb_filt = bdb_filt.loc[n_chains.isnull() | (n_chains <= 1)]

        # Exclude mutant/variant targets, based on the target name's bracket notation
        is_mutant = bdb_filt['Target Name'].apply(self._has_mutation)
        bdb_filt = bdb_filt.loc[~is_mutant]

        return bdb_filt.reset_index(drop=True)

    def _compute_pchembl(self, bdb_dat: pd.DataFrame, pChEMBL_thres: float) -> pd.DataFrame:
        """
        Computes a pChEMBL-equivalent value for each non-null activity column
        - Ki (nM)
        - IC50 (nM)
        - Kd (nM)
        - EC50 (nM)

        Unlike ChEMBL, BindingDB reports each of these as a single string that may carry a
        leading relation operator directly embedded in it (e.g. '>10000', '<0.5') rather than
        a separate relation column. The operator is extracted BEFORE the value is cleaned to a
        plain number, and handled exactly as in ChEMBL's _filter_database_pre: since these are
        always concentrations (BindingDB doesn't report association constants or p-scale values
        in these columns), relation '>' -> weaker/upper bound (true pChEMBL < value), relation
        '<' -> stronger/lower bound (true pChEMBL > value). A '>'-censored row only counts as a
        confirmed negative if the bound itself is <= pChEMBL_thres; a '<'-censored row only
        counts as a confirmed positive if the bound clears pChEMBL_thres. Anything else (an
        uninformative bound, e.g. 'Ki > 1 nM') is dropped rather than misclassified. An exact
        ('=') value is kept regardless of which side of the threshold it lands on - the
        threshold classifies known values, it doesn't gate them out.

        Parameters
        ----------
        bdb_dat : DataFrame
            bdb database
        pChEMBL_thres : float
            pChEMBL value used to classify a resolved interaction as positive (above) or
            negative (at or below)
        
        Returns
        -------
        DataFrame
            One row per (compound-protein record, measurement type) that passed filtering,
            with pChEMBL_eq / pChEMBL_lt / pChEMBL_gt and standard_type (Ki/IC50/Kd/EC50) columns.
        """
        cols = ['Ki (nM)','IC50 (nM)','Kd (nM)','EC50 (nM)']
        results = []

        for col in cols:
            sub = bdb_dat.dropna(subset=[col]).copy()
            if len(sub) == 0:
                continue

            raw_val = sub[col].astype(str)
            # Extract the relation operator before the value gets stripped down to digits.
            # No prefix (a bare number) defaults to '='.
            relation = raw_val.str.extract(r'^\s*(>=|<=|>|<)?')[0].fillna('=')
            numeric_str = raw_val.str.replace(r'[^0-9.]', '', regex=True)
            value_nM = pd.to_numeric(numeric_str, errors='coerce')

            valid = value_nM.notnull() & (value_nM != 0)
            sub = sub.loc[valid].copy()
            if len(sub) == 0:
                continue
            relation = relation.loc[valid]
            value_nM = value_nM.loc[valid]

            derived = -np.log10(value_nM * 1e-9)
            is_gt = relation.isin(['>', '>='])
            is_lt = relation.isin(['<', '<='])
            is_eq = relation == '='

            sub['pChEMBL_eq'] = derived.where(is_eq)
            sub['pChEMBL_lt'] = derived.where(is_gt)
            sub['pChEMBL_gt'] = derived.where(is_lt)
            sub['standard_type'] = col.split(' ')[0]

            results.append(sub.drop(columns=cols, errors='ignore'))

        if len(results) == 0:
            return bdb_dat.iloc[0:0].drop(columns=cols, errors='ignore').assign(
                pChEMBL_eq=pd.Series(dtype=float), pChEMBL_lt=pd.Series(dtype=float),
                pChEMBL_gt=pd.Series(dtype=float), standard_type=pd.Series(dtype=object))

        bdb_act = pd.concat(results, ignore_index=True)

        # A '>'-censored row is only genuine negative evidence if the bound itself rules out an
        # active result; a loose bound (e.g. 'Ki > 1 nM') doesn't confirm inactivity at all.
        true_negatives = bdb_act['pChEMBL_lt'].notnull() & (bdb_act['pChEMBL_lt'] <= pChEMBL_thres)
        # A '<'-censored row only counts as positive evidence if its lower bound clears the threshold.
        censored_positive = bdb_act['pChEMBL_gt'].notnull() & (bdb_act['pChEMBL_gt'] > pChEMBL_thres)
        # An exact value is kept regardless of which side of threshold it falls on.
        exact_value = bdb_act['pChEMBL_eq'].notnull()

        bdb_act = bdb_act.loc[exact_value | true_negatives | censored_positive].reset_index(drop=True)

        return bdb_act

    def _restructure_output(self, bdb_act: pd.DataFrame, compound_id_cols: list, protein_id_cols: list) -> pd.DataFrame:
        """
        Restructures filtered, row-per-activity-record BindingDB data into the shared
        row-per-passed-interaction shape used across CPIExtract databases:

        compound_identifier cols | protein_identifier cols | standard_type | standard_units
        | pChEMBL_eq | pChEMBL_lt | pChEMBL_gt | datasource

        Mirrors ChEMBL's _restructure_output. standard_units is always 'nM' here (BindingDB's
        four activity columns are all reported in nM), kept as an explicit column purely for
        schema consistency with ChEMBL's output ahead of the downstream cross-database merge.
        No aggregation happens here - the same compound-protein pair can appear on multiple
        rows (once per activity record/measurement type that passed filtering).

        Parameters
        ----------
        bdb_act : DataFrame
            Filtered, row-per-activity-record BindingDB data.
        compound_id_cols : list
            Columns that uniquely identify a compound (e.g. ['inchikey', 'CID']).
        protein_id_cols : list
            Columns that uniquely identify a protein (e.g. ['uniprot']).

        Returns
        -------
        DataFrame
            One row per activity record that passed filtering.
        """
        compound_id_cols = [c for c in compound_id_cols if c in bdb_act.columns]
        protein_id_cols = [c for c in protein_id_cols if c in bdb_act.columns]

        if len(bdb_act) == 0:
            return pd.DataFrame()

        out = bdb_act.copy()
        out['datasource'] = 'BindingDB'
        out['standard_units'] = 'nM'

        keep = compound_id_cols + protein_id_cols + ['standard_type', 'standard_units', 'pChEMBL_eq', 'pChEMBL_lt', 'pChEMBL_gt', 'datasource']
        keep = [c for c in keep if c in out.columns]
        return out[keep].reset_index(drop=True)


    def compounds(self, input_comp: pd.DataFrame, pChEMBL_thres: float=3.0) -> tuple[pd.DataFrame, str, pd.DataFrame]:
        """
        Retrieves proteins from bdb database interacting with compound passed as input.

        Steps
        -----
        - Filters bdb database to obtain proteins interacting with input compound \\
        Constraints:
            - Only Homo Sapiens interactions
            - Temp(°C) < 40 or null
            - 5 < pH < 9 or null
        - Uses activity columns to compute all possible pchembl values.
        
        Parameters
        ----------
        input_comp : dictionary
            dictionary of input compound data from which interacting proteins are found
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid
            
        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values: \\
            entrez, gene_type, hgnc_symbol, description, datasource (BindingDB), pChEMBL_eq
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all BindingDB info about the input compound
        """

        columns = ['inchikey','entrez','gene_type','hgnc_symbol','description','pChEMBL_eq','datasource']
        # Create an empty DataFrame with the specified columns
        bdb_act = pd.DataFrame(columns=columns)
        bdb_raw = pd.DataFrame(columns=columns)
        # Drop null values of inchikey, to make sure the first line of the dataframe has the inchi key 
        input_comp = input_comp.dropna(subset=['inchikey']).reset_index(drop=True)
        # Check if there are any input compounds remaining  
        if len(input_comp) > 0:
            input_comp_id = input_comp['inchikey_fb'][0]
            bdb_raw = self.data_manager.retrieve_raw_data('FirstBlock', input_comp_id)
            
            if len(bdb_raw) > 0:
                # Filter database
                bdb_filt = self._filter_database(bdb_raw)
                bdb_filt.drop(['Ligand SMILES','Ligand InChI'], axis=1, inplace=True)
                # Compute pchembl values
                bdb_act = self._compute_pchembl(bdb_filt, pChEMBL_thres)
                
                if len(bdb_act) > 0:
                    # Unify gene identifiers by Harmonizing IDs
                    ensembl = self.gene_server

                    # BindingDB uses Uniprot IDs
                    input_type='uniprotswissprot' 
                    attributes = ['uniprotswissprot', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                    names = ['uniprot','entrez','gene_type','hgnc_symbol','description']
                    bdb_targets = pd.DataFrame(columns=names)
                    
                    input_genes = list(bdb_act['UniProt (SwissProt) Primary ID of Target Chain 1'])
                    
                    bdb_targets = ensembl.subset_search(input_type, input_genes, attributes, names)
    
                    # For each compound, assign protein query column values to the ones from the original database
                    for index, row in bdb_act.iterrows():
                        # Find the compound using uniprot
                        S1 = bdb_targets.loc[bdb_targets['uniprot'] == row['UniProt (SwissProt) Primary ID of Target Chain 1']]
                        if len(S1) > 0:
                            bdb_act.loc[index, 'entrez'] = S1['entrez'].iloc[0]
                            bdb_act.loc[index, 'gene_type'] = S1['gene_type'].iloc[0]
                            bdb_act.loc[index, 'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                            bdb_act.loc[index, 'description'] = S1['description'].iloc[0]
                            bdb_act.loc[index, 'note'] ='Harmonized gene ID'
                        else:
                            bdb_act.loc[index, 'entrez'] = None
                            bdb_act.loc[index, 'gene_type'] = None
                            bdb_act.loc[index, 'hgnc_symbol'] = None
                            bdb_act.loc[index, 'description'] = None
                            bdb_act.loc[index, 'note'] ='Failed to harmonize gene ID'
                    bdb_act['datasource'] = 'BindingDB'
                    bdb_act['uniprot'] = bdb_act['UniProt (SwissProt) Primary ID of Target Chain 1']

                    # Restructure to the shared cross-database output shape
                    compound_id_cols = ['inchikey', 'CID']
                    protein_id_cols = ['uniprot', 'entrez', 'hgnc_symbol', 'gene_type', 'description']
                    bdb_act = self._restructure_output(bdb_act, compound_id_cols, protein_id_cols)

                    statement = 'completed'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input compound does not contain inchi key'

        return bdb_act, statement, bdb_raw
    

    def proteins(self, input_protein: pd.DataFrame, pChEMBL_thres: float=3.0) -> tuple[pd.DataFrame, str, pd.DataFrame]:
        """
        Retrieves compounds from bdb database interacting with proteins passed as input.

        Steps
        -----
        - Filters bdb database to obtain compounds interacting with input proteins:
        Constraints:
            - Only Homo Sapiens interactions
            - Temp(°C) < 40 or null
            - 5 < pH < 9 or null
        - Uses activity columns to compute all possible pchembl values. \\
        - Uses Pubchempy to obtain compounds info (modified with data from original bdb database) to return,
        searching with BindingDB ID.
        
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
            inchi, inchikey, smiles, iupac_name, datasource (BindingDB), pChEMBL_eq, standard_type (Ki/IC50/Kd/EC50)
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all BindingDB info about the protein
        """

        columns = ['inchi','inchikey','smiles','connectivity_smiles','iupac_name','datasource','pChEMBL_eq']
        # Create an empty DataFrame with the specified columns
        bdb_c1 = pd.DataFrame(columns=columns)
        bdb_raw = pd.DataFrame()
        # Drop duplicated and null values of uniprot, to make sure the first line of the dataframe has the uniprot ID
        input_protein= input_protein.dropna(subset=['uniprot'])
        input_protein= input_protein.drop_duplicates(subset='uniprot').reset_index(drop=True) 
        # Check if there are any input proteins remaining       
        if len(input_protein) > 0:
            input_protein_id=input_protein['uniprot'].iloc[0]
            # Select only compounds interacting with the input protein
            bdb_raw = self.data_manager.retrieve_raw_data('UniProt (SwissProt) Primary ID of Target Chain 1', input_protein_id)
            # bdb_raw = BDB_data.loc[(BDB_data['UniProt (SwissProt) Primary ID of Target Chain 1']==input_protein_id)]
            if len(bdb_raw):
                # Filter database
                bdb_c = self._filter_database(bdb_raw)
                
                # Compute pchembl values
                bdb_act = self._compute_pchembl(bdb_c, pChEMBL_thres)
                # Check if at least one interacting compound has been found
                if len(bdb_act) > 0:
                    pc = PubChemServer()

                    extra_cols = [c for c in ['standard_type', 'pChEMBL_lt', 'pChEMBL_gt'] if c in bdb_act.columns]

                    if 'CID' in bdb_act:
                        # Retrieve compounds using cids
                        compounds = self._pubchem_search_cid(bdb_act, columns, pc)

                        if len(compounds) > 0:    
                            # Add additional values from activity dataframe
                            bdb_info = bdb_act[['CID', 'UniProt (SwissProt) Primary ID of Target Chain 1', 
                                                'BindingDB Ligand Name', 'pChEMBL_eq'] + extra_cols].\
                                        rename(columns={'UniProt (SwissProt) Primary ID of Target Chain 1': 'uniprot'}).\
                                                        drop_duplicates()
                            
                            bdb_c1 = pd.merge(compounds, bdb_info, on='CID', how='left')
                    else:
                        # Filter only columns from pubchempy to return
                        selected_columns = pc.get_columns(columns[:-2])
                        compounds = pd.DataFrame()
                        # Compute BindingDB ID and use it to search the compound from Pubchempy
                        bdb_act.loc[:, 'BindingDB MonomerID'] = bdb_act['BindingDB MonomerID'].apply(lambda x: 'BDBM' + str(x))
                        for _, row in bdb_act.iterrows():        
                            bdb_id = row['BindingDB MonomerID'] 
                            try:
                                comp = pc.get_compounds(bdb_id, selected_columns, namespace='name')
                                compounds = pd.concat([compounds, comp])
                            except:
                                None
                            time.sleep(0.5)
                        # Check if at least a match has been found
                        if len(compounds) > 0:
                            bdb_c1 = compounds.rename(columns=pc.properties)

                            # Add additional values from activity dataframe
                            bdb_info = bdb_act[['inchikey', 'UniProt (SwissProt) Primary ID of Target Chain 1', 
                                                'BindingDB Ligand Name', 'pChEMBL_eq'] + extra_cols].\
                                        rename(columns={'UniProt (SwissProt) Primary ID of Target Chain 1': 'uniprot'}).\
                                                        drop_duplicates()
                            
                            bdb_c1 = pd.merge(bdb_c1, bdb_info, on='inchikey', how='left')

                    if len(bdb_c1) > 0:
                        # Restructure to the shared cross-database output shape
                        compound_id_cols = ['inchikey', 'CID']
                        protein_id_cols = ['uniprot']
                        bdb_c1 = self._restructure_output(bdb_c1, compound_id_cols, protein_id_cols)

                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using PubChem'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein does not contain uniprot id'
            
        return bdb_c1, statement, bdb_raw