'''Loading,searching,filtering and preprocessing data from OTP.'''

import pandas as pd
import numpy as np
import requests
import time
import json

from ..utils.typing import Connection
from ..servers.BiomartServer import BiomartServer
from ..servers.ChEMBLServer import ChEMBLServer as chembl
from ..servers.PubChemServer import PubChemServer
from ..data_manager import *
from .Database import Database

class OTP(Database):

    def __init__(self, connection:Connection|None=None, database:pd.DataFrame|None=None, chembl_database:pd.DataFrame|None=None, merge_stereoisomers=False):
        """
        Initialize OTP database.
    
        Parameters
        ----------
        database : pd.DataFrame, optional
            Preprocessed OTP DataFrame with columns: chemblIds, targets, targetName,
            inchikey, FirstBlock, CID. If provided, will use local data instead of API.
        chembl_database : pd.DataFrame, optional
            Local ChEMBL database for faster lookups
        merge_stereoisomers : bool
            Whether to merge stereoisomers
        """
        super().__init__(merge_stereoisomers)
    
        # OTP uses either local DataFrame or API (no SQL connection)
        if connection is not None:
            raise ValueError('OTP does not support SQL connections. Use database parameter for local data or leave empty for API.')
    
        self.chembl_db = chembl_database
    
        # Use local data if provided
        if database is not None:
            self._local_data = database
        else:
            # No local data - will use API
            self._local_data = None
            funcs = {
                'met_comps' : self._retrieve_met_compounds,
                'proteins' : self._retrieve_proteins 
            }
            self.data_manager = APIManager(funcs)

    def interactions(self, input_comp: pd.DataFrame, chembl_ids: list, merge_stereoisomers: bool=False) -> tuple[pd.DataFrame, str, pd.DataFrame]:
        """
        Retrieves proteins from OTP interacting with compound passed as input.

        Steps
        -----
        - Retrieves Mechanism interactions using the chembl ids passed as input.
        - Uses Biomart to obtain proteins info (modified with data from first API search) to return,
        searching with ensembl peptide id.
        
        Parameters
        ----------
        input_comp : DataFrame
            Dataframe of input compounds from which interacting proteins are found
        chembl_ids : list
            list of chembl ids for the input compounds
        merge_stereoisomers : bool
            Whether to merge stereoisomers
            
        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values: \\
            entrez, gene_type, hgnc_symbol, description, datasource (OTP), pchembl_value (nan)
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all OTP info about the input compound
        """
                
        columns = ['entrez','gene_type','hgnc_symbol','description','pchembl_value','datasource']
        otp_act = pd.DataFrame(columns=columns)
        raw_columns = ['Chembl','mID','mLabel']
        otp_raw = pd.DataFrame(columns=raw_columns)
        
        # Use local data only
        if self._local_data is not None:
            # Get ChEMBL IDs if not provided
            if not chembl_ids:
                input_comp = input_comp.dropna(subset=['inchikey']).reset_index(drop=True)
                
                if len(input_comp) > 0:
                    if merge_stereoisomers:
                        # Use FirstBlock for stereoisomer matching
                        inchikey_fb = input_comp['inchikey_fb'][0]
                        if 'FirstBlock' in self._local_data.columns:
                            otp_matches = self._local_data[self._local_data['FirstBlock'] == inchikey_fb]
                            if len(otp_matches) > 0:
                                chembl_ids.extend(otp_matches['chemblIds'].dropna().unique().tolist())
                    else:
                        # Use exact inchikey matching
                        inchikey = input_comp['inchikey'][0]
                        if 'inchikey' in self._local_data.columns:
                            otp_matches = self._local_data[self._local_data['inchikey'] == inchikey]
                            if len(otp_matches) > 0:
                                chembl_ids.extend(otp_matches['chemblIds'].dropna().unique().tolist())
            
            # Query local data
            if len(chembl_ids) > 0 and None not in chembl_ids:
                otp_raw = self._query_local_moa(chembl_ids)
                
                if len(otp_raw) > 0:
                    otp_filtered = otp_raw
                    if len(otp_filtered) > 0:
                        otp_act = pd.concat([otp_act, otp_filtered])

                        # Get mapping directly from local OTP data
                        unique_chembl_ids = otp_act['Chembl'].unique().tolist()
                        otp_mapping = self._local_data[self._local_data['chemblIds'].isin(unique_chembl_ids)]
                        
                        if len(otp_mapping) > 0:
                            # Select relevant columns and drop duplicates
                            mapping_cols = ['chemblIds', 'inchikey']
                            if 'CID' in self._local_data.columns:
                                mapping_cols.append('CID')
                            
                            otp_map = otp_mapping[mapping_cols].drop_duplicates()
                            
                            # Merge into otp_act
                            otp_act = otp_act.merge(otp_map, left_on='Chembl', right_on='chemblIds', how='left')
                            otp_act = otp_act.drop(columns=['chemblIds'])
                        else:
                            otp_act['inchikey'] = None
                            if 'CID' not in otp_act.columns:
                                otp_act['CID'] = None

                        # Convert OTP targets using Biomart
                        ensembl = BiomartServer()
                        attributes = ['ensembl_gene_id', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                        names = ['ensembl_id','entrez','gene_type','hgnc_symbol','description']
                        input_type = 'ensembl_gene_id'
                        
                        otp_targets = pd.DataFrame(columns=names)
                        input_genes = list(otp_act['mID'])
                        otp_targets = ensembl.subset_search(input_type, input_genes, attributes, names)
                        
                        for index, row in otp_act.iterrows():
                            S1 = otp_targets.loc[otp_targets['ensembl_id']==row['mID']]
                            if len(S1) > 0:
                                otp_act.loc[index,'entrez'] = S1['entrez'].iloc[0]
                                otp_act.loc[index,'gene_type'] = S1['gene_type'].iloc[0]
                                otp_act.loc[index,'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                                otp_act.loc[index,'description'] = S1['description'].iloc[0]
                            else:
                                otp_act.loc[index,'entrez'] = None
                                otp_act.loc[index,'gene_type'] = None
                                otp_act.loc[index,'hgnc_symbol'] = None
                                otp_act.loc[index,'description'] = None
                        
                        otp_act['datasource'] = 'OTP'
                        otp_act['pchembl_value'] = np.nan
                        statement = 'completed'
                    else:
                        statement = 'Filter reduced interactions to 0'
                else:
                    statement = 'No interaction data in local OTP'
            else:
                statement = 'No ChEMBL ids found in local OTP'
        
        # Use API only
        else:
            # Get ChEMBL IDs if not provided
            if not chembl_ids:
                # Try ChEMBL database first if available
                if self.chembl_db is not None:
                    input_comp = input_comp.dropna(subset=['inchikey']).reset_index(drop=True)

                    if len(input_comp) > 0:
                        if merge_stereoisomers:
                            inchikey_fb = input_comp['inchikey_fb'][0]
                            chembl_matches = self.chembl_db[self.chembl_db['FirstBlock'] == inchikey_fb]
                            if len(chembl_matches) > 0:
                                chembl_ids.extend(chembl_matches['Molecule ChEMBL ID'].dropna().unique().tolist())
                        else:
                            inchikey = input_comp['inchikey'][0]
                            inchikey_col = 'standard_inchi_key' if 'standard_inchi_key' in self.chembl_db.columns else 'inchikey'
                            chembl_matches = self.chembl_db[self.chembl_db[inchikey_col] == inchikey]
                            if len(chembl_matches) > 0:
                                chembl_ids.extend(chembl_matches['Molecule ChEMBL ID'].dropna().unique().tolist())
            
                # Use ChEMBL API for chembl ids
                if not chembl_ids:
                    print("Using ChEMBL API to find ChEMBL IDs (slower)")
                    chembl_ids.extend(chembl.identify_chembl_ids(input_comp))

            if len(chembl_ids) > 0 and None not in chembl_ids:
                otp_raw = self.data_manager.retrieve_raw_data('met_comps', chembl_ids)

                if len(otp_raw) > 0:
                    otp_filtered = otp_raw
                    if len(otp_filtered) > 0:
                        otp_act = pd.concat([otp_act, otp_filtered])

                        # Add InChIKey and CID from ChEMBL database if available
                        if self.chembl_db is not None:
                            unique_chembl_ids = otp_act['Chembl'].unique().tolist()
                            chembl_inchikey_map = self.chembl_db[self.chembl_db['Molecule ChEMBL ID'].isin(unique_chembl_ids)]
                        
                            if len(chembl_inchikey_map) > 0:
                                inchikey_col = 'standard_inchi_key' if 'standard_inchi_key' in self.chembl_db.columns else 'inchikey'
                                chembl_map = chembl_inchikey_map[['Molecule ChEMBL ID', inchikey_col, 'CID']].drop_duplicates()
                                chembl_map = chembl_map.rename(columns={inchikey_col: 'inchikey'})
                            
                                otp_act = otp_act.merge(chembl_map, left_on='Chembl', right_on='Molecule ChEMBL ID', how='left')
                                otp_act = otp_act.drop(columns=['Molecule ChEMBL ID'])
                            else:
                                otp_act['inchikey'] = None
                                otp_act['CID'] = None
                        else:
                            otp_act['inchikey'] = None
                            otp_act['CID'] = None

                        # Convert OTP targets using Biomart
                        ensembl = BiomartServer()
                        attributes = ['ensembl_gene_id', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                        names = ['ensembl_id','entrez','gene_type','hgnc_symbol','description']
                        input_type = 'ensembl_gene_id'
                        
                        otp_targets = pd.DataFrame(columns=names)
                        input_genes = list(otp_act['mID'])
                        otp_targets = ensembl.subset_search(input_type, input_genes, attributes, names)
                        
                        for index, row in otp_act.iterrows():
                            S1 = otp_targets.loc[otp_targets['ensembl_id']==row['mID']]
                            if len(S1) > 0:
                                otp_act.loc[index,'entrez'] = S1['entrez'].iloc[0]
                                otp_act.loc[index,'gene_type'] = S1['gene_type'].iloc[0]
                                otp_act.loc[index,'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                                otp_act.loc[index,'description'] = S1['description'].iloc[0]
                            else:
                                otp_act.loc[index,'entrez'] = None
                                otp_act.loc[index,'gene_type'] = None
                                otp_act.loc[index,'hgnc_symbol'] = None
                                otp_act.loc[index,'description'] = None
                        
                        otp_act['datasource'] = 'OTP'
                        otp_act['pchembl_value'] = np.nan
                        statement = 'completed'
                    else:
                        statement = 'Filter reduced interactions to 0'
                else:
                    statement = 'No interaction data from API'
            else:
                statement = 'No ChEMBL ids found'
            
        return otp_act, statement, otp_raw

    def compounds(self, input_protein: pd.DataFrame, merge_stereoisomers: bool=False) -> tuple[pd.DataFrame, str, pd.DataFrame]:
        """
        Retrieves compounds from OTP database interacting with proteins passed as input.
        Uses local data first then use API if local is not available.
    
        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compounds are found
        merge_stereoisomers : bool
            Not used in this method

        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values: \\
            inchi, inchikey, isomeric_smiles, iupac_name, datasource (OTP), pchembl_value, notes (nan)
        String
            Statement describing outcome
        DataFrame
            Raw data
        """
        
        columns = ['inchi','inchikey','isomeric_smiles','iupac_name','datasource','pchembl_value']
        # Create an empty DataFrame with the specified columns
        otp_c1 = pd.DataFrame(columns=columns)
        otp_raw = pd.DataFrame()

        # Drop duplicated and null values of ensembl gene id
        input_protein = input_protein.dropna(subset=['ensembl_gene_id'])
        input_protein = input_protein.drop_duplicates(subset='ensembl_gene_id').reset_index(drop=True)    

        # Check if there are any input proteins remaining  
        if len(input_protein) > 0:
            input_protein_id = input_protein['ensembl_gene_id'][0]
        
            # Use local data only
            if self._local_data is not None:
                otp_raw = self._query_local_moa_by_target(input_protein_id)
            
                if len(otp_raw) > 0:
                    pc = PubChemServer()
                
                    # If local data has CID
                    if 'CID' in otp_raw.columns:
                        # Retrieve compounds using CIDs
                        compounds = self._pubchem_search_cid(otp_raw, columns, pc)
                    
                        if len(compounds) > 0:
                            otp_c1 = compounds
                
                    # Use InChIKey if no CID
                    if len(otp_c1) == 0 and 'inchikey' in otp_raw.columns:
                        # Filter only columns from pubchempy
                        selected_columns = pc.get_columns(columns[:-2])
                        compounds = pd.DataFrame()
                        ids = list(otp_raw['inchikey'].unique())
                    
                        for id in ids:
                            try:
                                comp = pc.get_compounds(id, selected_columns, namespace='inchikey')
                                compounds = pd.concat([compounds, comp])
                            except:
                                None
                            time.sleep(0.5)
                    
                        # Check if at least a match has been found
                        if len(compounds) > 0:
                            otp_c1 = compounds.rename(columns=pc.properties)
                    
                    # If still no results, try chemblIds
                    if len(otp_c1) == 0:
                        id_col = 'chemblIds' if 'chemblIds' in otp_raw.columns else None
                        if id_col:
                            otp_c1 = self._pubchem_search_chembl(otp_raw, id_col, columns, pc)

                    if len(otp_c1) > 0:
                        otp_c1['datasource'] = 'OTP'
                        otp_c1['pchembl_value'] = np.nan
                        otp_c1['notes'] = np.nan
                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using PubChem'
                else:
                    statement = 'No interaction data in local OTP'
            
            # Use API only
            else:
                otp_raw = self.data_manager.retrieve_raw_data('proteins', input_protein_id)
            
                if len(otp_raw) > 0:
                    pc = PubChemServer()
                    # Perform search using drugId and inchis
                    otp_c1 = self._pubchem_search_chembl(otp_raw, 'drugId', columns, pc)

                    if len(otp_c1) > 0:
                        otp_c1['datasource'] = 'OTP'
                        otp_c1['pchembl_value'] = np.nan
                        otp_c1['notes'] = np.nan
                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using PubChem'
                else:
                    statement = 'No interaction data from API'
        else:
            statement = 'Input protein does not contain ensembl id'

        return otp_c1, statement, otp_raw
    
    def _retrieve_met_compounds(self, chembl_ids):
        """Retrieve mechanism of action data from OTP API"""
        base_url = "https://api.platform.opentargets.org/api/v4/graphql"

        otp_met = pd.DataFrame(columns=['Chembl','mID','mLabel'])
        met_list = []

        for ChEMBLID in chembl_ids:
            # Set variables object of arguments for GraphQL API
            variables={"ChEMBL ID": ChEMBLID}
            # Set the Query for the GraphQL API to get the OTP bibliography
            # This inputs the ChEMBLID variable into the query_string at the %s spot
            query_string="""
            query {
                drug(chemblId: "%s"){
                    id
                    mechanismsOfAction{
                        rows{
                            mechanismOfAction
                            targetName
                                targets{
                                    id
                                    approvedSymbol
                                }
                                references {
                                    source
                                    urls
                                }
                        }
                    }
                }
            }
            """ % ChEMBLID 

            while True:
                try: # Perform request
                    dat=requests.post(base_url,json={"query": query_string,"variables": variables})
                    #print(dat.status_code) #200 code = Operating Normally
                    # Transform API response from JSON to dictionary
                    res=json.loads(dat.text) 
                    break
                except json.JSONDecodeError as e:
                    print(f"Error: {e}. Retrying in 10 seconds...")
                    time.sleep(10)

            # Check if a matching interaction has been found
            if res['data']['drug'] is not None:
                for mechanism_row in res['data']['drug']['mechanismsOfAction']['rows']:
                    # Check if the target name exists and if there's at least one target
                    if mechanism_row['targetName'] is not None and len(mechanism_row['targets']) > 0:
                        try: # Add protein to DataFrame
                            met_list.append({
                                'Chembl': ChEMBLID,
                                'mLabel': mechanism_row['targetName'],
                                'mID': mechanism_row['targets'][0]['id'],
                                'PMID': mechanism_row['references'][0]['urls'][0],
                                'mType': 'GP',
                                'strength': 'strong'
                            })
                        except:
                            continue

        # Create DataFrame
        otp_met = pd.concat([otp_met, pd.DataFrame(met_list)])
        otp_met = otp_met.drop_duplicates().reset_index(drop=True)

        return otp_met
    
    def _retrieve_proteins(self, input_protein_id) -> pd.DataFrame:
        # Create post request
        """Retrieve protein data from OTP API"""
        variables={"ensembl_gene_id": input_protein_id}
        base_url="https://api.platform.opentargets.org/api/v4/graphql"

        # Size 10000 is the maximum number of compounds that can be shown. 
        query = """
        query search{
            target(ensemblId:"%s"){
                approvedSymbol
                id
                knownDrugs(size:10000){ 
                    rows{
                        drugId
                        prefName
                        drugType
                    }
                }
            }
        }
        """% input_protein_id
        dat = requests.post(base_url,json={"query": query,"variables": variables})
        # Store data from API
        res = json.loads(dat.text)
        try: # Save relevant API data into dataframe
            otp_raw = pd.DataFrame(res['data']['target']['knownDrugs']['rows'])
            # Remove duplicates of drugId
            otp_raw = otp_raw.drop_duplicates(subset='drugId').reset_index(drop=True)
        except TypeError:
            # There was an error with the data types, return empty dataframe
            otp_raw = pd.DataFrame(columns=['inchi','inchikey','isomeric_smiles','iupac_name','datasource','pchembl_value'])

        return otp_raw

    def _expand_chembl_ids_to_stereoisomers(self, chembl_ids):
        """Expands ChEMBL IDs to include all stereoisomers by looking up first-blocks"""
        if self.chembl_db is None or not self.merge_stereoisomers:
            return chembl_ids
        
        chembl_lookup = self.chembl_db[self.chembl_db['Molecule ChEMBL ID'].isin(chembl_ids)]
        
        if len(chembl_lookup) == 0:
            return chembl_ids
        
        inchikeys = chembl_lookup['inchikey'].dropna().unique()
        pc = PubChemServer()
        first_blocks = [pc.get_inchikey_first_block(ik) for ik in inchikeys]
        first_blocks = list(set(first_blocks))

        expanded_chembl_ids = []
        for fb in first_blocks:
            matches = self.chembl_db[self.chembl_db['FirstBlock'] == fb]
            expanded_chembl_ids.extend(matches['Molecule ChEMBL ID'].tolist())
        
        return list(set(expanded_chembl_ids))
    
    def _query_local_moa(self, chembl_ids):
        """Query local mechanism of action data by ChEMBL IDs"""
        local_matches = self._local_data[self._local_data['chemblIds'].isin(chembl_ids)].copy()
        
        if len(local_matches) == 0:
            return pd.DataFrame(columns=['Chembl','mID','mLabel'])
        
        # Convert to format expected by the rest of the code
        otp_raw = pd.DataFrame()
        otp_raw['Chembl'] = local_matches['chemblIds']
        otp_raw['mID'] = local_matches['targets']
        otp_raw['mLabel'] = local_matches['targetName']
        
        return otp_raw

    def _query_local_moa_by_target(self, ensembl_gene_id):
        """Query local MOA data by target/protein Ensembl ID"""
        local_matches = self._local_data[self._local_data['targets'] == ensembl_gene_id].copy()
        
        if len(local_matches) == 0:
            return pd.DataFrame()
        
        # Build result with available columns
        result_cols = ['chemblIds']
        
        # Include inchikey and CID if available in local data
        if 'inchikey' in local_matches.columns:
            result_cols.append('inchikey')
        if 'CID' in local_matches.columns:
            result_cols.append('CID')
        
        otp_raw = local_matches[result_cols].copy()
        
        # Remove duplicates by chemblIds
        otp_raw = otp_raw.drop_duplicates(subset='chemblIds').reset_index(drop=True)
        
        return otp_raw
