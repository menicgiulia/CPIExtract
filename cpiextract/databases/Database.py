'''Template of a Database class.Loading,searching,filtering and preprocessing data from a specified database.'''

from abc import ABC, abstractmethod
from ..servers.ChEMBLServer import ChEMBLServer as chembl
from ..servers.PubChemServer import PubChemServer
from ..utils.helper import *
import pandas as pd
import time
import pubchempy as pcp

class Database(ABC):
'''Template of a Database class.Loading,searching,filtering and preprocessing data from a specified database.'''
    def __init__(self, merge_stereoisomers=False):
        """
        Initialize database with stereoisomer matching preference.
        
        Parameters
        ----------
        merge_stereoisomers : bool, default False
            Whether to match compounds using only the first block of inchikey.
            - False: Exact inchikey match
            - True: First-block match (all stereoisomers)
        """
        self.merge_stereoisomers = merge_stereoisomers

    @abstractmethod
    def interactions(self, *args, **kwargs):
        """
        Retrieves proteins from the database interacting with compound passed as input.
        """
        raise NotImplementedError

    @abstractmethod
    def compounds(self, *args, **kwargs):
        """
        Retrieves compounds from the database interacting with proteins passed as input.
        """
        raise NotImplementedError
    
    def _pubchem_search_cid(self, db_act: pd.DataFrame, columns: list[str], pc: PubChemServer) -> pd.DataFrame:
	'''Search componds from pubchem according to CID.'''
        # Drop null values of CID
        db_act = db_act.dropna(subset=['CID'])
        cids = db_act['CID'].astype(int).unique()
        # Create a list from cids to make batch search using Pubchempy
        cids = [str(x) for x in cids if 'error' not in str(x) and x != 0]
        db_comps = pd.DataFrame(columns=columns)
        selected_columns = pc.get_columns(columns[:-2])
        if len(cids) > 0:
            compounds = pd.DataFrame()
            # Divide the list into batches of 1000 ids
            if len(cids) > 1000:
                for start, end in generate_subsets(len(cids), 1000):
                    subset = cids[start:end]
                    comps = pc.get_compounds(subset, selected_columns, namespace='cid')
                    compounds = pd.concat([compounds, comps])
                    time.sleep(0.5)
            else:
                compounds = pc.get_compounds(cids, selected_columns, namespace='cid')
            db_comps = compounds.rename(columns=pc.properties)

        return db_comps

    def _pubchem_search_inchikey(self, db_act: pd.DataFrame, columns: list[str], pc: PubChemServer) -> pd.DataFrame:
        db_act = db_act.dropna(subset=['inchikey']) # Drop null values of inchikey
        if self.merge_stereoisomers: # Extract first blocks for matching all stereoisomers
            db_act['inchikey_search'] = db_act['inchikey'].apply(pc.get_inchikey_first_block)
            inchikeys = db_act['inchikey_search'].unique()
        else:
            inchikeys = db_act['inchikey'].unique() # Use full inchikey for exact matching
        # Create a list from inchikeys to make batch search using Pubchempy
        inchikeys = [str(x) for x in inchikeys if pd.notna(x) and x != '']
        db_comps = pd.DataFrame(columns=columns)
        selected_columns = pc.get_columns(columns[:-2])
        if len(inchikeys) > 0:
            
            all_cids = []
            for inchikey in inchikeys:
                try: # Get CIDs for this inchikey (or first-block)
                    cids = pcp.get_cids(inchikey, namespace='inchikey', searchtype=None)
                    if cids:
                        all_cids.extend(cids)
                    time.sleep(0.5)
                except Exception as e:
                    pass
            
            all_cids = list(set(all_cids)) # Remove duplicates and convert to strings
            all_cids = [str(cid) for cid in all_cids if cid and cid != 0]
            compounds = pd.DataFrame()
            if len(all_cids) > 0: #batch-fetch the compound data for all CIDs
                if len(all_cids) > 1000:
                    for start, end in generate_subsets(len(all_cids), 1000):
                        subset = all_cids[start:end]
                        comps = pc.get_compounds(subset, selected_columns, namespace='cid')
                        compounds = pd.concat([compounds, comps])
                        time.sleep(0.5)
                else:
                    compounds = pc.get_compounds(all_cids, selected_columns, namespace='cid')
                if len(compounds) > 0:
                    db_comps = compounds.rename(columns=pc.properties)

        return db_comps

    def _pubchem_search_chembl(self, db_act: pd.DataFrame, id_name: str, columns: list[str], pc: PubChemServer) -> pd.DataFrame:
	'''Retrieve compounds using compound id (chembl ids).'''
        pubchem_chembl = []
        existing_id = []
        inchis = []
        chembl_ids = []
        pubchem_inchi = []
        missing_id = []
        db_comps = pd.DataFrame(columns=columns) 
        selected_columns = pc.get_columns(columns[:-2])
        # Perform search of chembl id
        for ids in db_act[id_name].unique():
            try:
                # Retrieve compound using pubchempy
                compound = pc.get_compounds(ids, selected_columns, namespace='name')
                # Check if exactly one compound has been found
                if len(compound) == 1:
                    pubchem_chembl.append(compound)
                    existing_id.append(ids)
                else:
                    missing_id.append(ids)
            except:
                missing_id.append(ids)
            time.sleep(0.5)

        db_chembl = pd.DataFrame()
        # Store Pubchempy output to dataframe
        if len(pubchem_chembl) > 0:
            db_chembl = pd.concat(pubchem_chembl)
            db_chembl = db_chembl.rename(columns=pc.properties)
            db_chembl['id'] = existing_id
            db_comps = pd.concat([db_comps, db_chembl])

        # Utilize chembl molecule API for failed Pubchem ids
        if len(missing_id) > 0:
            chembl.molecule_search(missing_id, inchis, chembl_ids)
        # Search the inchi for the compound that couldn't be searched by chembl id
        if len(inchis) > 0:
            existing_id = []
            for ids in zip(inchis, chembl_ids):
                try:
                    # Use Pubchempy to find compounds based on inchi
                    pubchem_inchi.append(pc.get_compounds(ids[0], selected_columns, namespace='inchi'))
                    existing_id.append(ids[1])
                except:
                    None
                time.sleep(0.5)
        # Check if at least one compound has been found using Pubchempy with inchi
        if len(pubchem_inchi) > 0:
            db_inchi = pd.concat(pubchem_inchi)
            db_inchi = db_inchi.rename(columns=pc.properties)
            db_inchi['id'] = existing_id
            db_comps = pd.concat([db_comps, db_inchi])

        return db_comps


