import numpy as np
import pandas as pd
import re
from ..servers.BiomartServer import BiomartServer
from ..servers.PubchemServer import PubChemServer
from .Database import Database
from ..data_manager import *
import time

class BindingDB(Database):

    def __init__(self, connection=None, database=None):
        # if not connection and not database:
        #     raise ValueError('Either SQL connection or database should be not None')
        if database is not None:
            self.data_manager = LocalManager(database)
        else:
            self.data_manager = SQLManager(connection, 'BDB')

    def _filter_database(self, bdb_raw: pd.DataFrame):
        """
        Filters the bdb database with the following constraints:
            - Only Homo Sapiens interactions
            - Temp(°C) < 40 or null
            - 5 < pH < 9 or null
        
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
        return bdb_filt

    def _compute_pchembl(self, bdb_dat: pd.DataFrame, pChEMBL_thres: float):
        """
        Computes the pchembl value for each non-null activity columns
        - Ki (nM)
        - IC50 (nM)
        - Kd (nM)
        - EC50 (nM)

        Parameters
        ----------
        bdb_dat : DataFrame
            bdb database
        
        Returns
        -------
        DataFrame
            database with a row for each non-null activity column and respective pchembl value obtained 
        """
        cols = ['Ki (nM)','IC50 (nM)','Kd (nM)','EC50 (nM)']
        # Collect pChEMBL values
        pchem=[] 
        # Collect the activity type used 
        act=[] 
        # Convert the activity values to numbers
        bdb_dat[cols] = bdb_dat[cols].map(lambda x: re.sub(r'[^0-9.]', '', x) if type(x) is str else x)
        bdb_dat[cols] = bdb_dat[cols].map(pd.to_numeric, errors='coerce') 
        # print(len(bdb_dat))
        # Create bdb_act dataframe for export
        bdb_act = pd.DataFrame(columns=bdb_dat.columns).astype(bdb_dat.dtypes)
        # Iterate over each activity type column
        for col in cols: 
            # Select only non null values for current activity
            bdb_val = bdb_dat.dropna(subset=[col])
            if len(bdb_val) > 0:
                bdb_act = pd.concat([bdb_act, bdb_val])
                # Compute pchembl values of entire column
                pchems = bdb_val[col].to_numpy()
                # Calculate the pChEMBL value
                pchems = (-np.log10((pchems) / 1e9))
                pchem.extend(pchems)
                # Store activity type from which pchembl has been obtained
                act.extend([col]*bdb_val[col].size)

        bdb_act.drop(cols, axis=1, inplace=True)
        bdb_act['pchembl_value'] = pchem
        bdb_act['notes'] = act

        # Filter interactions by pChEMBL value
        bdb_act=bdb_act.loc[bdb_act['pchembl_value'] > pChEMBL_thres].reset_index(drop=True) 

        return bdb_act
                


    def interactions(self, input_comp: pd.DataFrame, pChEMBL_thres: float=0):
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
        - Uses Biomart to obtain proteins info (modified with data from original bdb database) to return,
        searching with uniprotswissprot ID.
        
        Parameters
        ----------
        BDB_data : DataFrame
            Dataframe containing all bdb database info
        input_comp : dictionary
            dictionary of input compound data from which interacting proteins are found
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid
            
        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values: \\
            entrez, gene_type, hgnc_symbol, description, datasource (BindingDB), pchembl_value
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all BindingDB info about the input compound
        """

        columns = ['entrez','gene_type','hgnc_symbol','description','pchembl_value','datasource']
        # Create an empty DataFrame with the specified columns
        bdb_act = pd.DataFrame(columns=columns)
        bdb_raw = pd.DataFrame(columns=columns)
        # Drop null values of inchikey, to make sure the first line of the dataframe has the inchi key 
        input_comp = input_comp.dropna(subset=['inchikey']).reset_index(drop=True)
        # Check if there are any input compounds remaining  
        if len(input_comp) > 0:
            # Select compound from BindingDB data based on inchi key
            input_comp_id = input_comp['inchikey'][0]
            bdb_raw = self.data_manager.retrieve_raw_data('Ligand InChI Key', input_comp_id)
            # bdb_raw = BDB_data.loc[BDB_data['Ligand InChI Key']==input_comp_id].reset_index(drop=True)
            if len(bdb_raw) > 0:
                # Filter database
                bdb_filt = self._filter_database(bdb_raw)
                bdb_filt.drop(['Ligand SMILES','Ligand InChI'], axis=1, inplace=True)
                # Compute pchembl values
                bdb_act = self._compute_pchembl(bdb_filt, pChEMBL_thres)
                
                if len(bdb_act) > 0:
                    # Unify gene identifiers by Converting BindingDB with biomart
                    ensembl = BiomartServer()
                    # BindingDB uses Uniprot IDs
                    input_type='uniprotswissprot' 
                    attributes = ['uniprotswissprot', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                    names = ['uniprot','entrez','gene_type','hgnc_symbol','description']
                    bdb_targets = pd.DataFrame(columns=names)
                    
                    input_genes = list(bdb_act['UniProt (SwissProt) Primary ID of Target Chain'])
                    
                    bdb_targets = ensembl.subset_search(input_type, input_genes, attributes, names)
    
                    # For each protein, assign specific biomart column values to the ones from the original bdb database
                    for index, row in bdb_act.iterrows():
                        # Find the compound using uniprot
                        S1 = bdb_targets.loc[bdb_targets['uniprot'] == row['UniProt (SwissProt) Primary ID of Target Chain']]
                        if len(S1) > 0:
                            bdb_act.loc[index, 'entrez'] = S1['entrez'].iloc[0]
                            bdb_act.loc[index, 'gene_type'] = S1['gene_type'].iloc[0]
                            bdb_act.loc[index, 'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                            bdb_act.loc[index, 'description'] = S1['description'].iloc[0]
                        else:
                            bdb_act.loc[index, 'entrez'] = None
                            bdb_act.loc[index, 'gene_type'] = None
                            bdb_act.loc[index, 'hgnc_symbol'] = None
                            bdb_act.loc[index, 'description'] = None
                            Des='Failed to convert the gene ID'
                    bdb_act['datasource'] = 'BindingDB'
                    statement = 'completed'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input compound doesn\'t contain inchi key'

        return bdb_act, statement, bdb_raw
    

    def compounds(self, input_protein: pd.DataFrame, pChEMBL_thres: float=0):
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
        BDB_data : DataFrame
            Dataframe containing all bdb database info

        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values: \\
            inchi, inchikey, isomeric_smiles, iupac_name, datasource (BindingDB), pchembl_value, notes (activity type)
        """

        columns = ['inchi','inchikey','isomeric_smiles','iupac_name','datasource','pchembl_value']
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
            bdb_raw = self.data_manager.retrieve_raw_data('UniProt (SwissProt) Primary ID of Target Chain', input_protein_id)
            # bdb_raw = BDB_data.loc[(BDB_data['UniProt (SwissProt) Primary ID of Target Chain']==input_protein_id)]
            if len(bdb_raw):
                # Filter database
                bdb_c = self._filter_database(bdb_raw)
                
                # Compute pchembl values
                bdb_act = self._compute_pchembl(bdb_c, pChEMBL_thres)
                # Check if at least one interacting compound has been found
                if len(bdb_act) > 0:
                    pc = PubChemServer()

                    if 'CID' in bdb_act:
                        # Retrieve compounds using cids
                        compounds = self._pubchem_search_cid(bdb_act, columns, pc)

                        if len(compounds) > 0:    
                            # Add additional values from activity dataframe
                            bdb_info = bdb_act[['CID', 'UniProt (SwissProt) Primary ID of Target Chain', 
                                                'BindingDB Ligand Name', 'pchembl_value', 'notes']].\
                                        rename(columns={'UniProt (SwissProt) Primary ID of Target Chain': 'uniprotid'}).\
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
                            bdb_info = bdb_act[['Ligand InChI Key', 'UniProt (SwissProt) Primary ID of Target Chain', 
                                                'BindingDB Ligand Name', 'pchembl_value', 'notes']].\
                                        rename(columns={'Ligand InChI Key': 'inchikey', 
                                                        'UniProt (SwissProt) Primary ID of Target Chain': 'uniprotid'}).\
                                                        drop_duplicates()
                            
                            bdb_c1 = pd.merge(bdb_c1, bdb_info, on='inchikey', how='left')

                    if len(bdb_c1) > 0:
                        bdb_c1['datasource'] = 'BindingDB'
                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using PubChem'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein doesn\'t contain uniprot id'
        return bdb_c1, statement, bdb_raw
