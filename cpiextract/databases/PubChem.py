'''Loading,searching,filtering and preprocessing data from PubChem.'''

import requests
import time
import pandas as pd
import numpy as np
import pubchempy as pcp
import duckdb
import os

from ..utils.typing import Connection
from ..servers.BiomartServer import BiomartServer
from ..servers.PubChemServer import PubChemServer
from ..data_manager import *
from ..utils.helper import generate_subsets
from .Database import Database

class PubChem(Database):

    def __init__(self, connection: Connection| None=None, database: pd.DataFrame|None=None, merge_stereoisomers=False,
                 bioact_file=None):
        super().__init__(merge_stereoisomers)
        if connection is not None or database is not None:
            raise ValueError('No SQL connection or database expected for Pubchem')
    
        self.db_file = None
        self.use_local = False
    
        if bioact_file:
            # Check if it's already a .duckdb file
            if bioact_file.endswith('.duckdb'):
                self.db_file = bioact_file
            else:
                # Derive .duckdb path from .tsv.gz path
                self.db_file = bioact_file.replace('pc_bioactivities.tsv.gz', 'pubchem.duckdb')
        
            self.use_local = os.path.exists(self.db_file)
    
        funcs = {
            'comps' : self._retrieve_compounds,
            'targets' : self._retrieve_targets,
            'proteins' : self._retrieve_proteins
        }
        self.data_manager = APIManager(funcs)

    def _filter_database(self, pubchem_raw: pd.DataFrame, identifier: str, pChEMBL_thres: float) -> pd.DataFrame:
    
        # Handle different column names between API and bulk file
        activity_col = 'Activity Value' if 'Activity Value' in pubchem_raw.columns else 'Activity Value [uM]'
    
        # Drop duplicates, invalid and empty values of identifier and Activity Values
        pubchem_filt = pubchem_raw.dropna(subset=[identifier, activity_col])
        pubchem_filt = pubchem_filt[(~pubchem_filt[identifier].eq('')) & (~pubchem_filt[activity_col].eq(''))]
        # Remove Unspecified, Inactive and Inconclusive activities
        pubchem_filt = pubchem_filt[(~pubchem_filt['Activity Outcome'].isin(['Unspecified', 'Inactive', 'Inconclusive'])) &
                                    # Limit measures to Ki, Kd, IC50, EC50, Potency             
                                    (pubchem_filt['Activity Name']!='CD')]
        # Remove duplicates with same identifier and activity value
        pubchem_filt = pubchem_filt.drop_duplicates(subset=[identifier, activity_col]).copy()
        # Remove non-numeric characters from the Activity Value (things like > and <)
        pubchem_filt[activity_col] = pubchem_filt[activity_col].astype(str).str.replace(r'[^0-9\.]','',regex=True)
        # Convert strings into numeric values 
        pubchem_filt[activity_col] = pubchem_filt[activity_col].apply(pd.to_numeric,errors='coerce')
    
        # Handle Activity Unit - API has it in column name, bulk file has separate column
        if 'Activity Unit' in pubchem_filt.columns:
            # Filter for uM only (bulk file)
            pubchem_filt = pubchem_filt[pubchem_filt['Activity Unit'] == 'uM']
        # else: API data already has [uM] in column name, assume all values are uM
    
        # Generate pChEMBL values (value is already in uM)
        pubchem_filt = pubchem_filt[pubchem_filt[activity_col] > 0].copy() # Filter out zero values first to avoid log10(0)
        pubchem_filt['pchembl_value'] = -np.log10((pubchem_filt[activity_col]/1e6))
        # Filter interactions by pChEMBL value
        pubchem_filt = pubchem_filt.loc[pubchem_filt['pchembl_value'] > pChEMBL_thres].reset_index(drop=True) 

        return pubchem_filt

    def interactions(self, input_comp: pd.DataFrame, pChEMBL_thres: float=0, merge_stereoisomers: bool=False, verbose: bool=False) -> tuple[pd.DataFrame, str, pd.DataFrame]:
        """
        Retrieves proteins from PubChem interacting with compound passed as input.
        Uses local database when available, API otherwise.
        """

        # Print database info if verbose
        if verbose:
            if self.use_local:
                print(f"Using PubChem database: {self.db_file}")
            else:
                print("Using PubChem API")

        columns = ['entrez','gene_type','hgnc_symbol','description','pchembl_value','datasource','inchikey']
        pubchem_act = pd.DataFrame(columns=columns)
        pubchem_raw = pd.DataFrame()
    
        # Determine CID list
        if merge_stereoisomers:
            input_comp = input_comp.dropna(subset=['inchikey']).reset_index(drop=True)
            if len(input_comp) == 0:
                return pubchem_act, 'Input compound does not contain inchikey', pubchem_raw
        
            firstblock = input_comp['inchikey_fb'][0]
        
            # Try local database first
            if self.use_local:  # CHANGED: removed "and self.cid_inchikey_file"
                cid_list = self._get_cids_from_firstblock(firstblock)
                if cid_list is None:  # Database error, use API
                    try:
                        cid_list = pcp.get_cids(firstblock, namespace='inchikey', searchtype=None)
                        if not cid_list:
                            return pubchem_act, 'No CIDs found for first block inchikey', pubchem_raw
                    except Exception as e:
                        return pubchem_act, f'Error retrieving CIDs: {str(e)}', pubchem_raw
            else:  # No local database, use API
                try:
                    cid_list = pcp.get_cids(firstblock, namespace='inchikey', searchtype=None)
                    if not cid_list:
                        return pubchem_act, 'No CIDs found for first block inchikey', pubchem_raw
                except Exception as e:
                    return pubchem_act, f'Error retrieving CIDs: {str(e)}', pubchem_raw
        else:  # Use single CID
            input_comp = input_comp.dropna(subset=['cid']).reset_index(drop=True)
            if len(input_comp) > 0:
                cid_list = [int(input_comp['cid'][0])]
            else:
                return pubchem_act, 'Input compound does not contain CID', pubchem_raw
        
        # Retrieve bioactivities - local or API
        if self.use_local:
            # Use local files with DuckDB
            pubchem_act = self._retrieve_from_local(cid_list, pChEMBL_thres)
            
            if len(pubchem_act) == 0:
                return pubchem_act, 'No interaction data in local files', pubchem_raw
            
            statement = 'completed (local)'
            
        else:
            # Use API
            raw_cid_data = []
            for cid in cid_list:
                try:
                    raw_data = self.data_manager.retrieve_raw_data('comps', cid)
                    if len(raw_data) > 0:
                        raw_cid_data.append(raw_data)
                    time.sleep(0.5)
                except Exception as e:
                    continue
            
            if len(raw_cid_data) == 0:
                return pubchem_act, 'Failed to retrieve data from PubChem API', pubchem_raw
            
            pubchem_raw = pd.concat(raw_cid_data, ignore_index=True)
            
            if len(pubchem_raw) == 0:
                return pubchem_act, 'No interaction data', pubchem_raw
            
            # Filter database
            pubchem_act = self._filter_database(pubchem_raw, 'Target GeneID', pChEMBL_thres)
            
            # Get taxonomy via API
            gene_list = list(pubchem_act['Target GeneID'].unique())
            gene_tax = self.data_manager.retrieve_raw_data('targets', gene_list)
            pubchem_act = pubchem_act.merge(gene_tax, left_on='Target GeneID', right_on='GeneID', how='left')
            pubchem_act = pubchem_act.loc[pubchem_act['TaxonomyID']==9606].reset_index(drop=True)
            
            if len(pubchem_act) == 0:
                return pubchem_act, 'Filter reduced interactions to 0', pubchem_raw
            
            statement = 'completed (API)'
        
        # Get InChIKeys for CIDs (common to both paths)
        unique_cids = pubchem_act['CID'].dropna().unique().tolist()
        if len(unique_cids) > 0:
            pc = PubChemServer()
            inchikey_columns = ['inchikey']
            selected_columns = pc.get_columns(inchikey_columns)

            try:
                cid_inchikey_map = pd.DataFrame()

                if len(unique_cids) > 1000:
                    for start, end in generate_subsets(len(unique_cids), 1000):
                        subset = [str(cid) for cid in unique_cids[start:end]]
                        batch_result = pc.get_compounds(subset, selected_columns, namespace='cid')
                        cid_inchikey_map = pd.concat([cid_inchikey_map, batch_result])
                        time.sleep(0.5)
                else:
                    cid_list_str = [str(cid) for cid in unique_cids]
                    cid_inchikey_map = pc.get_compounds(cid_list_str, selected_columns, namespace='cid')

                if len(cid_inchikey_map) > 0:
                    cid_inchikey_map = cid_inchikey_map.rename(columns=pc.properties)
                    pubchem_act['CID'] = pubchem_act['CID'].astype(str)
                    cid_inchikey_map['CID'] = cid_inchikey_map['CID'].astype(str)
                    pubchem_act = pubchem_act.merge(cid_inchikey_map[['CID', 'inchikey']], 
                                                    on='CID', how='left')
                else:
                    pubchem_act['inchikey'] = None
            except Exception as e:
                pubchem_act['inchikey'] = None

        # Get gene annotations from Biomart
        ensembl = BiomartServer()
        input_type = 'entrezgene_id' 
        attributes = ['entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
        names = ['entrez','gene_type','hgnc_symbol','description']

        genelist = list(pubchem_act['Target GeneID'])
        input_genes = [int(i) for i in genelist]

        pubchem_targets = ensembl.subset_search(input_type, input_genes, attributes, names)

        # Merge biomart data
        for index, row in pubchem_act.iterrows():
            S1 = pubchem_targets.loc[pubchem_targets['entrez']==int(row['Target GeneID'])]
            if len(S1) > 0:
                pubchem_act.loc[index,'entrez'] = S1['entrez'].iloc[0]
                pubchem_act.loc[index,'gene_type'] = S1['gene_type'].iloc[0]
                pubchem_act.loc[index,'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                pubchem_act.loc[index,'description'] = S1['description'].iloc[0]
            else:
                pubchem_act.loc[index,'entrez'] = None
                pubchem_act.loc[index,'gene_type'] = None
                pubchem_act.loc[index,'hgnc_symbol'] = None
                pubchem_act.loc[index,'description'] = None
        
        pubchem_act['datasource'] = 'PubChem'

        return pubchem_act, statement, pubchem_raw

    # Keep existing _retrieve_compounds, _retrieve_targets, _retrieve_proteins methods
    def _retrieve_compounds(self, input_comp_id):
        pubchem_raw=pd.DataFrame()
        if input_comp_id == '':
            Des='No Data'
        else:
            # Use pug rest to get the information from PubChem
            url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/assaysummary/JSON'%(str(input_comp_id))
            response=requests.get(url)
            time.sleep(0.5) 
            data=response.json()
            try:
                Des=data['Table']
            except:
                Des='No Data'
        
        if Des!='No Data':
            # Convert API data into a dataframe
            raw_columns = data['Table']['Columns']['Column']
            rows = [row['Cell'] for row in data['Table']['Row']]
            pubchem_raw = pd.DataFrame(rows, columns=raw_columns)

        return pubchem_raw
    
    def _retrieve_targets(self, gene_list) -> pd.DataFrame:
        # Use the PUG REST Pubchem API to assign the tax_id for each gene
        url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/gene/geneid/summary/JSON'
        if len(gene_list) > 1000: # Prepare the request payload
            gene_tax = pd.DataFrame(columns=['GeneID', 'TaxonomyID'])
            for start, end in generate_subsets(len(gene_list), 1000):
                subset = gene_list[start:end]        
                payload = {"geneid": ",".join(map(str, subset))}
                try:
                    response=requests.post(url, data=payload)
                    data = response.json()
                    gene_summary_list = data['GeneSummaries']['GeneSummary']
                    genes = pd.DataFrame(gene_summary_list, columns=['GeneID', 'TaxonomyID'])
                    gene_tax = pd.concat([gene_tax, genes])
                except:
                    continue
                time.sleep(0.5)
        else:
            payload = {"geneid": ",".join(map(str, gene_list))}
            try:
                response=requests.post(url, data=payload)
                data = response.json()
                gene_summary_list = data['GeneSummaries']['GeneSummary']
                gene_tax = pd.DataFrame(gene_summary_list, columns=['GeneID', 'TaxonomyID'])
            except:
                gene_tax = pd.DataFrame(columns=['GeneID', 'TaxonomyID'])
        gene_tax['GeneID'] = gene_tax['GeneID'].astype(str)
        return gene_tax
    
    def _retrieve_proteins(self, input_protein_id) -> pd.DataFrame:
        pubchem_raw = pd.DataFrame()
        # Use PUG REST to get the information of gene activity from PubChem
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/gene/geneid/%s/concise/JSON' %(str(input_protein_id))
        response=requests.get(url) 
        # Retrieve and filter only relevant data
        data=response.json()
        try:
            Des=data['Table']
        except:
            Des='No Data'

        if Des != 'No Data':
            cols=data['Table']['Columns']['Column']
            rows=[]
            for row in data['Table']['Row']:
                rows.append(row['Cell'])
            # Create dataframe from pubchem information
            pubchem_raw = pd.DataFrame(rows, columns=cols)

        return pubchem_raw

    def compounds(self, input_protein: pd.DataFrame, pChEMBL_thres: float=0, merge_stereoisomers: bool=False, verbose: bool=False):
        """
        Retrieves compounds from PubChem interacting with proteins passed as input.
        Returns basic compound info - additional metadata can be added in postprocessing.
        """
        
        # Print database info if verbose
        if verbose:
            if self.use_local:
                print(f"Using PubChem database: {self.db_file}")
            else:
                print("Using PubChem API")

        columns = ['inchi','inchikey','isomeric_smiles','iupac_name','datasource','pchembl_value']
        pubchem_c1 = pd.DataFrame(columns=columns)
        pubchem_raw = pd.DataFrame()
        
        # Drop duplicated and null values of entrez
        input_protein = input_protein.dropna(subset=['entrez'])
        input_protein = input_protein.drop_duplicates(subset='entrez').reset_index(drop=True)
        
        if len(input_protein) == 0:
            pubchem_c1 = pubchem_c1[columns]
            return pubchem_c1, 'Input protein does not contain entrez', pubchem_raw
        
        # Get entrez ID
        input_protein_id = int(input_protein['entrez'][0])
        
        # Try local database first
        if self.use_local:
            try:
                con = duckdb.connect(self.db_file, read_only=True)
                con.execute("SET enable_progress_bar=false")
                
                # Query compounds that interact with this protein
                query = f"""
                    SELECT DISTINCT c.CID, c.InChI, c.InChIKey, 
                        b."Activity Value", b."Activity Name", b."Activity Unit"
                    FROM bioactivities b
                    INNER JOIN cid_inchikey c ON b.CID = c.CID
                    WHERE CAST(b."Gene ID" AS INTEGER) = {input_protein_id}
                    AND b."Activity Value" IS NOT NULL
                    AND b."Activity Outcome" NOT IN ('Unspecified', 'Inconclusive', 'Inactive')
                    AND b."Activity Name" != 'CD'
                    AND b."Activity Unit" = 'uM'
                """
                
                result = con.execute(query).df()
                con.close()
                
                if len(result) > 0:
                    # Rename columns
                    result = result.rename(columns={
                        'InChI': 'inchi',
                        'InChIKey': 'inchikey'
                    })
                    
                    # Calculate pChEMBL
                    result['Activity Value'] = result['Activity Value'].astype(str).str.replace(r'[^0-9\.]','',regex=True)
                    result['Activity Value'] = result['Activity Value'].apply(pd.to_numeric, errors='coerce')
                    result['pchembl_value'] = -np.log10((result['Activity Value']/1e6))
                    
                    # Filter by threshold
                    result = result[result['pchembl_value'] > pChEMBL_thres].copy()
                    
                    if len(result) > 0:
                        # Create output with all required columns
                        pubchem_c1 = result[['inchi', 'inchikey', 'pchembl_value']].drop_duplicates()
                        pubchem_c1['isomeric_smiles'] = None  # Will be filled in postprocessing if needed
                        pubchem_c1['iupac_name'] = None  # Will be filled in postprocessing if needed
                        pubchem_c1['datasource'] = 'PubChem'
                        
                        # Reorder columns
                        pubchem_c1 = pubchem_c1[columns]
                        
                        statement = 'completed (local)'
                    else:
                        statement = 'Filter reduced interactions to 0'
                else:
                    statement = 'No interaction data in local database'
                    
            except Exception as e:
                if verbose:
                    print(f"Local database query failed: {e}, falling back to API")
                self.use_local = False  # Fall through to API
        
        # Use API if local failed or not available
        if not self.use_local or len(pubchem_c1) == 0:
            pubchem_raw = self.data_manager.retrieve_raw_data('proteins', input_protein_id)

            if len(pubchem_raw) > 0:
                # Filter for high quality interactions
                pubchem_act = self._filter_database(pubchem_raw, 'CID', pChEMBL_thres)

                if len(pubchem_act) > 0:
                    # Get unique CIDs
                    unique_cids = pubchem_act['CID'].dropna().unique().tolist()
                    
                    if len(unique_cids) > 0:
                        pc = PubChemServer()
                        try:
                            cid_list = [str(cid) for cid in unique_cids]
                            selected_columns = pc.get_columns(['inchi', 'inchikey', 'isomeric_smiles', 'iupac_name'])
                            
                            if len(cid_list) > 1000:
                                compounds = pd.DataFrame()
                                for start, end in generate_subsets(len(cid_list), 1000):
                                    subset = cid_list[start:end]
                                    batch = pc.get_compounds(subset, selected_columns, namespace='cid')
                                    compounds = pd.concat([compounds, batch])
                                    time.sleep(0.5)
                            else:
                                compounds = pc.get_compounds(cid_list, selected_columns, namespace='cid')
                            
                            compounds = compounds.rename(columns=pc.properties)
                            
                            if len(compounds) > 0:
                                pubchem_act['CID'] = pubchem_act['CID'].astype(int)
                                compounds['CID'] = compounds['CID'].astype(int)
                                
                                pubchem_info = pubchem_act[['CID', 'pchembl_value']].drop_duplicates()
                                
                                pubchem_c1 = pd.merge(compounds, pubchem_info, on='CID', how='left')
                                pubchem_c1['datasource'] = 'PubChem'
                                
                                # Ensure all columns exist
                                for col in columns:
                                    if col not in pubchem_c1.columns:
                                        pubchem_c1[col] = None
                                
                                pubchem_c1 = pubchem_c1[columns]
                                statement = 'completed (API)'
                            else:
                                statement = 'Compounds not found using PubChem'
                        except Exception as e:
                            statement = f'Error retrieving compound info: {str(e)}'
                    else:
                        statement = 'No valid CIDs found'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        
        pubchem_c1 = pubchem_c1[columns]
        return pubchem_c1, statement, pubchem_raw

    def _get_cids_from_firstblock(self, firstblock):
        """Get all CIDs matching a firstblock from database"""
    
        if not self.use_local:
            return None
        con = None
    
        try:
            con = duckdb.connect(self.db_file, read_only=True)
            con.execute("SET enable_progress_bar=false")

            result = con.execute(f"""
                SELECT DISTINCT CID 
                FROM cid_inchikey 
                WHERE firstblock = '{firstblock}'
            """).df()
            con.close()
        
            if len(result) == 0:
                return []
            return result['CID'].tolist()
        
        #except Exception as e:
        #    return None
        finally:
            if con is not None:
                con.close()

    def _retrieve_from_local(self, cid_list, pChEMBL_thres):
        """Retrieve bioactivities from database with SQL JOIN for human genes"""

        con=None
        try:
            con = duckdb.connect(self.db_file, read_only=True)
            con.execute("SET enable_progress_bar=false")
    
            cid_str = ','.join(map(str, cid_list))
    
            # JOIN with gene_info to filter for human genes in SQL
            pubchem_raw = con.execute(f"""
                SELECT b.*
                FROM bioactivities b
                INNER JOIN gene_info g ON CAST(b."Gene ID" AS INTEGER) = g.GeneID
                WHERE b.CID IN ({cid_str})
                AND g.TaxonomyID = 9606
                AND b."Activity Value" IS NOT NULL
                AND b."Activity Outcome" NOT IN ('Unspecified', 'Inconclusive', 'Inactive')
                AND b."Activity Name" != 'CD'
                AND b."Gene ID" IS NOT NULL
                AND b."Activity Unit" = 'uM'
            """).df()

            if len(pubchem_raw) == 0:
                return pd.DataFrame()
    
            # Rename columns to match API format
            pubchem_raw['Gene ID'] = pubchem_raw['Gene ID'].astype(int)
            pubchem_raw = pubchem_raw.rename(columns={
                'Gene ID': 'Target GeneID',
                'Protein Accession': 'Target Accession'
            })
        
            # Apply pChEMBL filtering
            pubchem_filt = self._filter_database(pubchem_raw, 'Target GeneID', pChEMBL_thres)
            return pubchem_filt
        
        #except Exception as e:
        #    return pd.DataFrame()
        finally:
            if con is not None:
                con.close()
