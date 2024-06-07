import requests
import time
import pandas as pd
import numpy as np
from ..servers.BiomartServer import BiomartServer
from ..servers.PubchemServer import PubChemServer
from ..data_manager import *
from ..utils.helper import generate_subsets
from .Database import Database

class PubChem(Database):

    def __init__(self, connection=None, database=None):
        if connection is not None or database is not None:
            raise ValueError('No SQL connection or database expected for Pubchem')
        funcs = {
            'comps' : self._retrieve_compounds,
            'targets' : self._retrieve_targets,
            'proteins' : self._retrieve_proteins
        }
        self.data_manager = APIManager(funcs)

    def _filter_database(self, pubchem_raw: pd.DataFrame, identifier: str, pChEMBL_thres: float):
        
        # Drop duplicates, invalid and empty values of identifier and Activity Values
        pubchem_filt = pubchem_raw.dropna(subset=[identifier, 'Activity Value [uM]'])
        pubchem_filt = pubchem_filt[(~pubchem_filt[identifier].eq('')) & (~pubchem_filt['Activity Value [uM]'].eq(''))]
        # Remove Unspecified, Inactive and Inconclusive activities
        pubchem_filt = pubchem_filt[(~pubchem_filt['Activity Outcome'].isin(['Unspecified', 'Inactive', 'Inconclusive'])) &
                                    # Limit measures to Ki, Kd, IC50, EC50, Potency             
                                    (pubchem_filt['Activity Name']!='CD')]
        # Remove duplicates with same identifier and activity value
        pubchem_filt = pubchem_filt.drop_duplicates(subset=[identifier, 'Activity Value [uM]']).copy()
        # Remove non-numeric characters from the Activity Value (things like > and <)
        pubchem_filt['Activity Value [uM]'] = pubchem_filt['Activity Value [uM]'].str.replace(r'[^0-9\.]','',regex=True)
        # Convert strings into numeric values 
        pubchem_filt['Activity Value [uM]'] = pubchem_filt['Activity Value [uM]'].apply(pd.to_numeric,errors='coerce')
        # Generate pChEMBL values 
        pubchem_filt['pchembl_value'] = -np.log10((pubchem_filt['Activity Value [uM]']/1e6)) 
        # Filter interactions by pChEMBL value
        pubchem_filt=pubchem_filt.loc[pubchem_filt['pchembl_value'] > pChEMBL_thres].reset_index(drop=True) 

        return pubchem_filt

    def interactions(self, input_comp: pd.DataFrame, pChEMBL_thres: float=0):
        """
        Retrieves proteins from PubChem API interacting with compound passed as input.

        Steps
        -----
        - Finds information of input compound using PubChem API and filters the resulting DataFrame \\
        Constraints:
            - Only entries with geneids
            - No Unspecified, Inactive or Inconclusive activities
            - Only Ki, Kd, IC50, EC50, Potency measures
        - Computes the pChembl value and identifies target proteins' taxonomy id, searching with gene ID \\
        Constraint:
            - Only Homo Sapiens genes
        - Uses Biomart to obtain proteins info (modified with data from first API search) to return,
        searching with Entrez ID.
        
        Parameters
        ----------
        input_comp : DataFrame
            Dataframe of input compounds from which interacting proteins are found
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid

        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values:
            - entrez 
            - gene_type 
            - hgnc_symbol 
            - description 
            - datasource (PubChem) 
            - pchembl_value

        String
            A statement string describing the outcome of the database search 
        DataFrame
            Raw Dataframe containing all PubChem info about the input compound
        """
        
        columns = ['entrez','gene_type','hgnc_symbol','description','pchembl_value','datasource']
        # Create an empty DataFrame with the specified columns
        pubchem_act=pd.DataFrame(columns=columns)
        # Create an empty Pubchem activity DataFrame
        pubchem_raw=pd.DataFrame()
        # Drop null values of CID, to make sure the first line of the dataframe has the CID
        input_comp = input_comp.dropna(subset=['cid']).reset_index(drop=True)
        # Check if there are any input compounds remaining  
        if len(input_comp) > 0:
            # Select compound from Pubchem based on CID
            input_comp_id = input_comp['cid'][0]
            
            pubchem_raw = self.data_manager.retrieve_raw_data('comps', input_comp_id)

            if len(pubchem_raw) > 0:
                # Filter Pubchem for high quality interactions and species
                pubchem_act = self._filter_database(pubchem_raw, 'Target GeneID', pChEMBL_thres)
                
                # Find remaining genes
                gene_list = list(pubchem_act['Target GeneID'].unique()) 
            
                gene_tax = self.data_manager.retrieve_raw_data('targets', gene_list)

                # Merge tax_id information into pubchem_act DataFrame
                pubchem_act = pubchem_act.merge(gene_tax, left_on='Target GeneID', right_on='GeneID', how='left')

                # Filter for rows where tax_id is 9606 (human)
                pubchem_act = pubchem_act.loc[pubchem_act['TaxonomyID']==9606].reset_index(drop=True) 
                
                # Check if there is at least one interaction remaining
                if len(pubchem_act) > 0:
                    # Unify gene identifiers by Converting PubChem with biomart
                    ensembl = BiomartServer()
                    # PubChem uses Entrez IDs
                    input_type = 'entrezgene_id' 
                    attributes = ['entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                    names = ['entrez','gene_type','hgnc_symbol','description']

                    pubchem_targets = pd.DataFrame(columns=names)
                    
                    genelist = list(pubchem_act['Target GeneID'])
                    # list of genes that need conversion
                    input_genes = [int(i) for i in genelist] 

                    pubchem_targets = ensembl.subset_search(input_type, input_genes, attributes, names)

                    # For each compound, assign specific biomart column values to the ones from the original CTD database    
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
                            Des='Failed to convert the gene ID'
                    pubchem_act['datasource']='PubChem'
                    statement = 'completed'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input compound doesn\'t contain CID'

        return pubchem_act, statement, pubchem_raw
    

    def compounds(self, input_protein: pd.DataFrame, pChEMBL_thres: float=0):
        """
        Retrieves compounds from PubChem database interacting with proteins passed as input.

        Steps
        -----
        - Retrieves gene data from PubChem PUG REST API to find compounds interacting with gene: \\
        Constraints:
            - Activity value = Active
        - Uses Pubchempy to obtain compounds info to return, searching with CID and merges information together.
        
        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compound are found
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid

        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values:
            - inchi
            - inchikey
            - isomeric smiles
            - iupac name
            - datasource (PubChem)
            - pchembl value (NaN)
            - notes (NaN)
        """

        columns = ['inchi','inchikey','isomeric_smiles','iupac_name','datasource','pchembl_value']
        # Create an empty DataFrame with the specified columns
        pubchem_c1 = pd.DataFrame(columns=columns)
        # Create dataframe from pubchem information
        pubchem_raw = pd.DataFrame()
        # Drop duplicated and null values of entrez
        input_protein = input_protein.dropna(subset=['entrez'])
        input_protein = input_protein.drop_duplicates(subset='entrez').reset_index(drop=True)
        # Check if there are any input proteins remaining         
        if len(input_protein) > 0:
            # Using entrez to search gene in pubchem
            input_protein_id = input_protein['entrez'][0]

            pubchem_raw = self.data_manager.retrieve_raw_data('proteins', input_protein_id)

            if len(pubchem_raw) > 0:
                # Filter Pubchem for high quality interactions and species
                pubchem_act = self._filter_database(pubchem_raw, 'CID', pChEMBL_thres)

                # Check if there are any compounds remaining
                if len(pubchem_act):
                    pc = PubChemServer()

                    # Retrieve compounds using cids
                    compounds = self._pubchem_search_cid(pubchem_act, columns, pc)

                    if len(compounds) > 0:
                        pubchem_act['CID'] = pubchem_act['CID'].astype(int)
                        pubchem_info = pubchem_act[['CID', 'Target Accession', 'pchembl_value']].\
                                        rename(columns={'Target Accession': 'Target Uniprot'}).drop_duplicates()
                        
                        # Add additional values from original request
                        pubchem_c1 = pd.merge(compounds, pubchem_info, on='CID', how='left')
                        pubchem_c1['datasource'] = 'PubChem'
                        pubchem_c1['notes'] = np.nan
                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using PubChem'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein doesn\'t contain entrez'
            
        return pubchem_c1, statement, pubchem_raw
    
    
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
            # Check if a match has been found
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
    
    def _retrieve_targets(self, gene_list):
        # Use the PUG REST Pubchem API to assign the tax_id for each gene
        url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/gene/geneid/summary/JSON'
        # Prepare the request payload
        if len(gene_list) > 1000:
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
    
    def _retrieve_proteins(self, input_protein_id):

        pubchem_raw = pd.DataFrame()
        # Use PUG REST to get the information of gene activity from PubChem
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/gene/geneid/%s/concise/JSON' %(str(input_protein_id))
        response=requests.get(url) 
        # Retrieve and filter only relevant data
        data=response.json()
        # Check if a match has been found
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
            