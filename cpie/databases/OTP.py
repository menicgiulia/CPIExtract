import pandas as pd
import numpy as np
import requests
import time
import json
from ..servers.BiomartServer import BiomartServer
from ..servers.ChEMBLServer import ChEMBLServer as chembl
from ..servers.PubchemServer import PubChemServer
from ..data_manager import *
from .Database import Database

class OTP(Database):

    def __init__(self, connection=None, database=None):
        if connection is not None or database is not None:
            raise ValueError('No SQL connection or database expected for OTP')
        funcs = {
            'bib_comps' : self._retrieve_bib_compounds,
            'met_comps' : self._retrieve_met_compounds,
            'proteins' : self._retrieve_proteins 
        }
        self.data_manager = APIManager(funcs)

    def interactions(self, input_comp: pd.DataFrame, chembl_ids: list, use_biblio: bool=False):
        """
        Retrieves proteins from OTP API interacting with compound passed as input.

        Steps
        -----
        - Retrieves Mechanism and Bibliography interactions, based on the use_biblio value, using the chembl
        ids passed as input.
        - Uses Biomart to obtain proteins info (modified with data from first API searchs) to return,
        searching with ensembl peptide id.
        
        Parameters
        ----------
        input_comp : DataFrame
            Dataframe of input compounds from which interacting proteins are found
        chembl_ids : list
            list of chembl ids for the input compounds
        use_biblio : bool
            value to determine which kind of interactions to retrieve:
            - False (default) - mechanism data only 
            - True - mechanism and bibliography \\ 
            False ensures real interactions, avoids keyword search matches within abstracts that can lead to false associations that are not real interactions
            
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
        # Create an empty DataFrame with the specified columns
        otp_act=pd.DataFrame(columns=columns)
        # Create an OTP activity DataFrame
        raw_columns = ['Chembl','PMID','mID','mLabel','mType','strength']
        otp_raw = pd.DataFrame(columns=raw_columns)
        
        # Check if chembl_ids has already been computed from the input compound, otherwise do so
        if not chembl_ids:
            chembl_ids.extend(chembl.identify_chembl_ids(input_comp))

        if len(chembl_ids) > 0 and None not in chembl_ids:
            # If True then takes both mechanisms and bibliography sections of OTP
            if use_biblio: 
                otp_bib = self.data_manager.retrieve_raw_data('bib_comps', chembl_ids)
            # Else only Mechanisms is collected from OTP, skipping bibliography
            otp_met = self.data_manager.retrieve_raw_data('met_comps', chembl_ids)
            
            # Merges bibliography and mechanism findings only if use_biblio is True
            if use_biblio:
                otp_raw = pd.concat([otp_bib,otp_met])
            else:
                otp_raw = otp_met

            if len(otp_raw) > 0:  
                # Filter interactions
                # Remove interactions without at least 1 source
                otp_filtered = otp_raw.loc[(otp_raw['PMID'].notnull()) &  
                                    # Removes interactions that are not gene interactions
                                    (otp_raw['mType']=='GP')].reset_index(drop=True) 
                # Note: OTP is human only, so no animal data here
                # Note: OTP does not collect activity data so we cannot calculate pChEMBL
                
                # Check if there is at least one interaction remaining 
                if len(otp_filtered) > 0:
                    otp_act = pd.concat([otp_act, otp_filtered])
                    # Unify gene identifiers by Converting OTP with biomart
                    ensembl = BiomartServer()
                    attributes = ['ensembl_gene_id', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                    names = ['ensembl_id','entrez','gene_type','hgnc_symbol','description']
                    # Search by ensembl peptide id match to biomart
                    input_type='ensembl_gene_id'
                    
                    otp_targets = pd.DataFrame(columns=names)
                    
                    input_genes = list(otp_act['mID'])
                    
                    otp_targets = ensembl.subset_search(input_type, input_genes, attributes, names)
                    
                    # For each compound, assign specific biomart column values to the ones from the original OTP database
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
                            Des='Failed to convert the gene ID'
                    otp_act['datasource'] = 'OTP'
                    otp_act['pchembl_value'] = np.nan
                    statement='completed'
                else:
                    statement='Filter reduced interactions to 0'
            else:
                statement='No interaction data'
        else:
            statement = 'No ChEMBL ids found'
            
        return otp_act, statement, otp_raw


    def compounds(self, input_protein: pd.DataFrame):
        """
        Retrieves compounds from OTP database interacting with proteins passed as input.

        Steps
        -----
        - Retrieves gene data from OTP with its API to find compounds interacting with gene \\
        Constraints:
            - None
        - Uses Pubchempy to retrieve the compound info (with additional data from chembl) to return, 
        searching with DrugId or with inchi (obtained from chembl API) if the first search is unsuccessful.
        
        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compound are found

        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values: \\
            inchi, inchikey, isomeric_smiles, iupac_name, datasource (OTP), pchembl_value, notes (nan)
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
            # Using ensembl gene id to search gene in otp
            input_protein_id = input_protein['ensembl_gene_id'][0]
            
            otp_raw = self.data_manager.retrieve_raw_data('proteins', input_protein_id)

            pc = PubChemServer()
            if len(otp_raw) > 0:    
                # Perform search using drugId and inchis if formers miss
                otp_c1 = self._pubchem_search_chembl(otp_raw, 'drugId', columns, pc)

                if len(otp_c1) > 0:
                    # Store into dataframe add additional data
                    otp_c1['datasource'] = 'OTP'
                    otp_c1['pchembl_value'] = np.nan
                    otp_c1['notes'] = np.nan
                    statement = 'completed'
                else:
                    statement = 'Compounds not found using Pubchem and ChEMBL'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein doesn\'t contain ensembl id'

        return otp_c1, statement, otp_raw
    

    def _retrieve_bib_compounds(self, chembl_ids):
        # Set base URL of GraphQL API endpoint
        base_url = "https://api.platform.opentargets.org/api/v4/graphql"
        # Note: bibliography of OTP only looks for concurrence of IDs in abstracts, it does not look for associations
        otp_bib = pd.DataFrame(columns=['Chembl','PMID','mID','mLabel','mType','strength'])
        bib_list = []

        for ChEMBLID in chembl_ids:          
            # Set variables object of arguments for GraphQL API
            variables = {"ChEMBL ID": ChEMBLID} 
            # Set the Query for the GraphQL API to get the OTP bibliography
            # This inputs the ChEMBLID variable into the query_string at the %s spot
            query_string="""
            query {
                drug(chemblId: "%s"){
                    literatureOcurrences{
                        rows{
                            pmid
                            sentences{
                                matches{
                                    mappedId
                                    matchedLabel
                                    matchedType
                                }
                            }
                        }
                    }
                }
            }
            """ % ChEMBLID 
            
            while True:
                try:
                    # Perform request
                    dat=requests.post(base_url,json={"query": query_string,"variables": variables}) 
                    # print(dat.status_code) #200 code = Operating Normally
                    # Transform API response from JSON to dictionary
                    res=json.loads(dat.text) 
                    break
                except json.JSONDecodeError as e:
                    print(f"Error: {e}. Retrying in 10 seconds...")
                    time.sleep(10)
            
            # Check if at least one interaction has been found
            if res['data']['drug'] is not None:
                    for row in res['data']['drug']['literatureOcurrences']['rows']:
                        for sentence in row['sentences']:
                            for match in sentence['matches']:
                                try:
                                    # Add protein to DataFrame
                                    bib_list.append({
                                        'Chembl': ChEMBLID,
                                        'PMID': row['pmid'],
                                        'mID': match['mappedId'],
                                        'mLabel': match['matchedLabel'],
                                        'mType': match['matchedType'],
                                        'strength': 'weak'
                                    })
                                except:
                                    continue
        # Create DataFrame
        otp_bib = pd.concat([otp_bib, pd.DataFrame(bib_list)])
        otp_bib = otp_bib.drop_duplicates().reset_index(drop=True)

        return otp_bib
    
    def _retrieve_met_compounds(self, chembl_ids):
        # Set base URL of GraphQL API endpoint
        base_url = "https://api.platform.opentargets.org/api/v4/graphql"

        otp_met = pd.DataFrame(columns=['Chembl','PMID','mID','mLabel','mType','strength'])
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
                try:
                    # Perform request
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
                        # Add protein to DataFrame
                        try:
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
    
    def _retrieve_proteins(self, input_protein_id):
        # Create post request
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
        try:
            # Save relevant API data into dataframe
            otp_raw = pd.DataFrame(res['data']['target']['knownDrugs']['rows'])
            # Remove duplicates of drugId
            otp_raw = otp_raw.drop_duplicates(subset='drugId').reset_index(drop=True)
        except TypeError:
            # There was an error with the data types, return empty dataframe
            otp_raw = pd.DataFrame(columns=['inchi','inchikey','isomeric_smiles','iupac_name','datasource','pchembl_value'])

        return otp_raw

            

