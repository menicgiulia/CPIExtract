import pandas as pd
import numpy as np
import time
import io
import requests
from ..servers.PubchemServer import PubChemServer
from ..servers.BiomartServer import BiomartServer
from ..data_manager import *
from ..utils.helper import generate_subsets
from .Database import Database

class Stitch(Database):
    
    def __init__(self, connection=None, database=None):
        # if not connection and not database:
        #     raise ValueError('Either SQL connection or database should be not None')
        if database is not None:
            self.data_manager = LocalManager(database)
        else:
            self.data_manager = SQLManager(connection, 'Stitch')

    def interactions(self, input_comp: pd.DataFrame, set_stereo: bool=True):
        """
        Retrieves proteins from Stitch database interacting with compound passed as input.

        Steps
        -----
        - Creates specific and unspecific identifiers from input CID.
        - Finds matches on the database, based on the set_stereo value \\
        Constraint:
            - Only experimental evidence
        - Links stitch to STRING to obtain updated protein ids and symbols.
        - Uses Biomart to obtain proteins info (modified with data from original DTC database) to return,
        searching with ensembl peptide id.
        
        Parameters
        ----------
        DTC_data : DataFrame
            Dataframe containing all DTC database info
        input_comp : DataFrame
            Dataframe of input compounds from which interacting proteins are found
        set_stereo : bool
            value to determine which kind of stereochemistry to consider for the input compound:
            - True (default) - specific stereochemistry
            - False - non-specific stereochemistry \\ 
            True ensures exact match to input_id

        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values:
            - entrez 
            - gene_type 
            - hgnc_symbol 
            - description 
            - datasource (Stitch) 
            - pchembl_value (NaN)
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all Stitch info about the input compound
        """

        columns = ['entrez','gene_type','hgnc_symbol','description','pchembl_value','datasource']
        # Create an empty DataFrame with the specified columns
        sttch_act = pd.DataFrame(columns=columns)
        # Create an empty Pubchem activity DataFrame
        sttch_raw = pd.DataFrame(columns=columns)
        # Drop null values of CID, to make sure the first line of the dataframe has the CID
        input_comp = input_comp.dropna(subset=['cid']).reset_index(drop=True)
        # Check if there are any input compounds remaining 
        if len(input_comp) > 0:
            # Find STITCH interactions for chemical by PubChem CID
            # Need to convert CID into stitch id, need a set number of leading zeroes (len = 8)
            CID = str(input_comp['cid'][0])
            # STITCH adds prefix CIDs for specific stereochemistry
            sttch1 = 'CIDs' + CID.zfill(8)
            # STITCH adds prefix CIDm for unspecific stereochemistry
            sttch2 = 'CIDm' + CID.zfill(8)
            
            # Merges all stereochemistries together
            if not set_stereo: 
                stitch_stereo = self.data_manager.retrieve_raw_data('chemical', sttch1)
                sttch_flat = self.data_manager.retrieve_raw_data('chemical', sttch2)

                sttch_raw = pd.concat([stitch_stereo,sttch_flat]).reset_index(drop=True)
            # Only uses the select input stereochemistry
            else: 
                sttch_raw = self.data_manager.retrieve_raw_data('chemical', sttch1)
        
            # Check if at least one match has been found
            if len(sttch_raw) > 0:
                # Filter STITCH to experimental evidence only
                sttch_act = sttch_raw.loc[sttch_raw['experimental'] > 0].reset_index(drop=True)
                # STITCH has no activity values so cannot calculate pChEMBL
                # The input STITCH data should already be filtered to human only interactions
                
                # Check if there is at least one interaction left
                if len(sttch_act) > 0:
                    # Linking Stitch to STRING to obtain updated protein IDs 
                    
                    # sttch_act['String_symbol'] = 'String Conversion Failure'
                    sttch_act['String_id'] = sttch_act['protein']

                    ids = list(sttch_act['protein'].unique())
                    
                    # Extract string ids using Stitch API
                    genes = self._extract_ids(ids)

                    if len(genes) > 0:
                        genes = genes[['output', 'symbol']].rename(columns={'output':'String_id',
                                                                            'symbol': 'String_symbol'})
                        sttch_act = pd.merge(sttch_act, genes, on='String_id', how='left')
                        
                    # Unify gene identifiers by Converting STITCH with biomart
                    ensembl = BiomartServer()
                    # Stitch uses ensembl peptide id
                    input_type = 'ensembl_peptide_id'
                    attributes = ['ensembl_peptide_id', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                    names = ['ensembl_peptide_id','entrez','gene_type','hgnc_symbol','description']

                    stitch_targets = pd.DataFrame(columns=names)    

                    # Removes the first 5 characters from each id
                    sttch_act.loc[:, ['String_id']] = sttch_act['String_id'].str[5:]
            
                    input_genes = list(sttch_act['String_id'])
                    
                    stitch_targets=ensembl.subset_search(input_type, input_genes, attributes, names)

                    # For each compound, assign specific biomart column values to the ones from the original Stitch database
                    for index, row in sttch_act.iterrows():
                        S1 = stitch_targets.loc[stitch_targets['ensembl_peptide_id']==row['String_id']]
                        if len(S1) > 0:
                            sttch_act.loc[index,'entrez'] = S1['entrez'].iloc[0]
                            sttch_act.loc[index,'gene_type'] = S1['gene_type'].iloc[0]
                            sttch_act.loc[index,'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                            sttch_act.loc[index,'description'] = S1['description'].iloc[0]
                        else:
                            sttch_act.loc[index,'entrez'] = None
                            sttch_act.loc[index,'gene_type'] = None
                            sttch_act.loc[index,'hgnc_symbol'] = None
                            sttch_act.loc[index,'description'] = None
                            Des = 'Failed to convert the gene ID'                   
                    sttch_act['datasource'] = 'Stitch'
                    sttch_act['pchembl_value'] = np.nan
                    statement = 'completed'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input compound doesn\'t contain CID'
        return sttch_act, statement, sttch_raw


    def compounds(self, input_protein: pd.DataFrame, set_stereo: bool=True):
        """
        Retrieves compounds from Stitch database interacting with proteins passed as input.

        Steps
        -----
        - Filters stitch database to find compounds interacting with input proteins: \\
        Constraints:
            - stereochemistry = specific
        - Uses Pubchempy to obtain compounds info to return, searching with CID.
        
        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compound are found
        stitch_data : DataFrame
            Dataframe containing all stitch database info        
        set_stereo : bool
            value to determine which kind of stereochemistry to consider for the input compound:
            - True (default) - specific stereochemistry
            - False - non-specific stereochemistry \\ 
            True ensures exact match to input_id

        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values:
            - inchi
            - inchikey
            - isomeric smiles
            - iupac name
            - datasource (Stitch)
            - pchembl value (NaN)
            - notes (NaN)
        """

        columns = ['inchi','inchikey','isomeric_smiles','iupac_name','datasource','pchembl_value']
        # Create an empty DataFrame with the specified columns
        stitch_c1 = pd.DataFrame(columns=columns)
        stitch_raw = pd.DataFrame()
        # Drop duplicated and null values of uniprot and ensembl peptide id
        input_protein = input_protein.dropna(subset=['uniprot', 'ensembl_peptide_id'])
        input_protein = input_protein.drop_duplicates(subset='uniprot')
        input_protein = input_protein.drop_duplicates(subset='ensembl_peptide_id').reset_index(drop=True)
        # Check if there are any input proteins remaining   
        if len(input_protein) > 0:
            input_protein_id = input_protein['ensembl_peptide_id'][0]
            input_protein_id = '9606.' + str(input_protein_id)
            # Remove homo sapiens identifier from the protein id column
            stitch_raw = self.data_manager.retrieve_raw_data('protein', input_protein_id)
            stitch_raw['protein'] = stitch_raw['protein'].str.replace('9606.', '', regex=True)
            # stitch_raw = stitch_data.loc[stitch_data['protein']==input_protein_id]

            if len(stitch_raw) > 0:
                if set_stereo:
                    # Select only specific (stereochemistry) compounds, not non specific
                    stitch_act = stitch_raw.loc[stitch_raw['chemical'].str.contains("s", case=False, na=False)].copy()
                # Remove CIDs and CIDm strings from chemical column
                stitch_act['CID'] = stitch_act['chemical'].str.replace(r'^(CIDs|CIDm)', '', regex=True).str.lstrip('0')
                stitch_act = stitch_act.dropna(subset=['CID'])
                # Select only non empty cids and experimental values > 0
                stitch_act = stitch_act.loc[(~stitch_act['CID'].eq('')) & 
                                        (stitch_act['experimental'] > 0)].reset_index(drop=True)
                # Check if there are any compounds remaining
                if len(stitch_act) > 0:

                    pc = PubChemServer()
                    # Retrieve compounds using cids
                    stitch_c1 = self._pubchem_search_cid(stitch_act, columns, pc)

                    if len(stitch_c1) > 0:
                        stitch_c1.loc[:, 'pchembl_value'] = np.nan        
                        stitch_c1.loc[:, 'datasource'] = 'Stitch'
                        stitch_c1.loc[:, 'notes'] = np.nan
                        statement = 'completed'
                    else:
                        statement = 'Compounds not found using PubChem'
                else:
                    statement = 'Filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein doesn\'t contain ensembl id'  

        return stitch_c1, statement, stitch_raw

    
    def _extract_ids(self, ids):

        names = ['input','value','output','taxid','species','symbol','description']
        params={
                # Protein list
                "identifiers" : None, 
                # species NCBI identifier (Human is being used)
                "species" : 9606, 
                # Only one (best) identifier per input protein
                "limit" : 1, 
                # See input identifiers in the output
                "echo_query" : 1, 
                # App name
                "caller_identity" : "NEU" 
            }
        # Construct URL and Call STRING
        string_api_url="https://version-11-5.string-db.org/api"
        output_format="tsv-no-header"
        method="get_string_ids"
        request_url = '/'.join([string_api_url,output_format,method])

        if len(ids) > 250:
            genes = pd.DataFrame(columns=names)
            for start, end in generate_subsets(len(ids), 250):
                subset = ids[start:end]
                params['identifiers'] = "\r".join([str(gene) for gene in subset])
                try:
                    response = requests.post(request_url, data=params)
                    # Read and parse the results
                    st = pd.read_csv(io.StringIO(response.text), sep='\t', header=None, 
                                        names=['input','value','output','taxid','species','symbol','description'])
                    genes = pd.concat([genes, st])
                except:
                    continue
                time.sleep(0.5)
        else:
            params['identifiers'] = "\r".join([str(gene) for gene in ids])
            try:
                response = requests.post(request_url, data=params)
                genes = pd.read_csv(io.StringIO(response.text), sep='\t', header=None, 
                            names=['input','value','output','taxid','species','symbol','description'])
            except:
                genes = pd.DataFrame(columns=names)
        genes['output'] = genes['output'].astype(str)
        return genes
