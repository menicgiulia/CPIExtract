import pandas as pd
import numpy as np
from ..servers.BiomartServer import BiomartServer
from ..servers.ChEMBLServer import ChEMBLServer as chembl
from ..servers.PubchemServer import PubChemServer
from .Database import Database
from ..data_manager import *


class DTC(Database):

    def __init__(self, connection=None, database=None):
        # if not connection and not database:
        #     raise ValueError('Either SQL connection or database should be not None')
        if database is not None:
            self.data_manager = LocalManager(database)
        else:
            self.data_manager = SQLManager(connection, 'DTC')


    def _filter_database(self, DTC_raw: pd.DataFrame):
        # Filter DTC_data to measurements and units for calculation of pChEMBL
        valid_standard_types = ['IC50','KI','EC50','KD','Kd','ED50','XC50','AC50','Ki','ED50','CC50','AVERAGEIC50',
                            'ACTIVITYEC50','SC50','AD50','DC50','XI50','GC50','GI50','LOGEC50','-LOGIC50','PIC50',
                            'LOGIC50','LOGKD','LOGKI','LOGKI','logIC50','LOGIC50','-LOGKD','PEC50','LOGEC50',
                            'PIC50(CALC)','LOG(10^6/IC50)','-LOGED50','LOG1/IC50','PED50','LOGKD','LOG1/KD',
                            'LOG(1/IC50)','1/KI','PXC50']
        valid_unit_types = ['NM','M','/NM','NMOL/L','MG.KG-1','UG.ML-1','UMOL.KG-1','P.P.M.','PPM',
                        'UG KG-1','MG KG-1','MG/ML','NG ML-1',"10'-4UMOL/L",'UG ML-1','MM','NMOL/MG',
                        "10'13NM","10'-10M","10'8NM",'10^5/M','MOL/G',"10'-8M",'NM G-1','M-1',"10'-5M",
                        "10'7NM",'/UM',"10'-7M","10'-9M","10'-9L/MOL","10'-10L/MOL",'PM','UM L-1','UG/ML']
        
                            # Keep only valid standard types
        DTC_filt = DTC_raw.loc[DTC_raw['standard_type'].isin(valid_standard_types) & 
                            # Keep only valid standard units
                            DTC_raw['standard_units'].isin(valid_unit_types) &
                            # Remove entries with no measurement
                            DTC_raw['standard_value'].notnull() &
                            # Remove entries without target info
                            DTC_raw['target_id'].notnull() &
                            # Remove entries without compound id
                            DTC_raw['compound_id'].notnull()
                            ].drop_duplicates().copy().reset_index(drop=True)  
        return DTC_filt
        

    def _standardize_database(self, DTC_raw: pd.DataFrame, dtc_mutated: bool, pChEMBL_thres: float):

                    # Only use interactions with reported sources
        DTC_act = DTC_raw[(DTC_raw['doc_type'].notnull()) & 
                        # Remove inconclusive and undetermined activities
                        (~DTC_raw['activity_comment'].isin(['inconclusive', 'Not Active']))].reset_index(drop=True)
        if not dtc_mutated:
            # Remove all mutated data
            DTC_act = DTC_act[(DTC_act['wildtype_or_mutant'] != 'mutated') & 
                        (DTC_act['mutation_info'].isnull())].reset_index(drop=True)
        # Converts all the measured values into nM for the calculation of pChEMBL
        DTC_act = self._standard_converter(DTC_act)
        # Convert infinite values to nan
        DTC_act['standardized_val'].replace([np.inf, -np.inf], np.nan, inplace=True)
        # Remove empty values
        DTC_act = DTC_act[DTC_act['standardized_val'].notnull()]
        # Remove non-numeric characters from the Activity Value (things like > and <)
        #DTC_act['standardized_val']=DTC_act['standardized_val'].str.replace(r'[^0-9\.]','',regex=True)
        # Convert strings into numeric values 
        DTC_act['standardized_val'] = DTC_act['standardized_val'].apply(pd.to_numeric,errors='coerce')
        # Generate pChEMBL values
        DTC_act['pchembl_value'] = -np.log10(DTC_act['standardized_val']/1e9)  
        # Filter interactions by pChEMBL value
        DTC_act=DTC_act.loc[DTC_act['pchembl_value'] > pChEMBL_thres].reset_index(drop=True) 

        return DTC_act
            

    def interactions(self, input_comp: pd.DataFrame, chembl_ids: list, dtc_mutated: bool=False, pChEMBL_thres: float=0):
        """
        Retrieves proteins from DTC database interacting with compound passed as input.

        Steps
        -----
        - Filters input database from unnecessary data \\
        Constraints:
            - Only valid standard types 
            - Only valid standard units
            - Entries with measurement only
            - Entiries with target info only \\
        - Finds matches on the database for all compounds ids passed as input and applies a second filter \\
        Constraints:
            - Only interactions with reported sources
            - No inconclusive and undetermined activities
            - No mutated data \\
        - Applies a standard conversion for all values required to compute pChembl, removing non valid values.
        - Uses Biomart to obtain proteins info (modified with data from original DTC database) to return,
        searching with uniprotswissprot.
        
        Parameters
        ----------
        DTC_data : DataFrame
            Dataframe containing all DTC database info
        input_comp : DataFrame
            Dataframe of input compounds from which interacting proteins are found
        chembl_ids : list
            list of chembl ids for the input compounds
        dtc_mutated: bool
            bool to select whether to included mutated targets interactions
        pChEMBL_thres : float
            minimum pChEMBL value necessary for interaction to be considered valid
        
        Returns
        -------
        DataFrame
            Dataframe of interacting proteins, containing the following values: \\
            entrez, gene_type, hgnc_symbol, description, datasource (DTC), pchembl_value (nan)
        String
            A statement string describing the outcome of the database search
        DataFrame
            Raw Dataframe containing all DTC info about the input compound
        """

        columns = ['entrez','gene_type','hgnc_symbol','description','pchembl_value','datasource']
        # Create an empty DataFrame with the specified columns
        DTC_act = pd.DataFrame(columns=columns)
        DTC_raw = pd.DataFrame(columns=columns)
        
        # Check if chembl_ids has already been computed from the input compound, otherwise do so
        if not chembl_ids:
            chembl_ids.extend(chembl.identify_chembl_ids(input_comp))
        
        if len(chembl_ids) > 0 and None not in chembl_ids:
            # Filter DTC by compound of interest using ChEMBL IDs
            DTC_raw = self.data_manager.retrieve_raw_data('compound_id', chembl_ids)
            # DTC_raw = DTC_data.loc[DTC_data['compound_id'].isin(chembl_ids)].reset_index(drop=True)
            
            # Check if at least one match has been found
            if len(DTC_raw) > 0:

                DTC_filt = self._filter_database(DTC_raw)

                DTC_filt['molecular_weight'] = input_comp['molecular_weight'].iloc[0]

                # Filter database
                DTC_act = self._standardize_database(DTC_filt, dtc_mutated, pChEMBL_thres)
                
                # Note: DTC has no tax_id or species column, can only filter for human via biomart conversion
                if len(DTC_act) > 0:
                    # Unify gene identifiers by Converting DTC with biomart
                    ensembl = BiomartServer()
                    attributes = ['uniprotswissprot', 'entrezgene_id', 'gene_biotype', 'hgnc_symbol', 'description']
                    names = ['uniprot','entrez','gene_type','hgnc_symbol','description']
                    # Search by uniprot id match to biomart
                    input_type = 'uniprotswissprot'
                    
                    dtc_targets = pd.DataFrame(columns=names)
                    
                    input_genes = list(DTC_act['target_id'])
                    
                    dtc_targets = ensembl.subset_search(input_type, input_genes, attributes, names)

                    # For each compound, assign specific biomart column values to the ones from the original DTC database
                    for index, row in DTC_act.iterrows():
                        S1 = dtc_targets.loc[dtc_targets['uniprot']==row['target_id']]
                        if len(S1) > 0:
                            DTC_act.loc[index,'entrez'] = S1['entrez'].iloc[0]
                            DTC_act.loc[index,'gene_type'] = S1['gene_type'].iloc[0]
                            DTC_act.loc[index,'hgnc_symbol'] = S1['hgnc_symbol'].iloc[0]
                            DTC_act.loc[index,'description'] = S1['description'].iloc[0]
                        else:
                            DTC_act.loc[index, 'entrez'] = None
                            DTC_act.loc[index, 'gene_type'] = None
                            DTC_act.loc[index, 'hgnc_symbol'] = None
                            DTC_act.loc[index, 'description'] = None  
                            Des='Failed to convert the gene ID'
                    DTC_act['datasource']='DTC'
                    statement='completed'
                else:
                    statement='Filter reduced interactions to 0'
            else:
                statement='No interaction data'
        else:
            statement = 'No ChEMBL ids found'

        return DTC_act, statement, DTC_raw
    

    def compounds(self, input_protein: pd.DataFrame, dtc_mutated: bool=False, pChEMBL_thres: float=0):
        """
        Retrieves compounds from DTC database interacting with proteins passed as input.

        Steps
        -----
        - Retrieves gene data from DTC database to find compounds interacting with gene
        Constraints:
            - Standard type and units not null
            - pChEMBL value > threshold
        - Uses Pubchempy to retrieve the compound info (with additional data from chembl) to return, 
        searching with compound ID or with inchi (obtained from chembl API) if the first search is unsuccessful.
        - Applies a standard conversion for all values required to compute pChembl.
        
        Parameters
        ----------
        input_protein : DataFrame
            Dataframe of input proteins from which interacting compound are found
        DTC_data : DataFrame
            Dataframe containing all DTC database info
        dtc_mutated: bool
            bool to select whether to included mutated targets interactions
        pChEMBL_thres : integer
            Threshold for pChEMBL value to be valid

        Returns
        -------
        DataFrame
            Dataframe of interacting compounds, containing the following values: \\
            inchi, inchikey, isomeric_smiles, iupac_name, datasource (DTC), pchembl_value, notes (nan)
        """

        columns = ['inchi','inchikey','isomeric_smiles','iupac_name','molecular_weight','datasource','pchembl_value']
        # Create an empty DataFrame with the specified columns
        DTC_c1 = pd.DataFrame(columns=columns) 
        DTC_raw = pd.DataFrame()
        # Drop duplicated and null values of uniprot
        input_protein = input_protein.dropna(subset=['uniprot'])
        input_protein = input_protein.drop_duplicates(subset='uniprot').reset_index(drop=True)  
        # Check if there are any input proteins remaining      
        if len(input_protein) > 0:
            input_protein_id = input_protein['uniprot'][0]

            # Filter DTC by protein of interest using uniprot ID
            DTC_raw = self.data_manager.retrieve_raw_data('target_id', input_protein_id)
            # Check if there are any compounds remaining
            if len(DTC_raw) > 0:

                DTC_filt = self._filter_database(DTC_raw)

                if len(DTC_filt) > 0:

                    pc = PubChemServer()
                    DTC_c = pd.DataFrame()
                    
                    if 'CID' in DTC_filt.columns:
                        # Retrieve compounds using cids
                        compounds = self._pubchem_search_cid(DTC_filt, columns, pc)

                        if len(compounds) > 0:
                            # Add additional information to interacting compounds and standardize its values       
                            DTC_c = pd.merge(compounds, DTC_filt, on='CID', how='left')
                    else:
                        # Retrieve compounds using compound id (chembl ids) and inchis if formers miss
                        compounds = self._pubchem_search_chembl(DTC_filt, 'compound_id', columns, pc)

                        if len(compounds) > 0:
                            # Add additional information to interacting compounds and standardize its values       
                            DTC_c = pd.merge(compounds, DTC_filt, left_on='id', right_on='compound_id', how='left')
                    
                    # Check if there are any compounds remaining
                    if len(DTC_c) > 0:
                        # Standardize database
                        DTC_act = self._standardize_database(DTC_c, dtc_mutated, pChEMBL_thres)

                        if len(DTC_act) > 0:
                            # Extend compounds information
                            DTC_c1 = DTC_act[columns[:-2]+['CID', 'compound_id', 'standard_type', 'standard_value', 'standard_units', 
                                                'pchembl_value', 'activity_comment']].drop_duplicates()
                            # Add additional information
                            DTC_c1['notes'] = np.nan
                            DTC_c1['datasource'] = 'DTC'
                            statement = 'completed'
                        else:
                            statement = 'Standardization reduced interactions to 0'
                    else:
                        statement = 'Compounds not found using Pubchem and ChEMBL'
                else:
                    statement = 'First filter reduced interactions to 0'
            else:
                statement = 'No interaction data'
        else:
            statement = 'Input protein doesn\'t contain uniprot id'                         
        
        return DTC_c1, statement, DTC_raw


    def _standard_converter(self, DTC_act: pd.DataFrame):
        """
        Converts all types from DTC database to those needed for pCHEMBL by storing the respective
        standardized value. Also converts all units to nM as it is required to compute pCHEMBL.

        Parameters
        ----------
        DTC_act : DataFrame
            DTC dataframe with compounds whose types and units need to be standardized

        Returns
        -------
        DataFrame
            DTC dataframe of compounds with standardized types and unit
        """    

        std_types = ['IC50','KI','EC50','KD','Kd','ED50','XC50','AC50','Ki','ED50','CC50','AVERAGEIC50',
                'ACTIVITYEC50','SC50','AD50','DC50','XI50','GC50','GI50']
        log_types = ['LOGEC50','LOGIC50','LOGKD','LOGKI','logIC50','LOGIC50','LOGEC50','LOGKD']
        nlog_types = ['-LOGIC50','PIC50','-LOGKD','PEC50','PIC50(CALC)','-LOGED50','PED50','PXC50']
        replog_types = ['LOG1/IC50','LOG1/KD','LOG(1/IC50)']
        
        # Converting standard types
        DTC_act['standardized_val'] = None
        D1 = DTC_act.loc[DTC_act['standard_type'].isin(std_types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = DTC_act['standard_value'][h]
        D1 = DTC_act.loc[DTC_act['standard_type'].isin(nlog_types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = 10**(-DTC_act['standard_value'][h])
        D1 = DTC_act.loc[DTC_act['standard_type'].isin(log_types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = 10**(DTC_act['standard_value'][h])
        D1 = DTC_act.loc[DTC_act['standard_type']=='LOG(10^6/IC50)']
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = 10**(6-(DTC_act['standard_value'][h]))
        D1 = DTC_act.loc[DTC_act['standard_type'].isin(replog_types)] 
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = 10**(-DTC_act['standard_value'][h])
        D1 = DTC_act.loc[DTC_act['standard_type']=='1/KI'] 
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = (1/DTC_act['standard_value'][h])
        
        # Converting all standard units to nM
        types = ['NM','NMOL/L']
        D1 = DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]
        types = ['M']
        D1 = DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = DTC_act['standardized_val'][h]*(1e9)
        types = ['/NM']
        D1 = DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val'] = 1/DTC_act['standardized_val'][h]
        types = ['MM']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e9/1e3)
        types=['10^5/M']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=((1/(DTC_act['standardized_val'][h]*(1e5)))*1e9)
        types=['NM G-1']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e3)
        types=['M-1']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=((1/DTC_act['standardized_val'][h])*1e9)
        types=['/UM']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=(1/(DTC_act['standardized_val'][h]))*(1e9/1e6)
        types=['PM']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e9/1e12)
        types=['MOL/G']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e3/1)*(1e9/1)
        types=['NMOL/MG']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e3/1e-3)
        types=['UM L-1','UMOL.KG-1']  
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e9/1e6)
        
        types=['UG.ML-1','UG ML-1','UG/ML']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            if D1['molecular_weight'][h]=='':
                DTC_act.loc[h,'standardized_val']=np.nan
            else:
                DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e3/1)*(1/1e6)*(1/float(D1['molecular_weight'][h]))*(1e9/1)
        types=['MG/ML']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            if D1['molecular_weight'][h]=='':
                DTC_act.loc[h,'standardized_val']=np.nan
            else:
                DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e3/1)*(1/1e3)*(1/float(D1['molecular_weight'][h]))*(1e9/1)
        types=['NG ML-1']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            if D1['molecular_weight'][h]=='':
                DTC_act.loc[h,'standardized_val']=np.nan
            else:
                DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1e3/1)*(1/1e9)*(1/float(D1['molecular_weight'][h]))*(1e9/1)
        types=['MG.KG-1','MG KG-1']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            if D1['molecular_weight'][h]=='':
                DTC_act.loc[h,'standardized_val']=np.nan
            else:
                DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1/1e3)*(1/float(D1['molecular_weight'][h]))*(1e9/1)
        types=['UG KG-1']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            if D1['molecular_weight'][h]=='':
                DTC_act.loc[h,'standardized_val']=np.nan
            else:
                DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1/1e6)*(1/float(D1['molecular_weight'][h]))*(1e9/1)
        types=['P.P.M.','PPM']
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            if D1['molecular_weight'][h]=='':
                DTC_act.loc[h,'standardized_val']=np.nan
            else:
                DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(1/float(D1['molecular_weight'][h]))*(1/1e3)*(1e9/1)
        types=["10'13NM","10'8NM","10'7NM"]
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            E=DTC_act['standard_units'][h][DTC_act['standard_units'][h].find("'")+1:DTC_act['standard_units'][h].find("N")]
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(10**int(E))
        types=["10'-10M","10'-8M","10'-5M","10'-7M","10'-9M"]
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            E=DTC_act['standard_units'][h][DTC_act['standard_units'][h].find("'")+1:DTC_act['standard_units'][h].find("M")]
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(10**int(E))*(1e9)
        types=["10'-4UMOL/L"]
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            E=DTC_act['standard_units'][h][DTC_act['standard_units'][h].find("'")+1:DTC_act['standard_units'][h].find("U")]
            DTC_act.loc[h,'standardized_val']=DTC_act['standardized_val'][h]*(10**int(E))*(1e9/1e6)
        types=["10'-9L/MOL","10'-10L/MOL"]
        D1=DTC_act.loc[DTC_act['standard_units'].isin(types)]
        for h in D1.index:
            E=DTC_act['standard_units'][h][DTC_act['standard_units'][h].find("'")+1:DTC_act['standard_units'][h].find("L")]
            DTC_act.loc[h,'standardized_val']=(1/(DTC_act['standardized_val'][h]*(10**int(E))))*(1e9)
        
        return DTC_act
