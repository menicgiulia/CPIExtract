# C-PIE (Compound-Protein Interaction Extract)
## A software package to collect and harmonize small molecule and protein interactions
##### Authors: Andrea Piras, Shi Chenghao, Michael Sebek, Giulia Menichetti (giulia.menichetti@channing.harvard.edu)
## Introduction
The binding interactions between small molecules and proteins are the basis of cellular functions.
Yet, experimental data available regarding compound-protein interaction (CPI) is not harmonized into a single
entity but rather scattered across multiple institutions, each maintaining databases with different formats.
Extracting information from these multiple sources remains challenging due to data heterogeneity.

C-PIE interactively extracts, filters and harmonizes CPI data from 9 databases, providing the output in tabular format (csv).\
The package provides two separate pipelines:
- Comp2Prot: Extract protein target interactions for compounds provided as input
- Prot2Comp: Extract compound interactions for proteins provided as input
### Comp2Prot
The pipeline extracts compound information from an input identifier using the [PubChem REST API](https://pubchempy.readthedocs.io/en/latest/), accessed with the [`PubChemPy`](https://github.com/mcs07/PubChemPy) Python package. Then, it uses the information to perform compound matching for each database, extracting raw interaction data. This data is then filtered for each database with a custom filter that ensures only high-quality interactions are returned in the output. Finally, proteins data are extracted using the [`Biomart`](https://github.com/sebriois/biomart) Python package, harmonizing the information from the 9 databases. The collected output is finally returned to the user in a `.csv` file. \
An exemplative pipeline workflow is depicted in the figure below. An equivalent output network is also shown.

![Comp2Prot pipeline workflow example](/images/pipeline.png)
### Prot2Comp
The pipeline follows a workflow similar to the other pipeline. It extracts protein information from an input identifier using the Biomart Python package. It performs the same database matching and filtering, and then harmonizes the interaction data extracting compounds information with the PubChemPy Python package. The output is finally returned to the user in a `.csv` file. 
## Getting started 

### Setting up a work environment

#### I. Without installing the package

1. Ensure you have Python installed.
  
2. Copy the project to your local or remote machine:

   ```bash
      git clone https://github.com/menicgiulia/C-PIE.git
   ```
3. Navigate to the project directory:

   ```bash
      cd C-PIE
   ```

4. Installing the necessary dependencies:
   
   
##### Option A: working with Conda

Working with Conda is recommended, but it is not essential. If you choose to work with Conda, these are the steps you need to take:

- Ensure you have Conda installed.

- Create a new conda environment with the `environment.yml` file:

      conda env create -f environment.yml

- Activate your new conda environment:

      conda activate cpie
  
##### Option B: working without Conda

- Ensure the following dependencies are installed before proceeding:

      pip install numpy pandas mysql-connector-python biomart pubchempy chembl-webresource-client  
            
5. Set up your PYTHONPATH (Replace `/user_path_to/C-PIE/cpie` with the appropriate path of the package in your local/remote machine.):

    _On Linux/Mac_:
   
      ```bash
      export PYTHONPATH="/user_path_to/C-PIE/cpie":$PYTHONPATH
      ```
      
    _On Windows shell_:

    ```bash
      set PYTHONPATH="C:\\user_path_to\\C-PIE\\cpie";%PYTHONPATH%
    ```
      
    _On Powershell_:

    ```bash
      $env:PYTHONPATH = "C:\\user_path_to\\C-PIE\\cpie;" + $env:PYTHONPATH
    ```
    
      
#### II. With installing the package

1. Installing the necessary dependencies:
   
   
##### Option A: working with Conda

Working with Conda is recommended, but it is not essential. If you choose to work with Conda, these are the steps you need to take:

- Ensure you have Conda installed.

- Download the environment.yml and navigate to the directory of your local/remote machine where the file is located.

- Create a new conda environment with the `environment.yml` file:

      conda env create -f environment.yml

- Activate your new conda environment:

      conda activate cpie
  
##### Option B: working without Conda

- Ensure the following dependencies are installed before proceeding:

      pip install numpy pandas mysql-connector-python biomart pubchempy chembl-webresource-client 

2. Install the package:

      ```bash
      pip install c-pie
      ```
#### Using Jupyter Notebooks
We provide several Jupyer Notebooks to simplify the databases download and mantenaince. In order to use these notebooks, follow these steps:
- Make sure you have the `Jupyter` package installed.

  ```bash
      pip install jupyter
  ```
  
- Start the Jupyter Kernel

    a) If you are working on a local machine:
  
    ```bash
      jupyter notebook --browser="browser_of_choice"
    ```
    
  Replace browser_of_choice with your preferred browser (e.g., chrome, firefox). The browser window should pop up automatically. If it doesn't, copy and paste the link provided in the terminal into your browser. The link should look something like this:

    
   * http://localhost:8889/tree?token=5d4ebdddaf6cb1be76fd95c4dde891f24fd941da909129e6
    
       
    b) If you are working on a remote machine:
  
    ```bash
      jupyter notebook --no-browser
    ```
    
  Then copy and paste the link provided in the terminal in your local browser of choice; it should look something like this:

    
   * http://localhost:8888/?token=9feac8ff1d5ba3a86cf8c4309f4988e7db95f42d28fd7772
    
    
- Navigate to the selected notebook in the Jupyter Notebook interface and start executing the cells.

### Data Download
To operate the package, it is necessary to have downloaded databases `.csv` and `.tsv` files either in a data folder or stored in a mySQL server as tables. \
We provide the compressed preprocessed data used to obtain the results reported in the paper in the data folder. \
This data is to be used with the [Local Execution mode](#local-execution). 

Extract the files with the following code:
```bash
# replace root path with appropriate location (e.g., /home/username/)
cd root/path/C-PIE/data/
# unzip file
unzip Databases.zip
```
#### Data Update
Although the data used is up to date, each database periodically releases updated versions that will make the zipped data obsolete. 
For this reason, we suggest periodically redownloading the databases in order to have the latest CPI information available.
We also strongly recommend preprocessing the databases to obtain significantly faster pipelines execution times. \
To ease this process, we provide two notebooks to download ([db_download.ipynb](db_download.ipynb)) and preprocess ([db_preprocessing.ipynb](db_preprocessing.ipynb)) the databases.

### SQL Server
Users can also load the databases into a mySQL server. \
We suggest using [mySQL Workbench CE](https://dev.mysql.com/downloads/workbench/). Once setup, create a database schema named `cpie` and load the **downloaded and preprocessed** databases as tables into the SQL database, which can be done using the provided [SQL_load.ipynb](SQL_load.ipynb) notebook.

## Example execution
### Pipelines instantiation
The package provides two classes, one for each pipeline: Comp2Prot and Prot2Comp. 
The pipelines can be instanciated in two ways:
#### Local execution
Load locally the databases when executing the pipeline. \
To load the databases, follow the example notebook [Comp2Prot_example.ipynb](Comp2Prot_example.ipynb) at cell 2 (Load in Required Datasets).
```python 
from cpie.pipelines.Comp2Prot import Comp2Prot
from cpie.pipelines.Prot2Comp import Prot2Comp
import pandas as pd
import os

# Root data path
data_path = 'data'

#Downloaded from BindingDB on 3/30/2023
file_path=os.path.join(data_path, 'BindingDB.csv')
BDB_data=pd.read_csv(file_path,sep=',',usecols=['CID', 'Ligand SMILES','Ligand InChI','BindingDB MonomerID','Ligand InChI Key','BindingDB Ligand Name','Target Name Assigned by Curator or DataSource','Target Source Organism According to Curator or DataSource','Ki (nM)','IC50 (nM)','Kd (nM)','EC50 (nM)','pH','Temp (C)','Curation/DataSource','UniProt (SwissProt) Entry Name of Target Chain','UniProt (SwissProt) Primary ID of Target Chain'],on_bad_lines='skip')

#Downloaded from STITCH on 2/22/2023
file_path=os.path.join(data_path, 'STITCH.tsv')
sttch_data=pd.read_csv(file_path,sep='\t')

#Downloaded from ChEMBL on 2/01/2024
file_path=os.path.join(data_path, 'ChEMBL.csv')
chembl_data=pd.read_csv(file_path,sep=',')

file_path=os.path.join(data_path, 'CTD.csv')
CTD_data=pd.read_csv(file_path,sep=',')

#Downloaded from DTC on 2/24/2023
file_path=os.path.join(data_path, 'DTC.csv')
DTC_data=pd.read_csv(file_path,sep=',',usecols=['CID', 'compound_id','standard_inchi_key','target_id','gene_names','wildtype_or_mutant','mutation_info','standard_type','standard_relation','standard_value','standard_units','activity_comment','pubmed_id','doc_type'])

#Downloaded from DrugBank on 3/2/2022
file_path=os.path.join(data_path, 'DB.csv')
DB_data=pd.read_csv(file_path, sep=',')

#Downloaded from DrugCentral on 2/25/2024
file_path=os.path.join(data_path, 'DrugCentral.csv')
DC_data=pd.read_csv(file_path, sep=',')

# Data stored in pandas dataframes
data = {
    'chembl': chembl_data,
    'bdb': BDB_data,
    'stitch': sttch_data,
    'ctd': CTD_data,
    'dtc': DTC_data,
    'db': DB_data,
    'dc': DC_data
}

C2P = Comp2Prot(execution_mode='local', dbs=data)
P2C = Prot2Comp(execution_mode='local', dbs=data)
```
#### Server execution
The databases are stored on the mySQL server. 

```python
from cpie.pipelines import Comp2Prot
from cpie.pipelines import Prot2Comp

# Exmeplative dictionary containing the configuration info to connect to the mySQL server
info = {
    "host": "XXX.XX.XX.XXX",
    "user": "user",
    "password": "password",
    "database": "cpie",
}

C2P = Comp2Prot(execution_mode='server', server_info=info)
P2C = Prot2Comp(execution_mode='server', server_info=info)
```
### How to run the pipelines
Once instantiated, to run the pipelines call either `comp_interactions` or `comp_interactions_select` for Comp2Prot, and either `prot_interactions` or `prot_interactions_select` for Prot2Comp. 
#### Comp2Prot
Comp2Prot accepts the following parameters:
- input_id: the compound id
- pChEMBL_thresh: the minimum interaction pChEMBL value required to be added to the output file
- stitch_stereo: to select whether to consider the specific compound stereochemistry or group all stereoisomers interactions from STITCH
- otp_biblio: to select whether to include the *bibliography* data from OTP. This parameter is only available for Comp2Prot as OTP provides only known drug interactions for proteins
- dtc_mutated: to select whether also to consider interactions with mutated target proteins from DTC
- dc_extra: to select whether to include possibly non-Homo sapiens interactions

The output will include all the interactions found and a dataframe containing the statements for all the datasets for the specific input compound.
```python
# Chlorpromazine InChIKey
comp_id = 'ZPEIMTDSQAKGNT-UHFFFAOYSA-N'

interactions, db_states = C2P.comp_interactions(input_id=comp_id, pChEMBL_thres=0, stitch_stereo=True, otp_biblio=False, dtc_mutated=False, dc_extra=False)
```
To extract interactions only from selected databases, use the alternate function specifying which databases to use in a underscore separated string (to include all databases, which equates to using the previous function, use `'pc_chembl_bdb_stitch_ctd_dtc_otp_dc_db'`).
```python
# Chlorpromazine InChIKey
comp_id = 'ZPEIMTDSQAKGNT-UHFFFAOYSA-N'

# Interactions extracted from PubChem, ChEMBL, DB and DTC only.
interactions, db_states = C2P.comp_interactions_select(input_id=comp_id, selected_dbs='pc_chembl_db_dtc', pChEMBL_thres=0, stitch_stereo=True, otp_biblio=False, dtc_mutated=False, dc_extra=False)
```
#### Prot2Comp
Prot2Comp works similarly, with the exception that both functions only have the following additional parameters, as OTP will retrieve only known compounds interacting with the input protein:
- pChEMBL_thresh: the minimum interaction pChEMBL value required to be added to the output file
- stitch_stereo: to select whether to consider the specific compound stereochemistry or group all stereoisomers interactions from STITCH
- dtc_mutated: to select whether also to consider interactions with mutated target proteins from DTC
- dc_extra: to select whether to include possibly non-Homo sapiens interactions
```python
# HGNC symbol for Kallikrein-1
prot_id = 'KLK1'

interactions, db_states = P2C.prot_interactions(input_id=prot_id, pChEMBL_thres=0, stitch_stereo=True, dtc_mutated=False, dc_extra=False)

# Interactions extracted from PubChem, ChEMBL, DB and DTC only.
interactions, db_states = P2C.prot_interactions_select(input_id=prot_id, selected_dbs='pc_chembl_db_dtc', pChEMBL_thres=0, stitch_stereo=True, dtc_mutated=False, dc_extra=False)
```
## Package Structure
Root folder organization (```__init__.py``` files removed for simplicity):
```plaintext
│   .gitignore
│   environment.yml                                 // current conda env settings used
│   README.md
│   CPE_testing.ipynb                               // CPE pipeline testing notebook
│   PCE_testing.ipynb                               // PCE pipeline testing notebook
│
├───data                                            // database storage location
│   └───input                                       // pipeline input data location
│
└───cpie
    │   
    ├───data_manager                                
    │   ├───APIManager.py                           // to extract raw data from databases' APIs
    │   ├───DataManager.py                          // Abstract data manager class
    │   ├───LocalManager.py                         // to extract raw data from downloaded databases
    │   └───SQLManager.py                           // to extract raw data from SQL server
    │
    ├───databases
    │   ├───BindingDB.py                            // BDB database class
    │   ├───ChEMBL.py                               // ChEMBLB database class
    │   ├───CTD.py                                  // CTD database class
    │   ├───Database.py                             // Abstract database class
    │   ├───DrugBank.py                             // DB database class
    │   ├───DrugCentral.py                          // DC database class
    │   ├───DTC.py                                  // DTC database class
    │   ├───OTP.py                                  // OTP database class
    │   ├───PubChem.py                              // PubChem database class
    │   └───Stitch.py                               // STITCH database class
    │
    ├───pipelines                                  
    │   ├───CompProtExtract.py                      // CPE pipeline class
    │   ├───Pipeline                                // Abstract pipeline class
    │   └───ProtCompExtract.py                      // PCE pipeline class
    │
    ├───servers                                     
    │   ├───BiomartServer.py                        // to connect to Biomart API
    │   ├───ChEMBLServer.py                         // to connect to ChEMBL API
    │   └───PubChemServer.py                        // to connect to PubChem API
    │
    ├───sql_server  
    │   ├───dbs_config.json                         // file with databases configuration for the mySQL server                                 
    │   └───sql_conncetion.py                       // to connect to the SQL server
    │
    └───utils                                        
        ├───helper.py                               // helper functions and classes
        └───identifiers.py                          // functions to extract identifiers for compounds and proteins
```
## Further information

- Details about each function (what is it used for, what are the input parameters, the possible values of the input parameters, what is the output) from the pipeline are available in the `cpie` folder, in the comments before each class function. 
- An example on the use of the implemented functions is available in the notebooks [Comp2Prot_example.ipynb](Comp2Prot_example.ipynb) and [Prot2Comp_example.ipynb](Prot2Comp_example.ipynb), which can be executed in order to test the proper installation of the package and it's functionalities.

## License

This project is licensed under the terms of the MIT license.