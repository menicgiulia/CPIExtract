### Data Download

To operate the package, and try the examples it is necessary to have downloaded databases `.csv` and `.tsv` files either in a local directory or stored in a mySQL server as tables. The compressed preprocessed data necessary to reproduce the results in the manuscript are provided in subdirectory `data` - `Databases.zip`. This data is intended for use with the [Local Execution mode](#local-execution). Please make sure to download the `.zip` file separately and not with the whole GitHub directory.

2025 Update

The databases are now unified by InChIKey as well as CID. This allows CPIExtract to collect interactions at a purely structural or a deeper stereochemical level by match to the first block of InChIKey only or matching to the full InChIKey. Such matching especially helps when the different datasets assign their interactions to a default stereochemsitry; for example flat stereochemistry. 

2026 Update

CPIExtract has been updated to include a variety of options. Users can 1) create their own compound id dataframe inputs to bypass PubChem, if say the compound is not on PubChem, 2) set Stitch option to set experimental threshold level of confidence, 3) set protein type for DrugBank curation, and 4) set pchembl_grouping options to either average pchembl across standard types or average per standard type. CPIExtract now outputs both an aggregated dataframe containing the average pchembl per compound-protein interaction pair and a total dataframe containing the filtered interactions per database without any averaging for traceability. The output now classifies interactions based on pChEMBL value using default thresholds commonly used in pharmaceutical industry (these can be user defined if needed as well). Stereochemistry is now locked to scaffold matching with an output column "structure match" that indicates if the match is a stereochemical match or a scaffold match. For package stability, retry logic has been applied to all API queries 

The data for the 2026 update are available at:
[https://www.dropbox.com/scl/fi/pnknnogqf72lluduavkit/Databases.zip?rlkey=d8upmtl8cjmo7tooohm84l4np&st=r5b39h4v&dl=0]

Create a local folder, name it `data` and extract the files with the following code, which will save the databases in `data/Databases` subdirectory(Replace `/user_path_to/data/` with the appropriate path of the package in your local/remote machine):

  _On Linux/Mac_:
   
  ```bash
  cd /user_path_to/data/
  mkdir -p Databases
  unzip Databases.zip -d Databases/
  ```
      
  _On Windows shell/Powershell_:

  ```bash
  cd user_path_to\\data
  mkdir Databases
  tar -xf Databases.zip -C Databases\\
  ```

#### Data Update

Although the data used is up to date, each database periodically releases updated versions that will make the zipped data obsolete. 
For this reason, we suggest periodically redownloading the databases to have the latest CPI information available.
We also strongly recommend preprocessing the databases to obtain significantly faster execution times. \
To ease this process, we provide a notebook to download and preprocess the data ([db_preprocessing.ipynb](db_preprocessing.ipynb)).

### SQL Server

Users can also load the databases into a mySQL server. 

1. We suggest using [mySQL Workbench CE](https://dev.mysql.com/downloads/workbench/). Once set up, create a database schema named `cpie`.
2. (If the package has been installed with pip) \
Download the [`dbs_config.json`](data/dbs_config.json) and save it into the data folder.
3. Load the **downloaded and preprocessed** databases as tables into the SQL database, which can be done using the provided [SQL_load.ipynb](SQL_load.ipynb) notebook.
