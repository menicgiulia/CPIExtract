from .sql_connection import connect_to_mysql
import json
import os

# Get the directory where the current script is located
package_dir = os.path.dirname(__file__)

# Construct the full path to the JSON file
json_file_path = os.path.join(package_dir, 'dbs_config.json')

# Load the JSON data
with open(json_file_path, 'r') as f:
    dbs_config = json.load(f)

__all__ = ['dbs_config', 'connect_to_mysql']