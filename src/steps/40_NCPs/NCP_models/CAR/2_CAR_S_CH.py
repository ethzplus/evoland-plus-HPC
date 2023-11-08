# coding=UTF-8
# -----------------------------------------------
#
# Model: InVEST Carbon Model
# Python version : 3.8 
# -----------------------------------------------
# This script runs the Carbon model from InVEST. It first defines the folders, where the biophysical tables and the land use files are. Then it creates a list containing the names
# for each region and each biophiscal table and land use file associated to it. Then it runs the module for each region.

# import modules

import natcap.invest.carbon, os, string, yaml

# Requires to have all natcap.invest dependencies installed (e.g. GDAL, rtree, numpy, etc). Cannot run on python 2!!!

# Load the parameters from ../40_NCPs_params.yml (relative to this file)
with open(os.path.join(os.path.dirname(__file__), '..',
                       '40_NCPs_params.yml')) as stream:
    params = yaml.safe_load(stream)

# ------ Define paths to variables

# Path to a folder containing biophysical tables for each region+altitude
bp_tables_dir = params['CAR']['bp_tables_dir']
# Path to a folder containing LULC rasters for each region
lu_dir = params['CAR']['lu_dir']
# Path to output folder containing each InVEST model
output_dir = params['CAR']['output_dir']

# Process: browsing through folders and defining names

# Lists all files in the biophysical tables folder.
listdir = os.listdir(bp_tables_dir)
list_names = []

print("List of regions:")

# Extracts the region name from the filename by removing the last 4 characters (file extension) and adding them to list_names
for i in range(0, len(listdir)):
    name = listdir[i]
    name = name[:-4]
    list_names.append(name)
    print(name)

# Process: Individually running InVEST Carbon model on each raster map from "lu_dir"

# specifying filenames for getting biophysical table & land use and specifying the output file path
for i in range(0, len(list_names)):
    bptable = bp_tables_dir + "\\" + list_names[i] + ".csv"
    lu = lu_dir + "\\" + list_names[i] + ".tif"
    out = output_dir + "\\" + list_names[i]
    tag = list_names[i]

    args = {
        # Dictionary containing the input parameters required to run the InVEST Carbon Model for the current region.
        'calc_sequestration': False,
        'carbon_pools_path': bptable,
        'do_redd': False,
        'do_valuation': False,
        'lulc_cur_path': lu,
        'results_suffix': tag,
        'workspace_dir': out,
    }

    # Executing the model calculation for each region
    if __name__ == '__main__':
        natcap.invest.carbon.execute(args)
    print("InVEST model " + tag + " :done")

print(".........................................")
print("script 2, year 2018 done!")
