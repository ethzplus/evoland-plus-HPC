# coding=UTF-8
# -----------------------------------------------
#
# Model: InVEST Carbon Model
# Python version : 3.8 
# -----------------------------------------------
# This script runs the Carbon model from InVEST. It first defines the folders,
# where the biophysical tables and the land use files are.
# Then it creates a list containing the names for each region and each
# biophiscal table and land use file associated to it.
# Then it runs the module for each region.
import sys
from os import listdir, makedirs
from os.path import join, basename, dirname

import natcap.invest.carbon

sys.path.append(join(dirname(__file__), '..'))  # bad practice
from load_params import load_params

# Requires to have all natcap.invest dependencies installed
# (e.g. GDAL, rtree, numpy, etc). Cannot run on python 2!

params = load_params(
    check_params=[
        ['CAR', 'bp_tables_dir'],
        ['run_params', 'NCP_RUN_YEAR'],
        ['run_params', 'NCP_RUN_OUTPUT_DIR'],
        ['run_params', 'NCP_RUN_SCRATCH_DIR'],
    ]
)

# ------ Define paths to variables

# Path to a folder containing biophysical tables for each region+altitude
bp_tables_dir = params['CAR']['bp_tables_dir']
# Path to a folder containing the prepared LULC rasters for each region
lulc_clip_dir = join(params['run_params']['NCP_RUN_SCRATCH_DIR'], 'CAR',
                     'lulc_clip', params['run_params']['NCP_RUN_YEAR'])
# Path to output folder containing each InVEST model
model_dir = join(params['run_params']['NCP_RUN_OUTPUT_DIR'], 'CAR',
                 params['run_params']['NCP_RUN_YEAR'], 'Invest_models')


def _get_region_names(dir_path: str) -> list[str]:
    """
    This function returns of all files in a folder without the file extension.
    :param dir_path: Path to the folder containing the files
    :return: List of file names without the file extension
    """
    return [basename(f).split('.')[0] for f in listdir(dir_path)]


if __name__ == '__main__':

    list_names = _get_region_names(bp_tables_dir)
    print(f"CAR: Starting Carbon model for {len(list_names)} regions...\n"
          f"     Regions: {list_names}")

    # create output folder
    makedirs(model_dir, exist_ok=True)

    for region in list_names:
        args = {
            'calc_sequestration': False,
            'carbon_pools_path': join(bp_tables_dir, f"{region}.csv"),
            'do_redd': False,
            'do_valuation': False,
            # 'lulc_cur_path': lu_map,  # join(lu_dir, f"{region}.tif"),
            'lulc_cur_path': join(lulc_clip_dir, f"{region}.tif"),
            'results_suffix': region,
            'workspace_dir': join(model_dir, region),
        }
        natcap.invest.carbon.execute(args)

    print("CAR: ...done!")
