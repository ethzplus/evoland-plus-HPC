# Path: src/steps/40_NCPs/NCP_models/load_params.py
# Script to load the NCP parameters

# Load libraries
import logging
import sys
from os import cpu_count

from natcap.invest.utils import LOG_FMT
from ruamel.yaml import YAML
from os.path import join, dirname
from os import environ, makedirs


def _get_param(params, param):
    """
    Function to get a parameter from a dictionary using a list of keys.

    >>> _get_param({'a': {'b': 1}}, ['a'])
    {'b': 1}
    >>> _get_param({'a': {'b': 1}}, ['a', 'c'])
    KeyError
    >>> _get_param({'a': {'b': {'c': 1}}}, ['a', 'b', 'c'])
    1
    >>> _get_param({'a': 1}, [])
    Exception
    >>> _get_param({}, ['a'])
    KeyError

    :param params: dictionary to get the parameter from
    :type params: dict
    :param param: list of keys to use to get the parameter
    :type param: List[str]
    :return: parameter
    :rtype: any

    :exception: KeyError if the parameter does not exist
    :exception: Exception if the key list is empty
    """
    # Check that the parameter exists
    if len(param) == 0:
        raise Exception("Key list cannot be empty")
    if param[0] not in params:
        raise KeyError("Key %s not found in dictionary %s", param[0], params)
    # Get the parameter
    if len(param) == 1:
        # If this is the last key in the list, return the parameter
        return params[param[0]]
    else:
        # Recursively get the parameter
        return _get_param(params[param[0]], param[1:])


def _add_run_params(params: dict) -> None:
    """
    Function that adds the run parameters to the parameters dictionary in place.
    For each run there are three parameters that are dynamically obtained from the
    environment variables:
    - NCP_RUN_YEAR:        year of the run
    - NCP_RUN_LULC_MAP:    path to the land use/land cover map
    - NCP_RUN_OUTPUT_DIR:  path to the output directory
    - NCP_RUN_SCRATCH_DIR: path to the scratch directory

    These parameters are added to the subdictionary 'run_params'.

    :param params: parameters dictionary
    :return: None
    """
    # Get the run parameters from the environment variables
    params['run_params'] = {
        env_var: environ[env_var]
        for env_var in
        ['NCP_RUN_SCENARIO_ID', 'NCP_RUN_YEAR', 'NCP_RUN_INPUT_DIR',
         'NCP_RUN_OUTPUT_DIR', 'NCP_RUN_SCRATCH_DIR']
    }
    # Set the LULC_M_EI_LAYER_DIR
    params['run_params']['LULCC_M_EI_LAYER_DIR'] = join(
        environ['LULCC_CH_HPC_DIR'],
        environ['LULCC_M_EI_LAYER_DIR'],
        "Future_EI",
        f"EI_ID{params['run_params']['NCP_RUN_SCENARIO_ID']}",
        params['run_params']['NCP_RUN_YEAR']
    )
    ## Set lulc path
    # {NCP_RUN_INPUT_DIR}/{NCP_RUN_SCENARIO_ID}/simulated_LULC_simID_{NCP_RUN_SCENARIO_ID}_year_{NCP_RUN_YEAR}.tif
    params['data']['lulc'] = join(
        params['run_params']['NCP_RUN_INPUT_DIR'],
        params['run_params']['NCP_RUN_SCENARIO_ID'],
        f"simulated_LULC_simID_{params['run_params']['NCP_RUN_SCENARIO_ID']}_year_"
        f"{params['run_params']['NCP_RUN_YEAR']}.tif"
    )

    ## Set NCP_RUN_RCP
    # Get current RCP from $LULCC_M_SIM_CONTROL_TABLE
    from pandas import read_csv
    df = read_csv(environ['LULCC_M_SIM_CONTROL_TABLE'],
                  usecols=["Simulation_ID.string", "Climate_scenario.string",
                           "Scenario_end.real", "Step_length.real"])
    # get the Climate_scenario.string only where
    # NCP_RUN_SCENARIO_ID==Simulation_ID.string, should be unique
    df = df[df['Simulation_ID.string'] == int(
        params['run_params']['NCP_RUN_SCENARIO_ID'])]
    if df.empty:
        raise Exception(
            "Simulation_ID.string not found in LULCC_M_SIM_CONTROL_TABLE")
    if len(df) > 1:
        raise Exception(
            "Simulation_ID.string not unique in LULCC_M_SIM_CONTROL_TABLE")
    params['run_params']['NCP_RUN_RCP'] = df['Climate_scenario.string'].values[
        0].upper()
    # read Scenario_end.real,Step_length.real for use later
    end_step = int(df['Scenario_end.real'].values[0])
    step_length = int(df['Step_length.real'].values[0])

    # Number workers
    # slurm variables should take precedence
    if params['other']['n_workers'] != -1:
        if 'SLURM_CPUS_PER_TASK' in environ:
            params['other']['n_workers'] = int(environ['SLURM_CPUS_PER_TASK'])
        else:
            params['other']['n_workers'] = cpu_count()

    residential_year = (
        params['run_params']['NCP_RUN_YEAR']
        if int(params['run_params']['NCP_RUN_YEAR']) < end_step
        else end_step - step_length  # for the last, use second to last step
    )
    ## Rural residential path
    params['data']['rur_res'] = join(
        params['run_params']['NCP_RUN_INPUT_DIR'],
        params['run_params']['NCP_RUN_SCENARIO_ID'],
        f"rur_res_simID_{params['run_params']['NCP_RUN_SCENARIO_ID']}_year_"
        f"{residential_year}.tif"
    )
    ## Urban residential path
    params['data']['urb_res'] = join(
        params['run_params']['NCP_RUN_INPUT_DIR'],
        params['run_params']['NCP_RUN_SCENARIO_ID'],
        f"urb_res_simID_{params['run_params']['NCP_RUN_SCENARIO_ID']}_year_"
        f"{residential_year}.tif"
    )

    # Erosivity - based on current RCP, always 2060 (Erosivity_2060_RCP26.tif)
    # 2060_RCP26 and 2060_RCP45 (latter also for RCP85)
    # TODO: change for different time periods with new data if needed
    erosivity_rcp = "26" if params['run_params'][
                                'NCP_RUN_RCP'] == "RCP26" else "45"
    params['data']['erosivity_path'] = join(
        params['data']['erosivity_path'],
        f"Erosivity_2060_RCP{erosivity_rcp}.tif"
    )

    # Evapotranspiration "etp_{year}_{rcp}.tif"
    params['data']['eto'] = join(
        params['data']['eto'],
        f"etp_{params['run_params']['NCP_RUN_YEAR']}_"
        f"{params['run_params']['NCP_RUN_RCP']}.tif"
    )

    # Yearly precipitation "Prec_{year}_{rcp}.tif"
    params['data']['yearly_precipitation'] = join(
        params['data']['yearly_precipitation'],
        f"Prec_{params['run_params']['NCP_RUN_YEAR']}_"
        f"{params['run_params']['NCP_RUN_RCP']}.tif"
    )

    # Monthly precipitation


def setup_params(cache_file: str, check_params=None) -> None:
    """
    Function to load and cache the parameters and save them in the cache
    directory.
    It will load the parameters from environment variable NCP_PARAMS_YML if set,
    otherwise from ./40_NCPs_params.yml,
    and save them to the cache file.

    NCPS_PARAMS_YML is set to the cache directory,
    so further calls to load_params will use the cached parameters.

    :param cache_file: cache file path
    :type cache_file: str
    :param check_params: list of key lists to check for existence
                         lists can be of any length > 0
    :type check_params: List[List[str]]
    :return: None
    :rtype: None
    """
    params = load_params(check_params=check_params, add_run_params=True)

    # Create the output directory if it does not exist
    # makedirs(params['run_params']['NCP_RUN_OUTPUT_DIR'], exist_ok=True)
    makedirs(dirname(cache_file), exist_ok=True)
    # Cache the parameters
    with open(cache_file, 'w') as stream:
        YAML().dump(params, stream)


def load_params(check_params=None, add_run_params=False):
    """
    Function to load the parameters from environment variable NCP_PARAMS_YML if
    set, otherwise from ./40_NCPs_params.yml, and return them as a dictionary.

    >>> params = load_params(check_params=[
    >>>     ['CAR', 'bp_tables_dir'],
    >>>     ['CAR', 'water', 'depth'],
    >>> ])

    :param check_params: list of key lists to check for existence
                         lists can be of any length > 0
    :type check_params: List[List[str]]
    :param add_run_params: whether to add the run parameters to the parameters
    :type add_run_params: bool
    :return: parameters
    :rtype: dict
    """
    # Try to load the parameters from the environment
    if check_params is None:
        check_params = list()

    yaml = YAML()

    # Load the parameters from env var NCP_PARAMS_YML if set, otherwise from
    # ./40_NCPs_params.yml
    if environ['NCP_PARAMS_YML'] == "":
        print("NCP_PARAMS_YML environment variable is not set. "
              "Loading from ./40_NCPs_params.yml")
        with open(join(dirname(__file__), '40_NCPs_params.yml')) as stream:
            params = yaml.load(stream)
    else:
        print(f"Loading NCP parameters from {environ['NCP_PARAMS_YML']}")
        with open(environ['NCP_PARAMS_YML']) as stream:
            params = yaml.load(stream)

    # Add the run parameters to the parameter dictionary
    if add_run_params:
        _add_run_params(params)

    # Check that all required parameters are present
    for param in check_params:
        # Check that the parameter exists
        try:
            _get_param(params, param)
        except KeyError:
            raise Exception("Parameter %s not found in parameter file.", param)

    return params


# Logging for natcap.invest models

# Creates a logger object for the current module, which can be used to log
# messages specific to this module.
LOGGER = logging.getLogger(__name__)

# Creates a root logger object, which can be used to log messages at the
# application level.
root_logger = logging.getLogger()

# Creates a logging handler that sends log messages to the standard output (
# console).
handler = logging.StreamHandler(sys.stdout)

# Defines the format for logging messages using the format provided by the
# natcap.invest toolkit. This format includes details like timestamps,
# log levels, and the actual message content.
formatter = logging.Formatter(
    fmt=LOG_FMT,
    datefmt='%m/%d/%Y %H:%M:%S ')

# Sets the formatter for the logging handler.
handler.setFormatter(formatter)

# Configures the basic settings for logging. It sets the log level to INFO,
# meaning only messages of level INFO and above will be logged. It also adds
# the handler
logging.basicConfig(level=logging.INFO, handlers=[handler])

# If this script is run as the main program, it will cache the parameters
# and save them in the cache directory.
# The Cache file name needs to be passed as an argument
if __name__ == '__main__':
    setup_params(sys.argv[1])
    logging.info(f"Cached parameters to {sys.argv[1]}")
