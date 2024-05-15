# Path: src/steps/40_NCPs/NCP_models/load_params.py
# Script to load the NCP parameters

# Load libraries
from yaml import safe_load
from os.path import join, dirname
from os import environ


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
        ['NCP_RUN_YEAR', 'NCP_RUN_LULC_MAP', 'NCP_RUN_OUTPUT_DIR',
         'NCP_RUN_SCRATCH_DIR']
    }


def load_params(check_params=None):
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
    :return: parameters
    :rtype: dict
    """
    # Try to load the parameters from the environment
    if check_params is None:
        check_params = list()

    # Load the parameters from env var NCP_PARAMS_YML if set, otherwise from
    # ./40_NCPs_params.yml
    if environ['NCP_PARAMS_YML'] == "":
        print("NCP_PARAMS_YML environment variable is not set. "
              "Loading from ./40_NCPs_params.yml")
        with open(join(dirname(__file__), '40_NCPs_params.yml')) as stream:
            params = safe_load(stream)
    else:
        print("Loading NCP parameters from %s", environ['NCP_PARAMS_YML'])
        with open(environ['NCP_PARAMS_YML']) as stream:
            params = safe_load(stream)

    # Add the run parameters to the parameters dictionary
    _add_run_params(params)

    # Check that all required parameters are present
    for param in check_params:
        # Check that the parameter exists
        try:
            _get_param(params, param)
        except KeyError:
            raise Exception("Parameter %s not found in parameter file.", param)

    return params
