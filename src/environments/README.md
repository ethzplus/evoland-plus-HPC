# Environments

This folder collects the environments, their initial setup and frozen environments.

## Conda

For conda envirionments, the `{xx}_{env_name}.sh` scripts are initial installation 
instructions and can usually be executed to set up the whole environment. For 
clarity read the comments in the scripts. To make the scripts executable, run
`chmod +x {xx}_{env_name}.sh`, or for all in the current folder `chmod +x *.sh`.
If needed there is a `{xx}_{env_name}_requirements.txt` file that is used to
resolve the dependencies, referred in the setup scripts.
Finally there is a `{xx}_{env_name}_env.yml` file that is the frozen environment
with the exact build versions and can be used to recreate the exact environment with
`conda env create -f {xx}_{env_name}_env.yml`.