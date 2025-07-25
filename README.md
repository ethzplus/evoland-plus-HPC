[//]: # (Workflow: Quarto Pages Guide)
[//]: # (Workflow: Linting shell check)
[//]: # (R lifecycle batch)
[![Quarto: Pages Guide](https://github.com/ethzplus/Future-EI/actions/workflows/publish.yml/badge.svg?branch=main)](https://ethzplus.github.io/Future-EI/)
[![Shellcheck](https://github.com/ethzplus/Future-EI/actions/workflows/lint-sh.yml/badge.svg?branch=main)](https://github.com/ethzplus/Future-EI/actions/workflows/lint-sh.yml)
[![Lintr](https://github.com/ethzplus/Future-EI/actions/workflows/lint-R.yml/badge.svg?branch=main)](https://github.com/ethzplus/Future-EI/actions/workflows/lint-R.yml)
[![Hadolint](https://github.com/ethzplus/Future-EI/actions/workflows/lint-docker.yml/badge.svg?branch=main)](https://github.com/ethzplus/Future-EI/actions/workflows/lint-docker.yml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)


# About Future-EI

To explore future ecosystem services and nature's contributions to people (NCP)
in the context of climate scenarios, Future-EI (future ecological infrastructure) 
couples the
[land use land cover change (LULCC) model for Switzerland](https://github.com/blenback/LULCC-CH),
state-of-the-art
[species distribution modelling (N-SDM)](https://github.com/N-SDM/N-SDM),
and NCP calculations.

The repository contains the following:

1. Guide for the preparation and use of Future-EI in the
   [`FutureEiGuide`](FutureEiGuide) folder.
2. SLURM batch jobs as `bash` scripts. Find code in the [`src`](src) folder.

It specifically does not include the data and models used to run the workflow,
but an example is provided in the guide.

## Setup

For the setup of Future-EI, please refer to
the [setup guide](https://ethzplus.github.io/Future-EI/setup).
Generally, a shared configuration file can be found at
[`src/config.yml`](src/config.yml).

### Dependencies

For every step in the Future-EI pipeline, dependencies might differ.
Each step has its environment activated and in each step ([`steps`](src/steps)).
...

## Acknowledgements

This project was funded under the [ValPar.CH](https://valpar.ch/index_en.php?page=home_en)
project, funded as part of a pilot project of the "Action Plan for the Swiss 
Biodiversity Strategy (AP SBS)" by the Federal Office for the Environment (FOEN).

[N-SDM](https://github.com/N-SDM/N-SDM) has been developed within the
[Ecospatial Ecology Group (Ecospat)](https://www.unil.ch/ecospat/en/home.html)
at the [University of Lausanne](https://www.unil.ch/central/en/home.html)

[LULCC-CH](https://github.com/blenback/LULCC-CH)...

## References

- Adde, Antoine, Pierre-Louis Rey, Philipp Brun, Nathan Külling, Fabian Fopp, 
  Florian Altermatt, Olivier Broennimann, et al. 2023. “N-SDM : A High-Performance 
  Computing Pipeline for Nested Species Distribution Modelling.” Ecography 2023 (6): 
  e06540. https://doi.org/10.1111/ecog.06540.
- Black, Benjamin, Maarten J. Van Strien, Antoine Adde, and
  Adrienne Grêt-Regamey. 2023. “Re-Considering the Status Quo: Improving Calibration 
  of Land Use Change Models Through Validation of Transition Potential Predictions.” 
  Environmental Modelling & Software 159 (January): 105574.
  https://doi.org/10.1016/j.envsoft.2022.105574.
- Mayer, Paula, Sven-Erik Rabe, and Adrienne Grêt-Regamey. 2023. “Operationalizing 
  the Nature Futures Framework for Ecological Infrastructure.” Sustainability Science,
  July. https://doi.org/10.1007/s11625-023-01380-7.
