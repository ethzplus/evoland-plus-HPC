# Natures Contribution to People (NCPs)

This document describes the scripts used to calculate the NCPs indicators.
The categories of the NCPs are taken from the [IPBES](https://www.ipbes.net/).

- [CAR: Regulation of climate](#car-regulation-of-climate)
- [FF: Food and feed](#ff-food-and-feed)
- [HAB: Habitat creation and maintenance](#hab-habitat-creation-and-maintenance)
- [NDR: Nutrient Delivery Ratio](#ndr-nutrient-delivery-ratio)
- [POL: Pollination and dispersal of seeds](#pol-pollination-and-dispersal-of-seeds)
- [SDR: Formation, protection and decontamination of soils](#sdr-formation-protection-and-decontamination-of-soils)
- [WY: Regulation of freshwater quantity, location and timing](#wy-regulation-of-freshwater-quantity-location-and-timing)

## CAR: Regulation of climate

#### Indicator: Carbon stored in biomass and soil

[`1_CAR_S_CH.R`](NCP_models/CAR/1_CAR_S_CH.R):
The aim of this script is to reclassify the Land-use map for the carbon
storage NCP mapping, based on altitude (DEM) and production region of
Switzerland.
The outputs are one land-use map per region/elevation.

Elevation (DEM) and land-use/land-cover (LULC) data are loaded,
the DEM reclassified into altitude classes, and intersects them with
production regions.
The script then reclassifies the LULC data into specific categories and
creates LULC rasters for each region based on the intersection of elevation
and production regions.

[`2_CAR_S_CH.py`](NCP_models/CAR/2_CAR_S_CH.py):
This script runs the Carbon model from InVEST.
First, the input folders are defined, where the biophysical tables and the land
use files are.
Then, a list containing the names for each region and each biophysical
table and land use file associated to it are created.
For each region, the model is run with the corresponding biophysical table and
land use file.

[`3_CAR_S_CH.R`](NCP_models/CAR/3_CAR_S_CH.R):
The aim of this script is to reunite the different InVEST models for each region
previously calculated in the other code files into one final map.

## FF: Food and feed

#### Indicator: Crop production potential (`ecocrop`)

[`0_FF_ecocrop`](NCP_models/FF/data_preparation/0_FF_ecocrop.R):
The Ecocrop model is applied to each crop of interest to generate
suitability maps.
The script reads precipitation, temperature, pH and crop datasets and
preprocesses them.
For each crop in the list, the script runs the Ecocrop model, which predicts
the suitability of the crop based on the provided environmental variables.
This only needs to be done once for each climate scenario, it is not
dependent on any LULC data.

[`1_FF_aggregate_masking.R`](NCP_models/FF/1_FF_aggregate_masking.R):
The aim of this script is to first, aggregate map outputs from ecocrop
models, and second, mask aggregated map with only agricultural land.

The code processes a land use/land cover raster to identify agricultural
areas and then uses this information to mask a stack of Ecocrop map outputs.
The masked raster values are then normalized to a range between $0$ and $1$.

## HAB: Habitat creation and maintenance

#### Indicator: Habitat quality index

[`1_HAB_S_CH.py`](NCP_models/HAB/1_HAB_S_CH.py):
This script runs the Habitat Quality model from InVEST.
Additionally, the parameters and input files are defined before the model is
run.

[`0_threat_layers_generation.R`](NCP_models/HAB/0_threat_layers_generation.R):
Creating threat raster layers for Habitat quality module InVEST model is the
aim of this script.
These layers are based on LULC maps and municipality polygons.

The Swiss municipality layer is processed to calculate population density and
identify rural residential areas.
The `threat_hab` function is then defined to create the threat layers.
This function classifies the land use/land cover raster into crop, rural
residential, and urban threat layers.

## NDR: Nutrient Delivery Ratio

#### Indicator: Annual nutrient retention by vegetation 

[`1_NDR_S_CH.py`](NCP_models/NDR/1_NDR_S_CH.py):
This script runs the Nutrient Delivery Ratio model from InVEST.
It defines the needed arguments (inputs) by specifying their values or the
path to the input files.
Additionally, before it creates a handler that specifies which log output is
sent to the console and how it is formatted.

## POL: Pollination and dispersal of seeds

#### Indicator: Habitat abundance for pollinators

[`1_POL_S_CH.py`](NCP_models/POL/1_POL_S_CH.py):
This script runs the Crop Pollination model from InVEST.
The parameters and input files are created before the model is run.

[`2_POL_S_CH_aggregating.R`](NCP_models/POL/2_POL_S_CH_aggregating.R):
The aim of this script is to aggregate the pollination information from the
InVEST model.

The code stacks all files of the data_folder with `pattern = "supply"` and then
calculates the sum and mean of the raster stack, aggregating the pollination
information.
Finally, a normalization to values between $0$ and $1$ is conducted.

## SDR: Formation, protection and decontamination of soils

#### Indicator: Erosion control by sediment retention

[`1_SDR_S_CH.py`](NCP_models/SDR/1_SDR_S_CH.py):
This script runs the Sediment Delivery Ratio model from InVEST.
Then it defines the needed arguments (inputs) by specifying their values or the
path to the input files.
Finally, it runs the InVEST model.
Additionally, before it creates a handler that specifies which log output is
sent to the console and how it is formatted.

## WY: Regulation of freshwater quantity, location and timing

#### Indicator: Annual water yield

[`1_WY_S_CH.py`](NCP_models/WY/1_WY_S_CH.py):
This script specifies the arguments for the InVEST Hydropower Water Yield model
and then runs it.
