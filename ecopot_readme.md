This README.txt file was generated on 2026-01-21 by Benjamin Black

GENERAL INFORMATION

1.  Title of Dataset: Bern Landscape Change Explorer Ecosystem Condition
    and Potential Results.

2.  Author Information A. Principal Investigator Contact Information
    Name: Isabel Nicholson Thomas Institution: Institute for Spatial and
    Landscape Development, Swiss Federal Institute of Technology (ETH)
    Address: HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich,
    Switzerland Email: inthomas@ethz.ch

    B. Additional researcher Contact Information Name: Benjamin Black
    Institution: Institute for Spatial and Landscape Development, Swiss
    Federal Institute of Technology (ETH) Address: HIL H 52.1,
    Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland Email:
    bblack@ethz.ch

C. Additional researcher Contact Information Name: Manuel Kurmann
Institution: Institute for Spatial and Landscape Development, Swiss
Federal Institute of Technology (ETH) Address: HIL H 52.1,
Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland Email:
mankurma@ethz.ch

3.  Date of data collection: 2025-06-01 \> 2026-01-1

4.  Geographic location of data: Switzerland

SHARING/ACCESS INFORMATION

1.  Licenses/restrictions placed on the data: Creative Commons
    Attribution 4.0 International

2.  Recommended citation for this dataset: Nicholson Thomas, I., Black,
    B., & Kurmann, M. (2026). Bern Landscape Change Explorer Ecosystem
    Condition and Potential Results (Version 1.0.0) \[Data set\]. ETH
    Zurich. https://doi.org/10.5281/zenodo.18327368

DESCRIPTION/REPRODUCIBILITX: - This dataset contains rasters of
ecosystem condition and potential for the Bern Landscape Change Explorer
project. - All rasters are projected in the CH1903+ / LV95 coordinate
reference system and share the same extent and resolution (1 hectare). -
For detailed descriptions of the data generation process, raster value
interpretation and suggested colour scheme, please refer to the
accompanying metadata files: -
BeschreibungMetadatenOekosystemzustand_Potenzial_DE.pdf -
BeschreibungMetadatenOekosystemzustand_Potenzial_ENG.pdf

This dataset contains raster data on ecosystem condition (EC) and the
potential for condition dependent improvements in selected ecosystem
services. The calculations were performed separately for three ecosystem
types: forests, grasslands and alpine pastures, and agricultural land.
The data sets on ecosystem condition are stored as ec_index and
separated by ecosystem type. ec_index_LU_Forest contains the ecosystem
condition index for forest areas, ec_index_LU_GrassPasture for
grasslands and alpine pastures, and ec_index_LU_IntensiveAg for
intensive agricultural areas. In all three cases, the pixel value
corresponds to the calculated, composite ecosystem condition index for
the respective land use group. In addition, there is the cluster_EC_all
dataset, which represents spatial patterns in ecosystem condition as
classes. The clusters are based on the same indicator layers that were
used to calculate the ecosystem condition index. Each grid cell is
assigned a cluster ID, which groups together cells that have similar
characteristics across several indicators, for example, simultaneously
high or low values in the same state characteristics. This allows areas
to be delineated that are similar in terms of the interaction of the
state indicators considered. The potential datasets are named after
ecosystem services and end in \_gap. CAR_gap describes the potential for
carbon storage, FF_gap for food and feed production, HAB_gap for habitat
quality, POL_gap for pollination, and REC_gap for recreational
potential. The pixel values indicate the difference between theoretical
and effective provision, i.e., the additional amount that could be
achieved with an improved ecosystem status.

Included files: BeschreibungMetadatenOekosystemzustand_Potenzial_DE.pdf
BeschreibungMetadatenOekosystemzustand_Potenzial_ENG.pdf CAR_gap.tif
ec_index_combined.tif ec_index_LU_13.tif ec_index_LU_15.tif
ec_index_LU_17.tif FF_gap.tif HAB_gap.tif POL_gap.tif REC_gap.tif
