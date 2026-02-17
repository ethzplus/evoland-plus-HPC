#' Data archiving
#'
#'
#'
#'
library(fs)
library(rmarkdown)

raster_dir <- "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform/raster_data"
metadata_dir <- "X:/CH_Kanton_Bern/03_Workspaces/08_Metadata/Metadaten"


### =========================================================================
### LULC data archiving
### =========================================================================

# copy all LULC rasters for final data archiving
lulc_archive_dir <- file.path(
  "X:/CH_Kanton_Bern/03_Workspaces/09_Data_archiving/blce-lulc-data-archive"
)
if (!dir.exists(lulc_archive_dir)) {
  dir.create(lulc_archive_dir, recursive = TRUE, showWarnings = FALSE)
}

all_lulc_rasters <- list.files(
  raster_dir,
  pattern = "^lulc-.*//.tif$",
  full.names = TRUE
)

# copy files to archive directory
file.copy(
  all_lulc_rasters,
  lulc_archive_dir,
  overwrite = TRUE
)

# Move the corresponding metadata file to the archive directory
# match on string LULC with extension pdf
lulc_metadata_files <- list.files(
  metadata_dir,
  pattern = "LULC.*//.pdf$",
  full.names = TRUE
)

file.copy(
  lulc_metadata_files,
  lulc_archive_dir,
  overwrite = TRUE
)

# produce a file tree of the archived data including only file names not directories
file_tree <- dir_tree(
  lulc_archive_dir,
  recurse = TRUE,
  type = "file"
)

# convert to basenames only
file_tree <- basename(file_tree)


# add this text to the top of the file tree
lulc_metadata_content <- c(
  "This README.txt file was generated on 2026-01-21 by Benjamin Black

GENERAL INFORMATION

1. Title of Dataset: Bern Landscape Change Explorer Land Use and Land Cover Change Simulations.    

2. Author Information
	A. Principal Investigator Contact Information
		Name: Benjamin Black
		Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH) 
		Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: bblack@ethz.ch

	B. Additional researcher Contact Information
		Name: Manuel Kurmann
		Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH) 
		Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: mankurma@ethz.ch

  C. Additional researcher Contact Information
		Name: Isabel Nicholson Thomas
    Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH)
    Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: inthomas@ethz.ch


3. Date of data collection: 2025-06-01 > 2026-01-1 

4. Geographic location of data: Switzerland 


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: Creative Commons Attribution 4.0 International 

2. Recommended citation for this dataset: Black, B., Kurmann, M., & Nicholson Thomas, I. (2026). Bern Landscape Change Explorer Land Use and Land Cover Change Scenarios (Version 1.0.0) [Data set]. ETH Zurich. https://doi.org/10.5281/zenodo.18324442

DESCRIPTION/REPRODUCIBILITX: 
- This dataset contains rasters of land use and land cover (LULC) change scenarios for the Bern Landscape Change Explorer project.
- All rasters are projected in the CH1903+ / LV95 coordinate reference system and share the same extent and resolution (1 hectare).
- For detailed descriptions of the data generation process, raster value interpretation and suggested colour scheme, please refer to the accompanying metadata files:",
  paste("- ", basename(lulc_metadata_files)),
  "",
  "The simulations were carried out in five-year increments
for the period 2025 to 2060. Three climate scenarios, five economic scenarios, three population
scenarios, and five land use scenarios were taken into account. The land use scenarios control the
future transition probabilities between land use classes. A total of 2,025 LULC files are available for
all combinations of year, climate, economy, population, and land use scenario.
The file names follow the pattern:
lulc-<year>-<climate scenario>-<economic scenario>-<population scenario>-<land use scenario>-
canton.tif
Example: lulc-2025-rcp45-ref_peri_urban-low-ei_cul-canton.tif",
  "",
  "Category                             Options
---------------------------------------------------------------------------------------------------------------
Year                                   2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060
Climate scenarios                      rcp26, rcp45, rcp85
Economic scenarios                     combined_urban, ecolo_urban, ecolo_central, ref_central, ref_peri_urban
Demographic scenarios                  low, ref, high
Land use change rates scenarios        bau, ei_cul, ei_nat, ei_soc, gr_ex",
  "",
  "Included files:",
  file_tree
)

lulc_readme_path <- file.path(lulc_archive_dir, "README.txt")
writeLines(lulc_metadata_content, con = lulc_readme_path)

# path to markdown readme file
lulc_md_path <- file.path(getwd(), "lulc_readme.md")

# Convert markdown to plain text
pandoc_convert(
  input = lulc_readme_path,
  from = "plain",
  to = "markdown",
  output = lulc_md_path,
  verbose = TRUE
)

system2(
  command = "pandoc",
  args = c(
    lulc_readme_path,
    "-f",
    "plain",
    "-t",
    "markdown",
    "-o",
    lulc_md_path
  )
)


### =========================================================================
### ES data archiving
### =========================================================================

# copy all ES rasters for final data archiving
es_archive_dir <- file.path(
  "X:/CH_Kanton_Bern/03_Workspaces/09_Data_archiving/blce-es-data-archive"
)
if (!dir.exists(es_archive_dir)) {
  dir.create(es_archive_dir, recursive = TRUE, showWarnings = FALSE)
}

es_abbreviations <- c('rec', 'car', 'ndr', 'pol', 'sdr', 'wy', 'hab', 'ff')

all_es_rasters <- c(sapply(es_abbreviations, function(abbr) {
  list.files(
    raster_dir,
    pattern = paste0("^", abbr, ".*//.tif$"),
    full.names = TRUE
  )
}))


# copy files to archive directory
file.copy(
  all_es_rasters,
  es_archive_dir,
  overwrite = TRUE
)

# Move the corresponding metadata file to the archive directory
es_metadata_files <- list.files(
  metadata_dir,
  pattern = "Oekosystemleistungen.*//.pdf$",
  full.names = TRUE
)

file.copy(
  es_metadata_files,
  es_archive_dir,
  overwrite = TRUE
)

# produce a file tree of the archived data including only file names not directories
file_tree_es <- dir_tree(
  es_archive_dir,
  recurse = TRUE,
  type = "file"
)

# convert to basenames only
file_tree_es <- basename(file_tree_es)

# add this text to the top of the file tree
es_metadata_content <- c(
  "This README.txt file was generated on 2026-01-21 by Benjamin Black

GENERAL INFORMATION

1. Title of Dataset: Bern Landscape Change Explorer Ecosystem Service Provision Simulations.    

2. Author Information
	A. Principal Investigator Contact Information
		Name: Benjamin Black
		Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH) 
		Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: bblack@ethz.ch

	B. Additional researcher Contact Information
		Name: Manuel Kurmann
		Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH) 
		Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: mankurma@ethz.ch

  C. Additional researcher Contact Information
		Name: Isabel Nicholson Thomas
    Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH)
    Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: inthomas@ethz.ch


3. Date of data collection: 2025-06-01 > 2026-01-1 

4. Geographic location of data: Switzerland 


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: Creative Commons Attribution 4.0 International 

2. Recommended citation for this dataset: Black, B., Kurmann, M., & Nicholson Thomas, I. (2026). Bern Landscape Change Explorer Ecosystem Service Provision Simulations (Version 1.0.0) [Data set]. ETH Zurich. https://doi.org/10.5281/zenodo.18325820

DESCRIPTION/REPRODUCIBILITX: 
- This dataset contains rasters of Ecosystem Service Provision scenarios for the Bern Landscape Change Explorer project.
- All rasters are projected in the CH1903+ / LV95 coordinate reference system and share the same extent and resolution (1 hectare).
- For detailed descriptions of the data generation process, raster value interpretation and suggested colour scheme, please refer to the accompanying metadata files:",
  paste("- ", basename(es_metadata_files)),
  "",
  "This dataset contains files on eight ecosystem services. The simulations were carried out in five-year
increments for the period 2025 to 2060. Three climate scenarios, five economic scenarios, three
population scenarios, and five land use scenarios were taken into account. The latter control the
future transition probabilities between different land use classes. All possible combinations were
calculated for each ecosystem service, resulting in a total of 2,025 files per ecosystem service.
The file names follow the pattern:
< Ecosystem service>-<Year>-<Climate scenario>-<Economic scenario>-<Population scenario>-
<land use scenario>-canton.tif
Example: car-2025-rcp45-ref_peri_urban-low-ei-cul-canton.tif",
  "",
  "Category                             Options
---------------------------------------------------------------------------------------------------------------
Ecosystem Service                      car, ff, hab, ndr, pol, rec, sdr, wy
Year                                   2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060
Climate scenarios                      rcp26, rcp45, rcp85
Economic scenarios                     combined_urban, ecolo_urban, ecolo_central, ref_central, ref_peri_urban
Demographic scenarios                  low, ref, high
Land use change rates scenarios        bau, ei_cul, ei_nat, ei_soc, gr_ex",
  "",
  "Included files:",
  file_tree_es
)

es_readme_path <- file.path(es_archive_dir, "README.txt")
writeLines(es_metadata_content, con = es_readme_path)

# path to markdown readme file
es_md_path <- file.path(getwd(), "es_readme.md")
# Convert markdown to plain text
pandoc_convert(
  input = es_readme_path,
  to = "markdown",
  output = es_md_path
)

### =========================================================================
### SDM data archiving
### =========================================================================

sdm_dir <- "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform/sdm_data/raster_data"
sdm_archive_dir <- file.path(
  "X:/CH_Kanton_Bern/03_Workspaces/09_Data_archiving/blce-sdm-data-archive"
)
if (!dir.exists(sdm_archive_dir)) {
  dir.create(sdm_archive_dir, recursive = TRUE, showWarnings = FALSE)
}

all_sdm_rasters <- list.files(
  sdm_dir,
  pattern = ".*//.tif$",
  full.names = TRUE
)

# copy files to archive directory
file.copy(
  all_sdm_rasters,
  sdm_archive_dir,
  overwrite = TRUE
)

# Move the corresponding metadata file to the archive directory
sdm_metadata_files <- list.files(
  metadata_dir,
  pattern = "Habitateignung.*//.pdf$",
  full.names = TRUE
)
file.copy(
  sdm_metadata_files,
  sdm_archive_dir,
  overwrite = TRUE
)

# produce a file tree of the archived data including only file names not directories
file_tree_sdm <- dir_tree(
  sdm_archive_dir,
  recurse = TRUE,
  type = "file"
)

# convert to basenames only
file_tree_sdm <- basename(file_tree_sdm)

# add this text to the top of the file tree
sdm_metadata_content <- c(
  "This README.txt file was generated on 2026-01-21 by Benjamin Black

GENERAL INFORMATION

1. Title of Dataset: Bern Landscape Change Explorer Species Habitat Suitability Simulations.    

2. Author Information
	A. Principal Investigator Contact Information
		Name: Manuel Kurmann
		Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH) 
		Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: mankurma@ethz.ch

	B. Additional researcher Contact Information
    Name: Benjamin Black
		Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH) 
		Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: bblack@ethz.ch

  C. Additional researcher Contact Information
		Name: Isabel Nicholson Thomas
    Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH)
    Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: inthomas@ethz.ch


3. Date of data collection: 2025-06-01 > 2026-01-1 

4. Geographic location of data: Switzerland 


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: Creative Commons Attribution 4.0 International 

2. Recommended citation for this dataset: Kurmann, M., Black, B., & Nicholson Thomas, I. (2026). Bern Landscape Change Explorer Species Habitat Suitability Simulations (Version 1.0.0) [Data set]. ETH Zurich. https://doi.org/10.5281/zenodo.18326626

DESCRIPTION/REPRODUCIBILITX: 
- This dataset contains rasters of Species habitat suitability scenarios for the Bern Landscape Change Explorer project.
- All rasters are projected in the CH1903+ / LV95 coordinate reference system and share the same extent and resolution (1 hectare).
- For detailed descriptions of the data generation process, raster value interpretation and suggested colour scheme, please refer to the accompanying metadata files:",
  paste("- ", basename(sdm_metadata_files)),
  "",
  "This dataset contains raster files on habitat suitability. The projections are available in five-year
increments for the period from 2025 to 2060. Each projection combines three climate scenarios and
five land use scenarios. The climate scenarios control future climatic conditions, while the land use
scenarios control the spatial distribution of land use classes in the future. Land use is then applied
as a filter by setting habitat suitability in raster cells to zero if the land use there is considered
unsuitable for the species in question. A total of 6600 files are available for all combinations of year,
climate scenario, land use scenario, species list, and taxonomic group.
The file names follow the pattern:
sdm-<year>-<climate scenario>-<land use scenario>-<species list>-<group>-canton.tif
Example: sdm-2025-ssp126-ei_cul-uzl-amphibians-canton.tif",
  "",
  "Category                             Options
---------------------------------------------------------------------------------------------------------------
Year                                   2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060
Climate scenarios                      rcp26, rcp45, rcp85
Economic scenarios                     combined_urban, ecolo_urban, ecolo_central, ref_central, ref_peri_urban
Demographic scenarios                  low, ref, high
Land use change rates scenarios        bau, ei_cul, ei_nat, ei_soc, gr_ex
Species list                           uzl, npa, uzl_npa
Taxonomic group                        all, amphibians, beetles, birds, butterflies, caddisflies,
                                       dragonflies, fungi, grasshoppers, lichens, mayflies, mammals,
                                       mosses, neuropterans, reptiles, snails_bivalves, stoneflies,
                                       stoneworts, vascular plants, wild bees",
  "",
  "Included files:",
  file_tree_sdm
)

sdm_readme_path <- file.path(sdm_archive_dir, "README.txt")
writeLines(sdm_metadata_content, con = sdm_readme_path)

# path to markdown readme file
sdm_md_path <- file.path(getwd(), "sdm_readme.md")
# Convert markdown to plain text
pandoc_convert(
  input = sdm_readme_path,
  to = "markdown",
  output = sdm_md_path
)

### =========================================================================
### Robustness data archiving
### =========================================================================

robustness_dir <- "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform/robustness"
robustness_archive_dir <- file.path(
  "X:/CH_Kanton_Bern/03_Workspaces/09_Data_archiving/blce-robustness-data-archive"
)
if (!dir.exists(robustness_archive_dir)) {
  dir.create(robustness_archive_dir, recursive = TRUE, showWarnings = FALSE)
}

all_robustness_files <- list.files(
  robustness_dir,
  pattern = ".*//.tif$",
  full.names = TRUE
)

# remove any files containing string 'rescaled' in the name
all_robustness_files <- all_robustness_files[
  !grepl("rescaled", all_robustness_files)
]

# copy files to archive directory
file.copy(
  all_robustness_files,
  robustness_archive_dir,
  overwrite = TRUE
)

# Move the corresponding metadata file to the archive directory
robustness_metadata_files <- list.files(
  metadata_dir,
  pattern = "Robustheit.*//.pdf$",
  full.names = TRUE
)

file.copy(
  robustness_metadata_files,
  robustness_archive_dir,
  overwrite = TRUE
)

# produce a file tree of the archived data including only file names not directories
file_tree_robustness <- dir_tree(
  robustness_archive_dir,
  recurse = TRUE,
  type = "file"
)

# convert to basenames only
file_tree_robustness <- basename(file_tree_robustness)

robustness_metadata_content <- c(
  "This README.txt file was generated on 2026-01-21 by Benjamin Black

GENERAL INFORMATION

1. Title of Dataset: Bern Landscape Change Explorer Landscape Robustness Results.    

2. Author Information
	A. Principal Investigator Contact Information
		Name: Benjamin Black
		Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH) 
		Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: bblack@ethz.ch

	B. Additional researcher Contact Information
		Name: Manuel Kurmann
		Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH) 
		Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: mankurma@ethz.ch

  C. Additional researcher Contact Information
		Name: Isabel Nicholson Thomas
    Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH)
    Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: inthomas@ethz.ch


3. Date of data collection: 2025-06-01 > 2026-01-1 

4. Geographic location of data: Switzerland 


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: Creative Commons Attribution 4.0 International 

2. Recommended citation for this dataset: Black, B., Kurmann, M., & Nicholson Thomas, I. (2026). Bern Landscape Change Explorer Landscape Robustness Results (Version 1.0.0) [Data set]. ETH Zurich. https://doi.org/10.5281/zenodo.18327228

DESCRIPTION/REPRODUCIBILITX: 
- This dataset contains rasters of landscape robustness for the Bern Landscape Change Explorer project.
- All rasters are projected in the CH1903+ / LV95 coordinate reference system and share the same extent and resolution (1 hectare).
- For detailed descriptions of the data generation process, raster value interpretation and suggested colour scheme, please refer to the accompanying metadata files:",
  paste("- ", basename(robustness_metadata_files)),
  "",
  "Two raster-based indicators are provided for each topic for the robustness evaluation: performance and
stability. The performance indicator is stored in files with the prefix Mean_sum_of_change_ and
describes the mean cumulative change value across all scenarios and future points in time. The
stability indicator is stored in files with the prefix Undesirable_deviation_sum_of_change_ and
describes the extent of undesirable negative deviations across all scenarios. Both indicators are
available for biodiversity (BD), ecosystem services (ES), and the combination of biodiversity and
ecosystem services (BD-ES). The robustness classes on the platform are derived from the joint
classification of these two indicator grids and are therefore not included as a separate grid layer.",
  "",
  "Included files:",
  file_tree_robustness
)

robustness_readme_path <- file.path(robustness_archive_dir, "README.txt")
writeLines(robustness_metadata_content, con = robustness_readme_path)
# path to markdown readme file
robustness_md_path <- file.path(getwd(), "robustness_readme.md")
# Convert markdown to plain text
pandoc_convert(
  input = robustness_readme_path,
  to = "markdown",
  output = robustness_md_path
)

### =========================================================================
### Ecosystem condition and potential data archiving
### =========================================================================

ecopot_dir <- c(
  "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform/restoration_potential/December/Updated_ES_Gap",
  "X:/CH_Kanton_Bern/03_Workspaces/05_Web_platform/restoration_potential/December/Updated_EC"
)
all_ecopot_files <- c()
for (dir in ecopot_dir) {
  files <- list.files(
    dir,
    full.names = TRUE,
    pattern = ".tif",
    recursive = TRUE
  )
  all_ecopot_files <- c(all_ecopot_files, files)
}


ecopot_archive_dir <- file.path(
  "X:/CH_Kanton_Bern/03_Workspaces/09_Data_archiving/blce-oekosystemzustand-potenzial-data-archive"
)
if (!dir.exists(ecopot_archive_dir)) {
  dir.create(ecopot_archive_dir, recursive = TRUE, showWarnings = FALSE)
}

# copy files to archive directory
file.copy(
  all_ecopot_files,
  ecopot_archive_dir,
  overwrite = TRUE
)

# Move the corresponding metadata file to the archive directory
ecopot_metadata_files <- list.files(
  metadata_dir,
  pattern = "Oekosystemzustand_Potenzial.*//.pdf$",
  full.names = TRUE
)

file.copy(
  ecopot_metadata_files,
  ecopot_archive_dir,
  overwrite = TRUE
)

# produce a file tree of the archived data including only file names not directories
file_tree_ecopot <- dir_tree(
  ecopot_archive_dir,
  recurse = TRUE,
  type = "file"
)
# convert to basenames only
file_tree_ecopot <- basename(file_tree_ecopot)

ecopot_metadata_content <- c(
  "This README.txt file was generated on 2026-01-21 by Benjamin Black

GENERAL INFORMATION

1. Title of Dataset: Bern Landscape Change Explorer Ecosystem Condition and Potential Results.    

2. Author Information
	A. Principal Investigator Contact Information
		Name: Isabel Nicholson Thomas
    Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH)
    Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: inthomas@ethz.ch

	B. Additional researcher Contact Information
		Name: Benjamin Black
		Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH) 
		Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: bblack@ethz.ch

  C. Additional researcher Contact Information
		Name: Manuel Kurmann
		Institution:  Institute for Spatial and Landscape Development, Swiss Federal Institute of Technology (ETH) 
		Address:  HIL H 52.1, Stefano-Franscini-Platz 5 CH-8093 Zürich, Switzerland
		Email: mankurma@ethz.ch

3. Date of data collection: 2025-06-01 > 2026-01-1 

4. Geographic location of data: Switzerland 

SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: Creative Commons Attribution 4.0 International 

2. Recommended citation for this dataset: Nicholson Thomas, I., Black, B., & Kurmann, M. (2026). Bern Landscape Change Explorer Ecosystem Condition and Potential Results (Version 1.0.0) [Data set]. ETH Zurich. https://doi.org/10.5281/zenodo.18327368

DESCRIPTION/REPRODUCIBILITX: 
- This dataset contains rasters of ecosystem condition and potential for the Bern Landscape Change Explorer project.
- All rasters are projected in the CH1903+ / LV95 coordinate reference system and share the same extent and resolution (1 hectare).
- For detailed descriptions of the data generation process, raster value interpretation and suggested colour scheme, please refer to the accompanying metadata files:",
  paste("- ", basename(ecopot_metadata_files)),
  "",
  "This dataset contains raster data on ecosystem condition (EC) and the potential for condition dependent
improvements in selected ecosystem services. The calculations were performed
separately for three ecosystem types: forests, grasslands and alpine pastures, and agricultural land.
The data sets on ecosystem condition are stored as ec_index and separated by ecosystem type.
ec_index_LU_Forest contains the ecosystem condition index for forest areas,
ec_index_LU_GrassPasture for grasslands and alpine pastures, and ec_index_LU_IntensiveAg for
intensive agricultural areas. In all three cases, the pixel value corresponds to the calculated,
composite ecosystem condition index for the respective land use group.
In addition, there is the cluster_EC_all dataset, which represents spatial patterns in ecosystem
condition as classes. The clusters are based on the same indicator layers that were used to calculate
the ecosystem condition index. Each grid cell is assigned a cluster ID, which groups together cells
that have similar characteristics across several indicators, for example, simultaneously high or low
values in the same state characteristics. This allows areas to be delineated that are similar in terms
of the interaction of the state indicators considered.
The potential datasets are named after ecosystem services and end in _gap. CAR_gap describes the
potential for carbon storage, FF_gap for food and feed production, HAB_gap for habitat quality,
POL_gap for pollination, and REC_gap for recreational potential.
The pixel values indicate the difference between theoretical and effective provision, i.e., the
additional amount that could be achieved with an improved ecosystem status.",
  "",
  "Included files:",
  file_tree_ecopot
)

ecopot_readme_path <- file.path(ecopot_archive_dir, "README.txt")
writeLines(ecopot_metadata_content, con = ecopot_readme_path)
# path to markdown readme file
ecopot_md_path <- file.path(getwd(), "ecopot_readme.md")
# Convert markdown to plain text
pandoc_convert(
  input = ecopot_readme_path,
  to = "markdown",
  output = ecopot_md_path
)
