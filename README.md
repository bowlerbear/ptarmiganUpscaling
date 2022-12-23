# ptarmiganUpscaling

Combining different data types to produce a national population estimate for the ptarmigan in Norway.

The repository contains the following files:

00_general functions.R - bunch of small functions used later, mainly formatting, processing stuff

01_formatting_occurrence_data.R - script to process and format the GBIF data for analysis
01_GBIF_data_retreival.R - script to retreive the GBIF data

02_environmental_data_grid.R - script to process all the environmental covariates for the 5 x 5 km grid
02_environmental_data_transect_buffers.R - scripts to process the covariates around each line-transect

03_mapping_linetransects_to_grids.R - script to align the occurrence data grid to the line-transect surveys

04_HPC_abundance_analysis.R - main script for analysis of the line-transect surveys
04_HPC_occurrence analysis.R - main script fot analysis of the occurrence data

05_crossvalidation_folds.R - script to define the spatial region for each fold

06_model_summaries.R - script to analysis the model outputs of both the line-transect abundance and occurrence data

07_HPC_simple_combined_analysis.R - script to combine the predictions from both models and predict total population size

08_uncertainity_analysis.R - script to examine the links between grid-level and national level uncertainty


Folders:

model: JAGS model code for both the abundance and occurrence models

model-outputs: model outputs (after running on the HPC) for both datasets (not pushed to GitHub)

ms: development of the paper (not pushed to GitHub)

plots: plots for the paper (not pushed to Github)

data: all kinds of data files - raw data as well as spatial environmental information (not pushed to Github)

old: a bunch of old files that can be ignored (not pushed to Github)
