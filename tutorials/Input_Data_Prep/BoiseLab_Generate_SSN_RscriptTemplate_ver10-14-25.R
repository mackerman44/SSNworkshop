################################################################
#
# Forest Service, RMRS, Boise Lab R script template for generating SSNs
# using the SSNbler R package
#
# Nagel, ver. 10/14/25
# Contact: david.nagel@usda.gov
# 
# This script can be used to generate an SSN using data
# formatting procedures developed at the FS RMRS Boise Lab.
# This version incorporates parallel processing to overcome previous
# limitations to the number of stream reaches.
#
# This procedure assumes that the user is will be utilizing the
# following inputs: Streamlines, observation points, and prediction
# points.
#
# For more information visit this website:
# https://research.fs.usda.gov/rmrs/projects/spatial-stream-network
#
################################################################

## ------- Install SSNbler and sf (if necessary) and load library ------

## Install packages from CRAN repository

install.packages("sf")
install.packages("SSNbler")

## Install two additional packages needed for SSNbler

install.packages(c("doParallel", "foreach"),
                 dependencies = TRUE)

## Note: If using RStudio, go to Packages tab and be sure SSNbler
## and sf are checked on.

## ------------- Load Libraries and Define Cores -------------

## Load libraries

library(SSNbler)
library(sf)

## Load two more libraries for parallel processing

library(tidyverse)
library(parallel)

## Set number of cores to use for parallel processing

detectCores() # How many cores do I have?
num.cores <- 10  ## Don't all the cores. Leave a few to do other things.


###############################

# ------------- Set file names and paths ----------------

## Note: Update paths, input, and output file names to reflect
## your personal preferences.

## Set the working directory, containing the input topologically
## corrected streamlines, observation points, and prediction
## points. 

setwd("D:/FolderName/FolderName/Foldername")

## Import the streams and point datasets as sf Objects

stream_net<- st_read("Streams_TopoCorrected_ForSSN.shp")
obs_points <- st_read("Observation_Points_ForSSN.shp")
pred_points <- st_read("Prediction_Points_ForSSN.shp")

## Set path and a name for LSN and SSN output folders

lsn_path1<- "D:/FolderName/FolderName/Foldername/LSN1"
ssn_path1<- "D:/FolderName/FolderName/Foldername/SSN1"

## Cast the streamline shapefile as linestring format

stream_net2 = st_cast(stream_net,"LINESTRING")  

# ------------- Check Network Topology ----------------

## Imports the strealines to edges
## Generate the LSN that will hold the lines, sites, and preds
## Check topology

edges<- lines_to_lsn(
  streams = stream_net2,
  lsn_path = lsn_path1,
  snap_tolerance = 0.09,
  check_topology = TRUE,
  topo_tolerance = 20,
  overwrite = TRUE,
  verbose = TRUE,
  remove_ZM = TRUE,
  use_parallel = TRUE,
  no_cores = num.cores)

# ------------- Generate SSN ----------------

## Import the observation sites

obs <- sites_to_lsn(
  sites = obs_points,
  edges = edges,
  snap_tolerance = 1,
  save_local = TRUE,
  lsn_path = lsn_path1,
  file_name = "sites.gpkg",
  overwrite = TRUE
)

## Import the prediction points

preds <- sites_to_lsn(
  sites = pred_points,
  edges = edges,
  snap_tolerance = 1,
  save_local = TRUE,
  lsn_path = lsn_path1,
  file_name = "preds.gpkg",
  overwrite = TRUE
)

# Generate upstream distance for the edges

edges <- updist_edges(
  edges = edges,
  lsn_path = lsn_path1,
  calc_length = TRUE,
  length_col = "Length",
  save_local = TRUE,
  overwrite = TRUE,
)

## Generate the sites list. This is a "named" list

sitelist = list(obs = obs,preds = preds)
names(sitelist)

## Compute updist on sites and preds

site.list <- updist_sites(
  sites = sitelist,
  edges = edges,
  length_col = "Length",
  lsn_path = lsn_path1,
  save_local = TRUE,
  overwrite = TRUE
)

## Generate additive function value for edges

## Note that a field name for weighting is required
## In this instance the field name is called
## TotDASqKM and represents total upstream drainage area
## for each reach in the dataset.

edges <- afv_edges(
  edges = edges,
  lsn_path = lsn_path1,
  infl_col = "TotDASqKM",
  # infl_col = "CUMDRAINAG",
  segpi_col = "areaPI",
  afv_col = "afvArea",
  save_local = TRUE,
  overwrite = TRUE
)

## Generate additive function value for sites and preds
## using the list

site.list<- afv_sites(
  sites = site.list,
  edges = edges,
  afv_col = "afvArea",
  save_local = TRUE,
  lsn_path = lsn_path1,
  overwrite = TRUE
)

## List for preds

names(site.list)

## Generate the SSN object

Output_ssn <- ssn_assemble(
  edges = edges,
  lsn_path = lsn_path1,
  obs_sites = site.list$obs, ##Put name of site.list element after $
  preds_list = site.list["preds"],
  ssn_path = ssn_path1,
  import = TRUE,
  check = TRUE,
  afv_col = "afvArea",
  overwrite = TRUE
)




