# PREPARE NSI HYDROGRAPHY
#
# Author(s): Mike Ackerman
#
# Purpose: Prep the NSI hydrography 17U (Pacific Northwest) dataset downloaded from: https://research.fs.usda.gov/rmrs/projects/national-stream-internet
#   We'll trim to SFSR using steelhead populations and prep for analysis.
#
# Created: October 2, 2026
#   Last modified: October 20, 2026
#
# NOTES: I'm using population polygons to establish a "study area", but any polygons (e.g., HUCs, watersheds, custom) could
#   be swapped and used instead.

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(sf)
library(here)

# --------------------------
# LOAD POPULATION BOUNDARIES

# load steelhead population polygons
load(here("data", "steelhead_pops.rda"))

# inspect
class(sth_pop)
names(sth_pop)
st_crs(sth_pop)

# all available populations
sort(unique(sth_pop$TRT_POPID))

# select populations of interest (used to establish study area)
pop_ids = c("SFMAI-s", "SFSEC-s")
pop_polys = sth_pop |>
  filter(TRT_POPID %in% pop_ids)

# map populations of interest
ggplot() +
  geom_sf(data = pop_polys, aes(fill = TRT_POPID)) +
  theme_bw()

# dissolve pop_polys into a single study area
study_area = pop_polys |>
  summarise()

# similar to...
# study_area = st_sf(
#   geometry = st_union(pop_polys)
# )

# map study area
ggplot() +
  geom_sf(data = study_area, fill = NA) +
  theme_bw()

# --------------
# LOAD NSI DATA
nsi_raw_dir = here("data", "NSI", "raw")
list.files(nsi_raw_dir, recursive = T)
list.files(
  nsi_raw_dir,
  pattern = "\\.shp$",
  recursive = TRUE,
  full.names = TRUE
)

# read NSI hydrography shapefile
nsi_flow = st_read(here(nsi_raw_dir, "VPU17_flowlines", "Flowline_PN17_NSI.shp"))

# inspect
dim(nsi_flow)
names(nsi_flow)
st_crs(nsi_flow)
st_crs(study_area)

# transform study area polygon to NSI CRS
study_area = st_transform(
  study_area,
  crs = st_crs(nsi_flow)
)
st_crs(study_area) == st_crs(nsi_flow)

# trim nsi to study area
nsi_flow_sfsr = nsi_flow |>
  st_filter(study_area)

# alternate method
#nsi_flow_sfsr = nsi_flow[st_intersects(nsi_flow, study_area, sparse = FALSE)[, 1],]

# how much did we trim?
nrow(nsi_flow)
nrow(nsi_flow_sfsr)

# map results
sfsr_p = ggplot() +
  geom_sf(data = nsi_flow_sfsr, linewidth = 0.3, color = "steelblue") +
  geom_sf(data = study_area, fill = NA, linewidth = 0.6) +
  theme_bw()
sfsr_p

# -----------------------
# SAVE PROCESSED NSI DATA

# create processed data directory, if needed
nsi_processed_dir = here("data", "NSI", "processed")
dir.create(
  nsi_processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# save selected NSI flowlines
st_write(
  nsi_flow_sfsr,
  file.path(nsi_processed_dir, "nsi_flow_sfsr.gpkg"),
  delete_dsn = TRUE
)

### END SCRIPT