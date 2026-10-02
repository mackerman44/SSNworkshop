# PREPARE NSI HYDROGRAPHY
#
# Author(s): Mike Ackerman
#
# Purpose: Prep the NSI hydrography 17U (Pacific Northwest) dataset downloaded from: https://research.fs.usda.gov/rmrs/projects/national-stream-internet
#   We'll trim to SFSR using steelhead populations and prep for analysis
#
# Created: October 2, 2026
#   Last modified:

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

sort(unique(sth_pop$TRT_POPID))   # all available populations
pop_ids = c("SFMAI-s", "SFSEC-s") # populations of interest
pop_polys = sth_pop |>
  filter(TRT_POPID %in% pop_ids)

# map
ggplot() +
  geom_sf(data = pop_polys, aes(fill = TRT_POPID)) +
  theme_bw()

# dissolve pop_polys into a single study area
study_area = pop_polys |>
  summarise()
study_area

# similar to...
# study_area = st_sf(
#   geometry = st_union(pop_polys)
# )

# map
ggplot() +
  geom_sf(data = study_area, fill = NA) +
  theme_bw()

# --------------
# LOAD NSI DATA
nsi_raw_dir = here("data", "hydrography", "NSI", "raw")
list.files(nsi_raw_dir, recursive = T)
list.files(
  nsi_raw_dir,
  pattern = "\\.shp$",
  recursive = TRUE,
  full.names = TRUE
)

# read flow and pred_point shapefiles
nsi_flow = st_read(here("data", "hydrography", "NSI", "raw", "VPU17_flowlines", "Flowline_PN17_NSI.shp"))
nsi_pred = st_read(here("data", "hydrography", "NSI", "raw", "VPU17_pred_points", "PredictionPoints_PN17_NSI.shp"))

# inspect
dim(nsi_flow)
dim(nsi_pred)
names(nsi_flow)
names(nsi_pred)
st_crs(nsi_flow)
st_crs(nsi_pred)
st_crs(study_area)

# transform study area polygon to NSI CRS
study_area = st_transform(
  study_area,
  st_crs(nsi_flow)
)
st_crs(study_area) == st_crs(nsi_flow)

# trim nsi to study area
nsi_flow_sfsr = nsi_flow[st_intersects(nsi_flow, study_area, sparse = FALSE)[, 1],]
nsi_pred_sfsr = nsi_pred[st_intersects(nsi_pred, study_area, sparse = FALSE)[, 1],]

# how much did we trim? PERFECT! down to 1335
nrow(nsi_flow)
nrow(nsi_flow_sfsr)
nrow(nsi_pred)
nrow(nsi_pred_sfsr)

# map results
sfsr_p = ggplot() +
  geom_sf(data = nsi_flow_sfsr, linewidth = 0.3, color = "steelblue") +
  geom_sf(data = nsi_pred_sfsr, size = 0.2, color = "pink") +
  geom_sf(data = study_area, fill = NA, linewidth = 0.6) +
  theme_bw()
sfsr_p

# -----------------------
# SAVE PROCESSED NSI DATA

# create processed data directory, if needed
nsi_processed_dir = here("data", "hydrography", "NSI", "processed")
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

# save selected NSI prediction points
st_write(
  nsi_pred_sfsr,
  file.path(nsi_processed_dir, "nsi_pred_sfsr.gpkg"),
  delete_dsn = TRUE
)

### END SCRIPT