# PREPARE PREDICTION POINTS
#
# Author(s): Mike Ackerman
#
# Purpose: Prepare alternative prediction point datasets for the South Fork
#   Salmon River using NSI-supplied points or custom points generated along
#   processed NSI hydrography.
#
# NSI-supplied points: https://research.fs.usda.gov/rmrs/projects/national-stream-internet#download-data
#
# Created: October 15, 2026
#   Last modified: October 20, 2026
#
# NOTES: Prepare two independent prediction datasets:
#   1. NSI prediction points downloaded with the original NSI hydrography.
#   2. Custom prediction points spaced approximately 1 km apart.
#   Custom points are spaced within individual NSI segments, not by
#   network-wide distance. At least one point is generated per segment.
#   Prediction IDs must not overlap with observation IDs.

# clear environment
rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(sf)

# ----------------------
# LOAD SPATIAL REFERENCE

# load processed NSI hydrography and study area
nsi_flow_sfsr = st_read(here("data", "NSI", "processed", "nsi_flow_sfsr.gpkg"))
study_area = st_read(here("data", "study_areas", "study_area_sfsr.gpkg"))

# load processed NorWeST observations
obs_sfsr = st_read(here("data", "norwest", "processed", "sfsr_norwest_obs.gpkg"))

# confirm all coordinate reference systems match
stopifnot(
  st_crs(study_area) == st_crs(nsi_flow_sfsr),
  st_crs(obs_sfsr)   == st_crs(nsi_flow_sfsr)
)

# confirm projected CRS with meter units
st_crs(nsi_flow_sfsr)
stopifnot(!st_is_longlat(nsi_flow_sfsr))
# NOTE: Confirm CRS linear units are meters before using 1000 below.

# -------------------------------
# OPTION 1: NSI PREDICTION POINTS

# load original NSI prediction points
nsi_pred_pts = st_read(here("data", "NSI", "raw", "VPU17_pred_points", "PredictionPoints_PN17_NSI.shp")
) |>
  st_transform(st_crs(nsi_flow_sfsr))

# inspect original NSI prediction point attributes
names(nsi_pred_pts)
nrow(nsi_pred_pts)

# trim NSI prediction points to study area
preds_nsi_sfsr = nsi_pred_pts |>
  st_filter(study_area)

# inspect
nrow(preds_nsi_sfsr)
head(preds_nsi_sfsr)

# ---------------------------------------
# OPTION 2: CUSTOM 1-KM PREDICTION POINTS

# convert multipart flowlines to individual lines
nsi_lines = nsi_flow_sfsr |>
  st_cast("LINESTRING", warn = FALSE) |>
  # calculate length (m) of each linestring
  mutate(length_m = as.numeric(st_length(geom))) |>
  # remove any linestrings that may have 0 length
  filter(as.numeric(st_length(geom)) > 0) |>
  # determine number of points per linestring; retain at least one prediction point per segment
  mutate(
    n_pts = pmax(1L, as.integer(round(length_m / 1000)))
  )

# sample points at regular intervals within each segment
# st_line_sample returns MULTIPOINT geometry
preds_1km_sfsr = nsi_lines

st_geometry(preds_1km_sfsr) = st_line_sample(
  st_geometry(preds_1km_sfsr),
  n = preds_1km_sfsr$n_pts,
  type = "regular"
)

# convert to individual points and retain COMID
preds_1km_sfsr = preds_1km_sfsr |>
  st_cast("POINT", warn = FALSE) |>
  select(COMID)

# ---------------------------
# ASSIGN PREDICTION POINT IDS

# NOTE: all OBSPRED_IDs need to be unique. I'm just retaining the original OBSPRED_IDs for the observation dataset
#   starting the first prediction IDs at 100000 and the second prediction IDs at 200000

# examine maximum observation ID
max_obs_id = max(obs_sfsr$OBSPRED_ID, na.rm = TRUE)

# assign unique IDs to NSI prediction points starting at 100000
preds_nsi_sfsr = preds_nsi_sfsr |>
  mutate(OBSPRED_ID = row_number() + 99999L)

# assign unique IDs to custom 1-km prediction points starting at 200000
preds_1km_sfsr = preds_1km_sfsr |>
  mutate(OBSPRED_ID = row_number() + 199999L)

# verify IDs are unique across all datasets
all_ids = c(
  obs_sfsr$OBSPRED_ID,
  preds_nsi_sfsr$OBSPRED_ID,
  preds_1km_sfsr$OBSPRED_ID
)
stopifnot(!anyNA(all_ids), !anyDuplicated(all_ids))

# inspect ID ranges
range(obs_sfsr$OBSPRED_ID)
range(preds_nsi_sfsr$OBSPRED_ID)
range(preds_1km_sfsr$OBSPRED_ID)

# -------------------------
# INSPECT PREDICTION POINTS

# compare numbers of prediction points
tibble(
  dataset = c("NSI supplied", "Custom 1 km"),
  n_points = c(nrow(preds_nsi_sfsr), nrow(preds_1km_sfsr))
)

# map alternative prediction datasets
ggplot() +
  geom_sf(data = nsi_flow_sfsr, color = "steelblue", linewidth = 0.3) +
  geom_sf(data = preds_nsi_sfsr, color = "purple", size = 1) +
  geom_sf(data = preds_1km_sfsr, color = "orange", size = 1) +
  geom_sf(data = obs_sfsr, color = "black", size = 2) + 
  geom_sf(data = study_area, fill = NA, linewidth = 0.6) +
  theme_bw() +
  labs(title = "NSI and Custom Prediction Points — SFSR")

# --------------------------------
# SAVE PROCESSED PREDICTION POINTS

# save NSI-supplied prediction points
st_write(
  preds_nsi_sfsr,
  here("output", "sfsr_examp", "preds_nsi_sfsr.gpkg"),
  delete_dsn = TRUE
)

# save custom 1-km prediction points
st_write(
  preds_1km_sfsr,
  here("output", "sfsr_examp", "preds_1km_sfsr.gpkg"),
  delete_dsn = TRUE
)

### END SCRIPT
