# PREPARE PREDICTION POINTS
#
# Author(s): Mike Ackerman
#
# Purpose: Create approximately 1-km prediction points on South Fork Salmon River NSI hydrography.
#
# Created: October 15, 2026
#   Last modified:
#
# NOTES: Use the processed NSI flowlines and assign OBSPRED_ID starting at 100000.
#        Points are spaced within individual NSI segments, not by network-wide distance.

# clear environment
rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(sf)

# --------------
# LOAD NSI DATA
sfsr_nsi = st_read(here("data", "hydrography", "NSI", "processed", "nsi_flow_sfsr.gpkg"))

# ------------------------
# CREATE PREDICTION POINTS

# ensure network is in a projected CRS with meter units
st_crs(sfsr_nsi)
stopifnot(!st_is_longlat(sfsr_nsi))
# Note: verify CRS linear units are meters before using 1000 below.

# convert multipart lines to individual lines, if necessary
nsi_lines = sfsr_nsi |>
  st_cast("LINESTRING", warn = FALSE) |>
  filter(as.numeric(st_length(geom)) > 0)

# generate points at approximately 1-km intervals
# use at least one point per segment so short NSI reaches aren't omitted
sfsr_pred = nsi_lines |>
  mutate(
    length_m = as.numeric(st_length(geom)),
    n_pts = pmax(1L, as.integer(round(length_m / 1000)))
  )

# st_line_sample returns a MULTIPOINT geometry for each line
st_geometry(sfsr_pred) = st_line_sample(
  st_geometry(sfsr_pred),
  n = sfsr_pred$n_pts,
  type = "regular"
)

sfsr_pred = sfsr_pred |>
  st_cast("POINT", warn = FALSE) |>
  select(COMID) |>
  # assign unique prediction IDs starting at 100000
  mutate(OBSPRED_ID = row_number() + 99999L)

# inspect
nrow(sfsr_pred)
head(sfsr_pred)
stopifnot(!anyDuplicated(sfsr_pred$OBSPRED_ID))

# check prediction IDs don't overlap with observations
sfsr_obs = st_read(here("data", "norwest", "processed", "sfsr_norwest_obs.gpkg"))
stopifnot(!any(sfsr_pred$OBSPRED_ID %in% sfsr_obs$OBSPRED_ID))

# plot prediction points and observation locations
ggplot() +
  geom_sf(data = sfsr_nsi, color = "steelblue", linewidth = 0.3) +
  geom_sf(data = sfsr_pred, color = "orange", size = 0.5) +
  geom_sf(data = sfsr_obs, color = "red", size = 0.8) +
  theme_bw()

# write prediction points to file
st_write(
  sfsr_pred,
  here("output", "ma_examp", "sfsr_pred_1km.gpkg"),
  delete_dsn = TRUE
)

### END SCRIPT
