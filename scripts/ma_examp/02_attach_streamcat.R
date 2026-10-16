# ATTACH STREAMCAT COVARIATES TO NSI
#
# Author(s): Mike Ackerman
#
# Purpose: Retrieve StreamCat landscape covariates for South Fork Salmon River NSI flowlines.
#
# Created: October 15, 2026
#   Last modified:
#
# NOTES: StreamCat API requests are temporarily paused. Uncomment this script
#        when ready to revisit the API; check returned columns before joining.

# clear environment
rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(sf)
library(StreamCatTools)

# --------------
# LOAD NSI DATA
sfsr_nsi = st_read(here("data", "hydrography", "NSI", "processed", "nsi_flow_sfsr.gpkg"))

# ----------------------------------
# ATTACH STREAMCAT COVARIATES TO NSI

# 1. Browse available metrics
sc_vars = sc_get_params(param = "variable_info")

# search metrics by keywords
sc_candidates = sc_vars |>
  filter(
    str_detect(
      str_to_lower(paste(category, metric, short_description, long_description)),
      "elev|forest|precip|baseflow|riparian|drainage"
    )
  ) |>
  select(category, metric, short_description, units, aoi, year)

View(sc_candidates)

# 2. Identify valid metric names
sc_metrics = sc_get_params(param = "metric_names")

# search valid names for selected variables
grep(
  "elev|pctconif|precip9120|Precip_Minus_EVT|pctfire",
  unlist(sc_metrics),
  value = TRUE,
  ignore.case = TRUE
)

# 3. Retrieve selected metrics for comids

# select valid metric names from step 2
selected_metrics = c(
  "elev",
  "precip9120"
)

# unique NSI COMIDs
nsi_comids = sfsr_nsi |>
  st_drop_geometry() |>
  filter(!is.na(COMID)) |>
  distinct(COMID) |>
  pull(COMID)

# retrieve metrics
# NOTE: Start with a small number of COMIDs to test API response and rate limits.
sc_data = sc_get_data(
  comid = head(nsi_comids, 5),
  metric = selected_metrics,
  aoi = "ws"
)

# inspect streamcat data
glimpse(sc_data)

# once the test works, retrieve the full set (possibly in batches)
sc_data = sc_get_data(
  comid = nsi_comids,
  metric = selected_metrics,
  aoi = "ws"
)

# join StreamCat covariates to NSI
# NOTE: Verify column names and one row per COMID in sc_data first.
sfsr_nsi_sc = sfsr_nsi |>
  mutate(COMID = as.character(COMID)) |>
  left_join(
    sc_data |>
      rename(COMID = comid) |>
      mutate(COMID = as.character(COMID)),
    by = "COMID",
    relationship = "many-to-one"
  )

# verify
nrow(sfsr_nsi)
nrow(sfsr_nsi_sc)

# write enriched NSI layer only after the join is verified
st_write(
  sfsr_nsi_sc,
  here("data", "hydrography", "NSI", "processed", "nsi_flow_sfsr_streamcat.gpkg"),
  delete_dsn = TRUE
)

### END SCRIPT
