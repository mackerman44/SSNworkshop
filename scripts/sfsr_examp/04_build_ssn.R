# BUILD SOUTH FORK SALMON RIVER SSN
#
# Author(s): Mike Ackerman
#
# Purpose: Use SSNbler to build a Spatial Stream Network (SSN) from NSI
#   hydrography, NorWeST temperature observations, and two sets of
#   prediction points for the South Fork Salmon River.
#
# Source: https://pet221.github.io/SSNbler/articles/introduction.html
#
# Created: October 15, 2026
#   Last modified: October 21, 2026
#
# NOTES: NSI flowlines have six StreamCat covariates attached.
#   Prediction datasets include original NSI prediction points and
#   custom approximately 1-km prediction points.
#   Preserve all observation records for subsequent model selection.

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(sf)
library(here)
library(SSN2)
library(SSNbler)

# -------------------
# LOAD PROCESSED DATA

nsi_flow_sfsr = st_read(here("data", "NSI", "processed", "nsi_flow_sfsr_streamcat.gpkg"))
obs_sfsr = st_read(here("data", "norwest", "processed", "sfsr_norwest_obs.gpkg"))

# load both prediction datasets as a named list
pred_list = list(
  pred_nsi_sfsr = st_read(here("output", "sfsr_examp", "preds_nsi_sfsr.gpkg")),
  pred_1km_sfsr = st_read(here("output", "sfsr_examp", "preds_1km_sfsr.gpkg"))
)

# --------------------------------
# EXPLORE PREPARED INPUT DATA

dim(nsi_flow_sfsr)
dim(obs_sfsr)
lapply(pred_list, dim)

# coordinate reference systems
st_crs(nsi_flow_sfsr)
st_crs(obs_sfsr)
lapply(pred_list, st_crs)

# variables
names(nsi_flow_sfsr)
names(obs_sfsr)
lapply(pred_list, names)

# verify projected CRS, point geometry, and unique site IDs
all_ids = c(
  obs_sfsr$OBSPRED_ID,
  unlist(lapply(pred_list, \(x) x$OBSPRED_ID))
)

stopifnot(
  !st_is_longlat(nsi_flow_sfsr),
  st_crs(nsi_flow_sfsr) == st_crs(obs_sfsr),
  all(vapply(
    pred_list,
    \(x) st_crs(x) == st_crs(nsi_flow_sfsr),
    logical(1)
  )),
  all(st_geometry_type(obs_sfsr) == "POINT"),
  all(vapply(
    pred_list,
    \(x) all(st_geometry_type(x) == "POINT"),
    logical(1)
  )),
  !anyNA(all_ids),
  !anyDuplicated(all_ids)
)

# inspect observation response availability
summary(obs_sfsr$meanDailyMean)
summary(obs_sfsr$meanDailyMax)

# --------------------------------
# IDENTIFY STREAMCAT COVARIATES

# all six StreamCat covariates attached in script 03
sc_covars = c(
  "elevcat", "elevws",
  "precip9120cat", "precip9120ws",
  "tmean9120cat", "tmean9120ws"
)

stopifnot(all(sc_covars %in% names(nsi_flow_sfsr)))

# ----------------
# PLOT INPUT DATA

ggplot() +
  geom_sf(data = nsi_flow_sfsr, linewidth = 0.3, color = "steelblue") +
  geom_sf(data = pred_list$pred_nsi_sfsr, size = 0.7, color = "orange") +
  geom_sf(data = pred_list$pred_1km_sfsr, size = 0.7, color = "purple") +
  geom_sf(data = obs_sfsr, size = 1.2, color = "red") +
  theme_bw() +
  labs(title = "SFSR SSN Inputs")

# ---------------------------------
# BUILD THE LANDSCAPE NETWORK (LSN)

# define output locations
lsn_path = here("output", "sfsr_examp", "lsn")
ssn_path = here("output", "sfsr_examp", "ssn", "sfsr_examp.ssn")

# convert MULTILINESTRING to LINESTRING
nsi_flow_ls_sfsr = st_cast(
  nsi_flow_sfsr,
  to = "LINESTRING",
  warn = FALSE
)

# verify geometry conversion preserved the number of records
nrow(nsi_flow_sfsr) == nrow(nsi_flow_ls_sfsr)

# create the landscape stream network
edges = lines_to_lsn(
  streams        = nsi_flow_ls_sfsr,
  lsn_path       = lsn_path,
  check_topology = TRUE,
  snap_tolerance = 0.05,
  topo_tolerance = 20,
  overwrite      = TRUE
)

# inspect output files and node types
list.files(lsn_path)

nodes = st_read(file.path(lsn_path, "nodes.gpkg"))
table(nodes$nodecat)

# visualize LSN topology
ggplot() +
  geom_sf(data = edges, color = "skyblue3", linewidth = 0.5) +
  geom_sf(data = nodes, aes(color = nodecat), size = 1.2) +
  theme_bw() +
  labs(title = "SFSR Landscape Stream Network", color = "Node Type")

# -------------------------------
# INCORPORATE SITES INTO THE LSN

# snap observations to the network; NorWeST locations may differ from the reconditioned NSI flowlines.
obs_lsn = sites_to_lsn(
  sites          = obs_sfsr,
  edges          = edges,
  lsn_path       = lsn_path,
  file_name      = "obs_sfsr",
  snap_tolerance = 20,
  save_local     = TRUE,
  overwrite      = TRUE
)

# snap each prediction dataset separately; use a small tolerance because predictions originate from NSI hydrography.
pred_lsn = lapply(names(pred_list), function(nm) {
  sites_to_lsn(
    sites          = pred_list[[nm]],
    edges          = edges,
    lsn_path       = lsn_path,
    file_name      = nm,
    snap_tolerance = 2,
    save_local     = TRUE,
    overwrite      = TRUE
  )
})

names(pred_lsn) = names(pred_list)

# inspect snapping distances
summary(obs_lsn$snapdist)
lapply(pred_lsn, \(x) summary(x$snapdist))

# compare site counts before and after snapping
cat("Observations:", nrow(obs_sfsr), "->", nrow(obs_lsn), "\n")

for (nm in names(pred_list)) {
  cat(
    nm, ":",
    nrow(pred_list[[nm]]), "->",
    nrow(pred_lsn[[nm]]), "\n"
  )
}

# stop if any sites were dropped
stopifnot(
  nrow(obs_lsn) == nrow(obs_sfsr),
  all(vapply(
    names(pred_list),
    \(nm) nrow(pred_lsn[[nm]]) == nrow(pred_list[[nm]]),
    logical(1)
  ))
)

# --------------------------------
# CALCULATE UPSTREAM DISTANCE

# calculate upstream distance for edges
edges = updist_edges(
  edges       = edges,
  save_local  = TRUE,
  lsn_path    = lsn_path,
  calc_length = TRUE
)

summary(edges$Length)
summary(edges$upDist)

# combine observation and prediction sites into one named list
site_list = c(
  list(obs_sfsr = obs_lsn),
  pred_lsn
)

# calculate upstream distance for all sites
site_list = updist_sites(
  sites      = site_list,
  edges      = edges,
  length_col = "Length",
  save_local = TRUE,
  lsn_path   = lsn_path
)

# inspect upstream distances
lapply(site_list, \(x) summary(x$upDist))

# map upstream distance
ggplot() +
  geom_sf(data = edges, aes(color = upDist), linewidth = 0.5) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(
    title = "Upstream Distance — SFSR",
    color = "Upstream\nDistance"
  )

# -----------------------------------------
# CALCULATE ADDITIVE FUNCTION VALUES (AFVs)

# verify watershed area is available and positive
stopifnot(
  "TotDASqKM" %in% names(edges),
  all(is.finite(edges$TotDASqKM)),
  all(edges$TotDASqKM > 0)
)

# calculate AFVs for edges
edges = afv_edges(
  edges     = edges,
  infl_col  = "TotDASqKM",
  segpi_col = "areaPI",
  afv_col   = "afvArea",
  lsn_path  = lsn_path
)

summary(edges$afvArea)

# calculate AFVs for all sites
site_list = afv_sites(
  sites      = site_list,
  edges      = edges,
  afv_col    = "afvArea",
  save_local = TRUE,
  lsn_path   = lsn_path
)

# inspect site AFVs
lapply(site_list, \(x) summary(x$afvArea))

# -------------------------------------
# ATTACH STREAMCAT COVARIATES TO SITES

# extract StreamCat covariates from network edges
edge_covars = edges |>
  st_drop_geometry() |>
  select(rid, all_of(sc_covars))

# verify one record per edge
stopifnot(!anyDuplicated(edge_covars$rid))

# attach covariates to observations and both prediction datasets
site_list = lapply(
  site_list,
  \(x) left_join(
    x,
    edge_covars,
    by = "rid",
    relationship = "many-to-one"
  )
)

# inspect covariate coverage for each site dataset
lapply(site_list, \(x) {
  x |>
    st_drop_geometry() |>
    summarise(
      n_sites = n(),
      across(
        all_of(sc_covars),
        ~ sum(is.na(.x)),
        .names = "missing_{.col}"
      )
    )
})

# ------------------------------------
# INSPECT OBSERVATIONS AND PREDICTIONS

# inspect missing temperature responses
site_list$obs_sfsr |>
  st_drop_geometry() |>
  summarise(
    n_total = n(),
    n_missing_mean = sum(!is.finite(meanDailyMean)),
    n_missing_max = sum(!is.finite(meanDailyMax))
  )

# preserve all observation records at this stage
# filter observations for individual model formulas in script 05

# inspect completeness of prediction covariates
lapply(site_list[names(pred_list)], \(x) {
  x |>
    st_drop_geometry() |>
    summarise(
      across(
        all_of(sc_covars),
        ~ sum(!is.finite(.x)),
        .names = "missing_{.col}"
      )
    )
})

# ------------------------
# ASSEMBLE THE SSN OBJECT

# assemble observations and both prediction datasets
sfsr_ssn = ssn_assemble(
  edges      = edges,
  lsn_path   = lsn_path,
  obs_sites  = site_list$obs_sfsr,
  preds_list = site_list[names(pred_list)],
  ssn_path   = ssn_path,
  import     = TRUE,
  check      = TRUE,
  afv_col    = "afvArea",
  overwrite  = TRUE
)

# ----------------------
# EXPLORE THE SSN OBJECT

class(sfsr_ssn)
names(sfsr_ssn)
summary(sfsr_ssn)

names(sfsr_ssn$edges)
names(sfsr_ssn$obs)
names(sfsr_ssn$preds)

# inspect observations
summary(sfsr_ssn$obs$meanDailyMean)
summary(sfsr_ssn$obs$meanDailyMax)

# inspect prediction datasets
lapply(sfsr_ssn$preds, dim)

# inspect StreamCat covariates
lapply(
  c(list(obs_sfsr = sfsr_ssn$obs), sfsr_ssn$preds),
  \(x) names(x)[names(x) %in% sc_covars]
)

# -------------------------
# PLOT ASSEMBLED SSN OBJECT

ggplot() +
  geom_sf(data = sfsr_ssn$edges, linewidth = 0.3, color = "grey60") +
  geom_sf(data = sfsr_ssn$preds$pred_nsi_sfsr, size = 0.7, color = "orange") +
  geom_sf(data = sfsr_ssn$preds$pred_1km_sfsr, size = 0.7, color = "purple") +
  geom_sf(data = sfsr_ssn$obs, aes(color = meanDailyMean), size = 1.5) +
  scale_color_viridis_c(na.value = "grey70") +
  theme_bw() +
  labs(
    title = "SFSR Spatial Stream Network",
    color = "August Mean\nTemperature (°C)"
  )

### END SCRIPT