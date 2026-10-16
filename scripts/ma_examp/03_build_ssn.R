# BUILD SOUTH FORK SALMON RIVER SSN
#
# Author(s): Mike Ackerman
#
# Purpose: Use SSNbler to create a Spatial Stream Network (SSN) object from
#   National Stream Internet (NSI) hydrography, NorWeST temperature observations,
#   and 1-km prediction points for the South Fork Salmon River.
#
# Source: https://pet221.github.io/SSNbler/articles/introduction.html
#
# Created: October 15, 2026
#   Last modified:
#
# NOTES: Use the NSI network trimmed to SFMAI-s and SFSEC-s steelhead populations.
#   StreamCat covariates have not yet been added. Observation and prediction
#   datasets were prepared in the preceding ma_examp scripts.

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
nsi_flow_sfsr = st_read(here("data", "hydrography", "NSI", "processed", "nsi_flow_sfsr.gpkg"))
sfsr_obs = st_read(here("data", "norwest", "processed", "sfsr_norwest_obs.gpkg"))
sfsr_pred = st_read(here("output", "ma_examp", "sfsr_pred_1km.gpkg"))

# ------------------------
# EXPLORE PREPPED INPUT DATA
dim(nsi_flow_sfsr)
dim(sfsr_obs)
dim(sfsr_pred)

# geometry types
st_geometry_type(nsi_flow_sfsr) |> table()
st_geometry_type(sfsr_obs) |> table()
st_geometry_type(sfsr_pred) |> table()

# coordinate reference systems
st_crs(nsi_flow_sfsr)
st_crs(sfsr_obs)
st_crs(sfsr_pred)

# variables
names(nsi_flow_sfsr)
names(sfsr_obs)
names(sfsr_pred)

# SSNbler requires:
#   streams          = LINESTRING geometry
#   sites            = POINT geometry
#   common projected = coordinate reference system
# Geographic coordinates (longitude/latitude) should not be used because
# SSNbler calculations rely on distances measured in map units.
stopifnot(
  !st_is_longlat(nsi_flow_sfsr),
  st_crs(nsi_flow_sfsr) == st_crs(sfsr_obs),
  st_crs(nsi_flow_sfsr) == st_crs(sfsr_pred),
  all(st_geometry_type(sfsr_obs) == "POINT"),
  all(st_geometry_type(sfsr_pred) == "POINT"),
  !anyDuplicated(sfsr_obs$OBSPRED_ID),
  !anyDuplicated(sfsr_pred$OBSPRED_ID),
  !any(sfsr_obs$OBSPRED_ID %in% sfsr_pred$OBSPRED_ID)
)

# inspect observation response availability
summary(sfsr_obs$meanDailyMean)
sum(!is.na(sfsr_obs$meanDailyMean))

# ----------------
# PLOT INPUT DATA
ggplot() +
  geom_sf(data = nsi_flow_sfsr, linewidth = 0.3, color = "steelblue") +
  geom_sf(data = sfsr_pred, size = 1, color = "pink") +
  geom_sf(data = sfsr_obs, size = 1, color = "red") +
  theme_bw() +
  labs(title = "SFSR SSN Inputs")

# ---------------------------------
# BUILD THE LANDSCAPE NETWORK (LSN)
# A LSN represents streams as directed edges connected by nodes. The direction
# of each edge represents downstream flow.
# Valid node types include:
#   source      = beginning of a headwater edge
#   outlet      = downstream end of a network
#   confluence  = multiple upstream edges converge
#   pseudonode  = one edge enters and one edge exits
# Correct topology is critical because SSN models use these relationships to
# determine hydrologic connectivity among locations.

# -------------------
# DEFINE LSN LOCATION
# SSNbler creates several intermediate files when building the LSN. Keep these
# separate from the original NSI hydrography.
lsn_path = here("output", "ssnbler", "ma_examp", "lsn")

# -------------------
# CREATE THE LSN
?lines_to_lsn

# convert MULTILINESTRING to LINESTRING for SSNbler
nsi_flow_ls_sfsr = st_cast(nsi_flow_sfsr, to = "LINESTRING")
nrow(nsi_flow_sfsr)
nrow(nsi_flow_ls_sfsr) # previously, each MULTILINESTRING contained one part
stopifnot(nrow(nsi_flow_sfsr) == nrow(nsi_flow_ls_sfsr))

# IMPORTANT: snap_tolerance and topo_tolerance depend on map units and network
# characteristics. The SSNbler vignette uses 0.05 and 20; try those for now.
edges = lines_to_lsn(
  streams        = nsi_flow_ls_sfsr,
  lsn_path       = lsn_path,
  check_topology = TRUE,
  snap_tolerance = 0.05,
  topo_tolerance = 20,
  overwrite      = TRUE
)

# -------------------------
# INSPECT LSN OUTPUT FILES
# lines_to_lsn() creates nodes.gpkg, edges.gpkg, nodexy.csv,
# noderelationships.csv, and relationships.csv.
# If node_errors.gpkg is generated, investigate before continuing.
list.files(lsn_path)
# Topology was previously checked; no additional resolution needed here.

nodes = st_read(file.path(lsn_path, "nodes.gpkg"))
table(nodes$nodecat)

# ------------------------
# VISUALIZE LSN TOPOLOGY
lsn_p = ggplot() +
  geom_sf(data = edges, color = "skyblue3", linewidth = 0.5) +
  geom_sf(data = nodes, aes(color = nodecat), size = 1.2) +
  theme_bw() +
  labs(title = "SFSR LSN", color = "Node Type")
lsn_p

# -------------------------------
# INCORPORATE SITES INTO THE LSN
# sites_to_lsn() snaps point locations to the nearest LSN edge and adds:
#   rid      = edge identifier where the site resides
#   ratio    = relative site position along the edge
#   snapdist = Euclidean distance the site was moved
# Run sites_to_lsn() separately for observations and prediction points.

# ----------------------
# ADD OBSERVATION POINTS
?sites_to_lsn
# NorWeST points may not coincide exactly with reconditioned NSI flowlines.
# Try 20 map units (meters, if CRS uses meters);
# review large snap distances and any unmatched sites before proceeding.
obs_lsn = sites_to_lsn(
  sites          = sfsr_obs,
  edges          = edges,
  lsn_path       = lsn_path,
  file_name      = "obs_sfsr",
  snap_tolerance = 20,
  save_local     = TRUE,
  overwrite      = TRUE
)
summary(obs_lsn$snapdist)
cat("Observations before snapping:", nrow(sfsr_obs), "\n")
cat("Observations after snapping:", nrow(obs_lsn), "\n")

# ----------------------
# ADD PREDICTION POINTS
# These points were sampled directly on the NSI network. Use a small tolerance
# initially, but inspect snapping results after LSN topology processing.
pred_lsn = sites_to_lsn(
  sites          = sfsr_pred,
  edges          = edges,
  lsn_path       = lsn_path,
  file_name      = "pred_sfsr",
  snap_tolerance = 1,
  save_local     = TRUE,
  overwrite      = TRUE
)
summary(pred_lsn$snapdist)
cat("Predictions before snapping:", nrow(sfsr_pred), "\n")
cat("Predictions after snapping:", nrow(pred_lsn), "\n")

# Pause if any sites were dropped; inspect their locations and snap distances.
# Do not simply increase tolerance without checking nearby stream branches.
stopifnot(nrow(obs_lsn) == nrow(sfsr_obs),
          nrow(pred_lsn) == nrow(sfsr_pred))

# ---------------------------
# CALCULATE UPSTREAM DISTANCE
# upDist is hydrologic distance from the network outlet to each feature.
# For edges, upDist represents distance from the outlet to the upstream end.
# For sites, upDist represents distance from the outlet to the site position.

# ------------------------
# UPSTREAM DISTANCE: EDGES
?updist_edges
edges = updist_edges(
  edges       = edges,
  save_local  = TRUE,
  lsn_path    = lsn_path,
  calc_length = TRUE
)
summary(edges$Length)
summary(edges$upDist)

# ------------------------
# UPSTREAM DISTANCE: SITES
?updist_sites
site_list = updist_sites(
  sites      = list(obs_sfsr = obs_lsn, pred_sfsr = pred_lsn),
  edges      = edges,
  length_col = "Length",
  save_local = TRUE,
  lsn_path   = lsn_path
)
summary(site_list$obs_sfsr$upDist)
summary(site_list$pred_sfsr$upDist)

# ----------------------------
# MAP UPSTREAM DISTANCE
ggplot() +
  geom_sf(data = edges, aes(color = upDist), linewidth = 0.5) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(title = "Upstream Distance", color = "Upstream\nDistance")

# -----------------------------------------
# CALCULATE ADDITIVE FUNCTION VALUES (AFVs)
# AFVs are used by tail-up covariance models to weight the relative influence
# of upstream branches at confluences.
# SSNbler calculates segment proportional influence (PI) and AFVs.
# A suitable influence variable should be numeric, complete, and preferably
# strictly positive. Cumulative watershed area is a common flow surrogate.

# -------------------------
# IDENTIFY INFLUENCE COLUMN
names(edges)
stopifnot("TotDASqKM" %in% names(edges))
summary(edges$TotDASqKM)
sum(is.na(edges$TotDASqKM))
sum(edges$TotDASqKM <= 0, na.rm = TRUE)

# Stop if watershed area is missing or nonpositive; investigate before AFVs.
stopifnot(all(is.finite(edges$TotDASqKM)), all(edges$TotDASqKM > 0))

# ----------------
# AFVs FOR EDGES
?afv_edges
edges = afv_edges(
  edges     = edges,
  infl_col  = "TotDASqKM",
  segpi_col = "areaPI",
  afv_col   = "afvArea",
  lsn_path  = lsn_path
)
summary(edges$afvArea)

# ----------------
# AFVs FOR SITES
?afv_sites
site_list = afv_sites(
  sites      = site_list,
  edges      = edges,
  afv_col    = "afvArea",
  save_local = TRUE,
  lsn_path   = lsn_path
)
summary(site_list$obs_sfsr$afvArea)
summary(site_list$pred_sfsr$afvArea)

# ------------------------
# ASSEMBLE THE SSN OBJECT
# ssn_assemble() converts the processed LSN into the .ssn structure used by SSN2.
# Observations include NorWeST August temperature summaries for subsequent models.
?ssn_assemble
ssn_path = here("output", "ma_examp", "ssn", "sfsr_example.ssn")
sfsr_ssn = ssn_assemble(
  edges      = edges,
  lsn_path   = lsn_path,
  obs_sites  = site_list$obs_sfsr,
  preds_list = site_list["pred_sfsr"],
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
summary(sfsr_ssn$obs$meanDailyMean)

# The assembled SSN object contains:
#   $edges = processed stream network
#   $obs   = observation sites and August temperature summaries
#   $preds = named list of prediction datasets
#   $path  = location of the .ssn directory
# ssn_assemble() also creates netID, pid, locID, and netgeom attributes.

# -------------------------
# PLOT ASSEMBLED SSN OBJECT
ggplot() +
  geom_sf(data = sfsr_ssn$edges, linewidth = 0.3) +
  geom_sf(data = sfsr_ssn$preds$pred_sfsr, size = 1, color = "pink") +
  geom_sf(data = sfsr_ssn$obs, aes(color = meanDailyMean), size = 2) +
  scale_color_viridis_c(na.value = "grey70") +
  theme_bw() +
  labs(title = "SFSR Spatial Stream Network", color = "August Mean\nTemperature (°C)")

### END SCRIPT
