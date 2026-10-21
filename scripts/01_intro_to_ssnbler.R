# SSNbler INTRODUCTION
#
# Author(s): Mike Ackerman
#
# Purpose: Learn the SSNbler workflow for creating a Spatial Stream Network (SSN) object using National Stream Internet (NSI) hydrography.
#   Here, I use the NSI trmmed to the South Fork Salmon River (SFMAI-s and SFSEC-s steelhead pops) prepped in 02_prep_nsi_4_sfsr.R
#
# Source: https://pet221.github.io/SSNbler/articles/introduction.html
#
# Created: October 2, 2026
#   Last modified:

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
nsi_processed_dir = here("data", "hydrography", "NSI", "processed")
nsi_flow_sfsr = st_read(file.path(nsi_processed_dir, "nsi_flow_sfsr.gpkg")) # NSI stream network
nsi_pred_sfsr = st_read(file.path(nsi_processed_dir, "nsi_pred_sfsr.gpkg")) # NSI prediction points

# ------------------------
# EXPLORE PREPPED NSI DATA
dim(nsi_flow_sfsr)
dim(nsi_pred_sfsr)
# geometry types
st_geometry_type(nsi_flow_sfsr) |> table()
st_geometry_type(nsi_pred_sfsr) |> table()
# coordinate reference systems
st_crs(nsi_flow_sfsr)
st_crs(nsi_pred_sfsr)
# variables
names(nsi_flow_sfsr)
names(nsi_pred_sfsr)

# SSNbler requires:
#   streams          = LINESTRING geometry
#   sites            = POINT geometry
#   common projected = coordinate reference system
#
# Geographic coordinates (longitude/latitude) should not be used because SSNbler calculations rely on distances measured in map units.

# ----------------
# PLOT INPUT DATA
ggplot() +
  geom_sf(data = nsi_flow_sfsr, linewidth = 0.3, color = "steelblue") +
  geom_sf(data = nsi_pred_sfsr, size = 0.5, color = "pink") +
  theme_bw()

# ---------------------------------
# BUILD THE LANDSCAPE NETWORK (LSN)

# A LSN represents streams as directed edges connected by nodes. The direction of each edge represents downstream flow.
# Valid node types include:
#   source      = beginning of a headwater edge
#   outlet      = downstream end of a network
#   confluence  = multiple upstream edges converge
#   pseudonode  = one edge enters and one edge exits
#
# Correct topology is critical because SSN models use these relationships to determine hydrologic connectivity among locations.

# -------------------
# DEFINE LSN LOCATION

# SSNbler creates several intermediate files when building the LSN. Keep these separate from the original NSI hydrography.
lsn_path = here("output", "ssnbler", "sfsr_example", "lsn")

# -------------------
# CREATE THE LSN
?lines_to_lsn

# convert MULTILINESTRING to LINESTRING for SSNbler
nsi_flow_ls_sfsr = st_cast(nsi_flow_sfsr, to = "LINESTRING")
nrow(nsi_flow_sfsr)
nrow(nsi_flow_ls_sfsr) # great, each MULTILINESTRING contained only one part

# IMPORTANT: snap_tolerance and topo_tolerance depend on the map units and characteristics
#   of the stream network. Inspect the CRS and NSI documentation before choosing final values.
#
# The SSNbler Middle Fork vignette uses:
#   snap_tolerance = 0.05
#   topo_tolerance = 20
#
# ...so we'll try those for now
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

# lines_to_lsn() creates at least:
#   nodes.gpkg
#   edges.gpkg
#   nodexy.csv
#   noderelationships.csv
#   relationships.csv
#
# if potential topology errors are detected: node_errors.gpkg, do not continue until potential topology errors have been investigated.

list.files(lsn_path)

# read in nodes
nodes = st_read(file.path(lsn_path, "nodes.gpkg"))
table(nodes$nodecat)

# ------------------------
# VISUALIZE LSN TOPOLOGY
lsn_p =  ggplot() +
  geom_sf(data = edges, color = "skyblue3", linewidth = 0.5) +
  geom_sf(data = nodes, aes(color = nodecat), size = 1.2) +
  theme_bw() +
  labs(
    title = "SFSR LSN",
    color = "Node Type"
  )
lsn_p

# -------------------------------
# INCORPORATE SITES INTO THE LSN
# sites_to_lsn() snaps point locations to the nearest LSN edge and adds:
#   rid      = edge identifier where the site resides
#   ratio    = relative site position along the edge
#   snapdist = Euclidean distance the site was moved
#
# sites_to_lsn() must be run separately for every observation and prediction dataset, even when the points already intersect stream edges.

# ----------------------
# ADD PREDICTION POINTS
?sites_to_lsn
pred_sfsr = sites_to_lsn(
  sites          = nsi_pred_sfsr,
  edges          = edges,
  lsn_path       = lsn_path,
  file_name      = "pred_sfsr",
  snap_tolerance = 0.05,
  save_local     = TRUE,
  overwrite      = TRUE
)
summary(pred_sfsr$snapdist) # check snapping distances

# ---------------------------
# CALCULATE UPSTREAM DISTANCE
# upDist is hydrologic distance from the network outlet to each feature. For edges, upDist represents distance from the outlet to the upstream end of the edge.
#
# For sites, upDist represents distance from the outlet to the site's position on its associated edge.

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
  sites = list(pred_sfsr = pred_sfsr),
  edges      = edges,
  length_col = "Length",
  save_local = TRUE,
  lsn_path   = lsn_path
)
summary(site_list$pred_sfsr$upDist)

# ----------------------------
# MAP UPSTREAM DISTANCE
ggplot() +
  geom_sf(data = edges, aes(color = upDist), linewidth = 0.5) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(
    title = "Upstream Distance",
    color = "Upstream\nDistance"
  )

# -----------------------------------------
# CALCULATE ADDITIVE FUNCTION VALUES (AFVs)

# AFVs are used by tail-up covariance models to weight the relative influence of upstream branches at confluences.
#
# SSNbler calculates:
#   1. segment proportional influence (PI)
#   2. additive function values (AFVs)
#
# SSN2 subsequently uses the AFVs to calculate spatial weights when fitting tail-up models.
# A suitable influence variable should:
#   - be numeric
#   - exist for every edge
#   - contain no missing values
#   - preferably contain no zeros
#
# Cumulative watershed area is commonly used as a surrogate for flow volume.

# -------------------------
# IDENTIFY INFLUENCE COLUMN

# inspect candidate NSI attributes
names(edges)

# Once a suitable cumulative watershed-area variable has been identified:
summary(edges$TotDASqKM)
sum(is.na(edges$TotDASqKM))
sum(edges$TotDASqKM == 0)

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

# AFVs should range from 0 to 1.
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
summary(site_list$pred_sfsr$afvArea)

# ------------------------
# ASSEMBLE THE SSN OBJECT

# ssn_assemble() converts the processed LSN into the .ssn structure used by SSN2.
# Observed sites are optional, so we can initially assemble this example using the stream network and prediction locations only.
?ssn_assemble
ssn_path = here("output", "ssnbler", "sfsr_example", "ssn", "sfsr_example.ssn")
sfsr_ssn = ssn_assemble(
  edges      = edges,
  lsn_path   = lsn_path,
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
names(sfsr_ssn$preds)

# The assembled SSN object contains:
#   $edges = processed stream network
#   $obs   = observation sites (NA when none are supplied)
#   $preds = named list of prediction datasets
#   $path  = location of the .ssn directory
#
# ssn_assemble() also creates important network attributes including:
#   netID
#   pid
#   locID
#   netgeom

# -------------------------
# PLOT ASSEMBLED SSN OBJECT
ggplot() +
  geom_sf(data = sfsr_ssn$edges, linewidth = 0.3) +
  geom_sf(data = sfsr_ssn$preds$pred_sfsr, size = 0.5) +
  theme_bw() +
  labs(title = "SFSR Spatial Stream Network")

### END SCRIPT