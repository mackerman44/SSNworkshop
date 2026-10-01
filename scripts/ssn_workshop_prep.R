# SSN2 WORKSHOP PREPARATION
#
# Author(s): Mike Ackerman
#
# Purpose: Prep and example code for the Spatial Stream Network (SSN) workshop, October 14-16, 2026, Boise,
#
#   Script follows examples from the SSN2 introductory vignette:
#   https://usepa.github.io/SSN2/articles/introduction.html
#
# Created: October 1, 2026
#   Last modified:

# clear environment
rm(list = ls())

# install packages, if needed (lines only need to be run once, so leave commented during normal use)
# install.packages("SSN2")
# install.packages("SSNbler")

# view package citations
citation("SSN2")
citation("SSNbler")

# Load packages
library(tidyverse)
library(SSN2)   
library(SSNbler)    
library(sf)         
library(here)
library(janitor)

# -------------------
# PROJECT DIRECTORIES

# root directory of the R project
here() 

# path to example middle fork .ssn dataset stored within this project
mf_ssn_path = here("data", "ssn", "MiddleForkTemperature", "MiddleForkTemperature.ssn")

# ensure that the directory exists
dir.exists(mf_ssn_path)


# ---------------
# IMPORT SSN DATA

# An SSN object contains:
#   $edges = stream network
#   $obs   = observation sites
#   $preds = prediction sites
#   $path  = path to the .ssn directory

# import the Middle Fork 2004 stream temperature dataset
?ssn_import
mf04p = ssn_import(
  path      = mf_ssn_path,
  predpts   = c("pred1km", "CapeHorn"),
  overwrite = TRUE
)

# ----------------------
# EXPLORE THE SSN OBJECT
summary(mf04p)   # overall summary
names(mf04p)     # examine the structure of the object
ssn_names(mf04p) # names of variables in observations and prediction datasets
mf04p$path       # path associated with the SSN object

# observation sites
head(mf04p$obs)
nrow(mf04p$obs)
names(mf04p$obs)
# View(mf04p$obs) # interactive view in RStudio, if desired

# prediction sites (stored as a named list)
names(mf04p$preds)
head(mf04p$preds$pred1km)
head(mf04p$preds$CapeHorn)
nrow(mf04p$preds$pred1km)
nrow(mf04p$preds$CapeHorn)

# stream network 
head(mf04p$edges)
nrow(mf04p$edges)
names(mf04p$edges)

# ----------------------------
# EXPLORE SPATIAL INFORMATION

# all three components are sf objects
class(mf04p$edges)
class(mf04p$obs)
class(mf04p$preds$pred1km)

# geometry types
st_geometry_type(mf04p$edges)
st_geometry_type(mf04p$obs) |> janitor::tabyl()
st_geometry_type(mf04p$preds$pred1km) |> janitor::tabyl()

# coordinate reference systems
st_crs(mf04p$edges)
st_crs(mf04p$obs)
st_crs(mf04p$preds$pred1km)

# NOTE: the observation layer in this copy of the example dataset imports without CRS metadata. The coordinates are in the same coordinate system as the
# stream network, so assign the network CRS to the observation layer. st_crs() assigns CRS metadata; it does NOT transform coordinates.
if (is.na(st_crs(mf04p$obs))) {
  st_crs(mf04p$obs) = st_crs(mf04p$edges)
}
st_crs(mf04p$obs) # confirm CRS

# -----------------------
# PLOT THE STREAM NETWORK
mf_p <- ggplot() +
  geom_sf(data = mf04p$edges) +
  geom_sf(data = mf04p$preds$pred1km, shape = 17, color = "blue") +
  geom_sf(data = mf04p$obs, color = "brown", size = 2) +
  labs(
    title = "Middle Fork 2004 SSN Example",
    subtitle = "Observation and 1-km prediction sites"
  ) +
  theme_bw()

mf_p

# ------------------------------
# EXPLORE SSN NETWORK ATTRIBUTES

# SSN data contain spatial attributes describing locations relative to the stream network.
mf04p$obs |>
  st_drop_geometry() |>
  select(netID, rid, pid, locID, upDist, ratio, netgeom) |>
  head()

# Key fields:
#   netID   = stream network identifier
#   rid     = reach/edge identifier
#   pid     = point/measurement identifier
#   locID   = unique spatial location identifier
#   upDist  = stream distance upstream from the network outlet
#   ratio   = relative position of a site along its edge
#   netgeom = protected representation of network topology/location

# -----------------------------
# EXPLORE THE RESPONSE VARIABLE
summary(mf04p$obs$Summer_mn) # Summer_mn = mean summer stream temperature (degrees C)

# distribution of summer mean stream temps
mn_summer_p = mf04p$obs |>
  ggplot(aes(x = Summer_mn)) +
  geom_histogram(
    bins  = 20,
    fill  = "steelblue",
    color = "white"
  ) +
  labs(
    x     = "Mean Summer Stream Temperature (°C)",
    y     = "Number of Sites",
  ) +
  theme_bw()

# plot summer mean temperature at observation sites
mf_temp_p = ggplot() +
  # stream network
  geom_sf(data = mf04p$edges) +
  # observation sites
  geom_sf(data = mf04p$obs, aes(color = Summer_mn), size = 2) +
  # set color scale; note min and max values (max below max observed in mf04p$obs$Summer_mn)
  scale_color_viridis_c(limits = c(0, 17), option = "H") +
  labs(
    title = "Summer Mean Stream Temperature",
    color = "Temperature (C)"
  ) +
  theme_bw()

mf_temp_p

### CONTINUE HERE

# ------------------------------------------------------------------------------
# 9. CREATE STREAM DISTANCE MATRICES
# ------------------------------------------------------------------------------

# SSN models require hydrologic distance matrices that preserve the
# direction and topology of the stream network.
#
# These are written to a "distance" directory inside the .ssn folder.

ssn_create_distmat(
  ssn.object = mf04p,
  predpts = c("pred1km", "CapeHorn"),
  among_predpts = TRUE,
  overwrite = TRUE
)


# ------------------------------------------------------------------------------
# 10. EXPLORE SPATIAL AUTOCORRELATION
# ------------------------------------------------------------------------------

# A Torgegram is analogous to a semivariogram but distinguishes spatial
# relationships unique to stream networks:
#
#   flowcon   = flow-connected sites
#   flowuncon = flow-unconnected sites
#   euclid    = ordinary Euclidean spatial relationships

tg <- Torgegram(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  type = c("flowcon", "flowuncon", "euclid")
)

plot(tg)


# ------------------------------------------------------------------------------
# 11. FIT AN EXAMPLE SPATIAL STREAM NETWORK MODEL
# ------------------------------------------------------------------------------

# Model mean summer stream temperature as a function of:
#
#   ELEV_DEM  = elevation
#   AREAWTMAP = area-weighted precipitation
#
# Spatial covariance components:
#
#   tail-up   = covariance among flow-connected locations
#   tail-down = covariance among flow-connected AND flow-unconnected locations
#   Euclidean = spatial covariance not explicitly tied to network topology
#   nugget    = spatially independent/local variation

ssn_mod <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_type = "gaussian",
  additive = "afvArea"
)


# ------------------------------------------------------------------------------
# 12. EXAMINE MODEL RESULTS
# ------------------------------------------------------------------------------

summary(ssn_mod)

# Fixed-effect coefficients
coef(ssn_mod)

# Additional model exploration can continue here during/after the workshop.


# ------------------------------------------------------------------------------
# 13. HELP AND DOCUMENTATION
# ------------------------------------------------------------------------------

# Metadata for the Middle Fork example dataset
?MiddleFork04.ssn

# Function documentation
?ssn_import
?ssn_create_distmat
?Torgegram
?ssn_lm
