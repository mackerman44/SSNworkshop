# SSN2 INTRODUCTION
#
# Author(s): Mike Ackerman
#
# Purpose: Prep and example code for the Spatial Stream Network (SSN) workshop, October 14-16, 2026, Boise
#
#   Script follows examples from the SSN2 introductory vignette:
#   https://usepa.github.io/SSN2/articles/introduction.html
#
# Created: October 1, 2026
#   Last modified:

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(SSN2)
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
?MiddleFork04.ssn # metadata for the Middle Fork example dataset

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
    x = "Mean Summer Stream Temperature (°C)",
    y = "Number of Sites",
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

# ----------------------------------
# CREATE STREAM DISTANCE MATRICES

# SSN models require hydrologic distance matrices that identify the direction and topology of the stream network
# by default, these are written to the "distance" directory inside the .ssn folder.
?ssn_create_distmat
ssn_create_distmat(
  ssn.object    = mf04p,
  predpts       = c("pred1km", "CapeHorn"),
  among_predpts = TRUE,
  overwrite     = TRUE
)

# -------------------------------
# EXPLORE SPATIAL AUTOCORRELATION

# torgegram distinguishes spatial relationships unique to stream networks:
#   flowcon   = flow-connected sites
#   flowuncon = flow-unconnected sites
#   euclid    = ordinary Euclidean spatial relationships
?Torgegram
tg = Torgegram(
  formula    = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04p,
  type       = c("flowcon", "flowuncon", "euclid")
)
plot(tg)


# ----------------------------------------
# FIT EXAMPLE SPATIAL STREAM NETWORK MODEL

# model mean summer stream temperature as a function of:
#   ELEV_DEM  = elevation
#   AREAWTMAP = area-weighted precipitation
#
# Spatial covariance components:
#   tail-up   = covariance among flow-connected locations
#   tail-down = covariance among flow-connected AND flow-unconnected locations
#   Euclidean = spatial covariance not explicitly tied to network topology
#   nugget    = spatially independent/local variation
?ssn_lm
ssn_mod = ssn_lm(
  formula       = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object    = mf04p,
  tailup_type   = "exponential",
  taildown_type = "spherical",
  euclid_type   = "gaussian",
  additive      = "afvArea"
)

# ----------------------
# EXAMINE MODEL RESULTS
summary(ssn_mod)               # model summary
coef(ssn_mod)                  # fixed-effect coefficients
varcomp(ssn_mod)               # proportion of variability attributed to fixed effects and each covariance component
tidy(ssn_mod, conf.int = TRUE) # tidy fixed-effect estimates and confidence intervals
glance(ssn_mod)                # overall model-fit statistics

# --------------------------------
# COMPARE COVARIANCE STRUCTURES

# fit a simpler model containing only a tail-down covariance component
ssn_mod2 = ssn_lm(
  formula       = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object    = mf04p,
  taildown_type = "spherical"
)
glances(ssn_mod, ssn_mod2) # compare model-fit statistics

# NOTE: ssn_lm() uses REML by default. AIC/AICc comparisons among REML models are appropriate when the fixed-effects structure is the same and only the
# covariance structure differs. If fixed effects differ among models, fit with ML for likelihood-based comparisons such as AIC/AICc.

# --------------------------------
# COMPARE FIXED-EFFECT STRUCTURES

# use maximum likelihood (ML) when comparing models with different fixed-effect structures
ml_mod = ssn_lm(
  formula       = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object    = mf04p,
  tailup_type   = "exponential",
  taildown_type = "spherical",
  euclid_type   = "gaussian",
  additive      = "afvArea",
  estmethod     = "ml"
)

# reduced fixed-effects model without elevation
ml_mod2 = ssn_lm(
  formula       = Summer_mn ~ AREAWTMAP,
  ssn.object    = mf04p,
  tailup_type   = "exponential",
  taildown_type = "spherical",
  euclid_type   = "gaussian",
  additive      = "afvArea",
  estmethod     = "ml"
)

glances(ml_mod, ml_mod2)

# --------------------------
# LEAVE-ONE-OUT CROSS VALIDATION

# cross-validation provides another way to compare predictive performance
loocv_mod = loocv(ssn_mod)
loocv_mod$RMSPE

loocv_mod2 = loocv(ssn_mod2)
loocv_mod2$RMSPE
# lower RMSPE indicates better out-of-sample predictive performance

# -----------------
# MODEL DIAGNOSTICS

# add fitted values, residuals, leverage, Cook's distance, etc.
aug_ssn_mod = augment(ssn_mod)
head(aug_ssn_mod)
names(aug_ssn_mod)

# standard diagnostic plot: fitted values vs standardized residuals
plot(ssn_mod, which = 1)

# SSN2 provides six standard diagnostic plots
# plot(ssn_mod, which = 1:6)

# optional: save diagnostics as a geopackage
# st_write(
#   aug_ssn_mod,
#   here("output", "aug_ssn_mod.gpkg"),
#   delete_dsn = TRUE
# )

# ----------------------
# PREDICTION (KRIGING)

# predict summer mean temperature at the 1-km prediction sites
pred_1km = predict(
  ssn_mod,
  newdata = "pred1km"
)
head(pred_1km)

# augment prediction points with model predictions
aug_preds = augment(
  ssn_mod,
  newdata = "pred1km"
)
head(aug_preds)

# map predicted summer mean stream temperature
mf_pred_p = ggplot() +
  geom_sf(data = mf04p$edges) +
  geom_sf(data = aug_preds, aes(color = .fitted), size = 2) +
  scale_color_viridis_c(option = "H") +
  labs(
    title = "Predicted Summer Mean Stream Temperature",
    color = "Temperature (°C)"
  ) +
  theme_bw()
mf_pred_p

# predictions can also be made for all prediction datasets
# predict(ssn_mod)
# predict(ssn_mod, newdata = "all")

# block prediction estimates the average response over a set of prediction locations rather than separate point predictions
predict(
  ssn_mod,
  newdata = "pred1km",
  block = TRUE,
  interval = "prediction"
)

# =====================================
# ADVANCED LINEAR SSN MODEL FEATURES
# =====================================

# ---------------------------
# FIX COVARIANCE PARAMETERS
# example: fix the Euclidean partial sill at 1
euclid_init = euclid_initial(
  "gaussian",
  de = 1,
  known = "de"
)
euclid_init

ssn_init = ssn_lm(
  formula        = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object     = mf04p,
  tailup_type    = "exponential",
  taildown_type  = "spherical",
  euclid_initial = euclid_init,
  additive       = "afvArea"
)
ssn_init

# --------------
# RANDOM EFFECTS (error here)

# random intercept for each stream network 
ssn_rand = ssn_lm(
  formula       = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object    = mf04p,
  tailup_type   = "exponential",
  taildown_type = "spherical",
  euclid_type   = "gaussian",
  additive      = "afvArea",
  random        = ~ as.factor(netID)
)
ssn_rand


# -----------------
# PARTITION FACTORS
# partition factors alter the covariance structure so observations in different groups are treated as uncorrelated
ssn_part = ssn_lm(
  formula          = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object       = mf04p,
  tailup_type      = "exponential",
  taildown_type    = "spherical",
  euclid_type      = "gaussian",
  additive         = "afvArea",
  partition_factor = ~ as.factor(netID)
)
ssn_part

# ===========================================
# GENERALIZED LINEAR SPATIAL STREAM MODELS
# ===========================================
summary(mf04p$obs$C16) # C16 = number of summer days stream temperature exceeded 16°C

# map C16
c16_p = ggplot() +
  geom_sf(data = mf04p$edges) +
  geom_sf(data = mf04p$obs, aes(color = C16), size = 2) +
  scale_color_viridis_c(option = "H") +
  labs(
    title = "Number of Days Exceeding 16°C",
    color = "Days"
  ) +
  theme_bw()
c16_p

# ------------------
# POISSON SSN MODEL
ssn_pois = ssn_glm(
  formula       = C16 ~ ELEV_DEM + AREAWTMAP,
  family        = "poisson",
  ssn.object    = mf04p,
  tailup_type   = "epa",
  taildown_type = "mariah",
  additive      = "afvArea"
)
summary(ssn_pois)

# --------------------------
# NEGATIVE BINOMIAL SSN MODEL

# negative binomial models can accommodate overdispersion in count data
ssn_nb = ssn_glm(
  formula       = C16 ~ ELEV_DEM + AREAWTMAP,
  family        = "nbinomial",
  ssn.object    = mf04p,
  tailup_type   = "epa",
  taildown_type = "mariah",
  additive      = "afvArea"
)
summary(ssn_nb)

# compare predictive performance
loocv_pois = loocv(ssn_pois)
loocv_pois$RMSPE

loocv_nb = loocv(ssn_nb)
loocv_nb$RMSPE

# ==========================
# SIMULATE DATA ON A NETWORK
# ==========================

# specify covariance parameters
tu_params = tailup_params(
  "exponential",
  de = 0.4,
  range = 1e5
)

td_params = taildown_params(
  "spherical",
  de = 0.1,
  range = 1e6
)

euc_params = euclid_params(
  "gaussian",
  de = 0.2,
  range = 1e3
)

nug_params = nugget_params(
  "nugget",
  nugget = 0.1
)

# ------------------------
# SIMULATE GAUSSIAN DATA
set.seed(2)
sims_gaussian = ssn_simulate(
  family          = "gaussian",
  ssn.object      = mf04p,
  network         = "obs",
  additive        = "afvArea",
  tailup_params   = tu_params,
  taildown_params = td_params,
  euclid_params   = euc_params,
  nugget_params   = nug_params,
  mean            = 0,
  samples         = 1
)
head(sims_gaussian)

# ------------------------
# SIMULATE BINOMIAL DATA
set.seed(2)
sims_binomial = ssn_simulate(
  family          = "binomial",
  ssn.object      = mf04p,
  network         = "obs",
  additive        = "afvArea",
  tailup_params   = tu_params,
  taildown_params = td_params,
  euclid_params   = euc_params,
  nugget_params   = nug_params,
  mean            = 0,
  samples         = 2
)
head(sims_binomial)

### END SCRIPT
