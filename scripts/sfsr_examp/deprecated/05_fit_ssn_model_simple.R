# FIT SOUTH FORK SALMON RIVER SSN MODEL
#
# Author(s): Mike Ackerman
#
# Purpose: Fit a preliminary spatial stream network model of August mean daily
#   temperature using NorWeST observations and predict at 1-km points.
#
# Source: https://usepa.github.io/SSN2/articles/introduction.html
#
# Created: October 15, 2026
#   Last modified:
#
# NOTES: This is a first test of the assembled SFSR SSN. StreamCat covariates
#   have not yet been added, so begin with an intercept-only spatial model.

# clear environment
rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(sf)
library(SSN2)

# -------------------
# IMPORT ASSEMBLED SSN
ssn_path = here("output", "ma_examp", "ssn", "sfsr_example.ssn")
stopifnot(dir.exists(ssn_path))

sfsr_ssn = ssn_import(
  path        = ssn_path,
  include_obs = TRUE,
  predpts     = "pred_sfsr",
  overwrite   = TRUE
)

# ----------------------
# EXPLORE THE SSN OBJECT
summary(sfsr_ssn)
names(sfsr_ssn)
ssn_names(sfsr_ssn)

# observation sites
head(sfsr_ssn$obs)
nrow(sfsr_ssn$obs)
names(sfsr_ssn$obs)

# prediction points (stored as a named list)
names(sfsr_ssn$preds$pred_sfsr)
head(sfsr_ssn$preds$pred_sfsr)

# stream network
head(sfsr_ssn$edges)
nrow(sfsr_ssn$edges)

# ---------------------------
# EXPLORE SPATIAL INFORMATION
st_crs(sfsr_ssn$edges)
st_crs(sfsr_ssn$obs)
st_crs(sfsr_ssn$preds$pred_sfsr)

# Check that the response and additive function values are available
stopifnot(
  "meanDailyMean" %in% names(sfsr_ssn$obs),
  "afvArea" %in% names(sfsr_ssn$obs),
  "afvArea" %in% names(sfsr_ssn$preds$pred_sfsr),
  sum(!is.na(sfsr_ssn$obs$meanDailyMean)) >= 3
)

# -----------------------
# PLOT OBSERVED TEMPERATURE
ggplot() +
  geom_sf(data = sfsr_ssn$edges, linewidth = 0.3, color = "steelblue") +
  geom_sf(data = sfsr_ssn$preds$pred_sfsr, color = "yellow", size = 0.5) +
  geom_sf(data = sfsr_ssn$obs, aes(color = meanDailyMax), size = 2) +
  scale_color_viridis_c(na.value = "grey70") +
  theme_bw() +
  labs(title = "Observed August Mean Daily Temperature",
       color = "Temperature (°C)")

# ------------------------------
# EXPLORE SSN NETWORK ATTRIBUTES

# SSN data contain spatial attributes describing locations relative to the stream network.
sfsr_ssn$obs |>
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

# ---------------------------
# EXPLORE A RESPONSE VARIABLE
summary(sfsr_ssn$obs$meanDailyMax) # the mean of daily maximum temperatures from site-year deployment during August 

# distribution of summer mean stream temps
response_p = sfsr_ssn$obs |>
  ggplot(aes(x = meanDailyMax)) +
  geom_histogram(bins  = 20, fill  = "steelblue", color = "white") +
  labs(
    x = "Mean Daily Max Aug Stream Temp (°C)",
    y = "Number of Sites",
  ) +
  theme_bw()
response_p

# ----------------------------------
# CREATE STREAM DISTANCE MATRICES
# Distances are needed between observations and between observations and predictions.
ssn_create_distmat(
  ssn.object    = sfsr_ssn,
  predpts       = "pred_sfsr",
  among_predpts = TRUE,     # needed if you want to do kriging later (FALSE only reduces computation time)
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
  formula    = meanDailyMax ~ elevws + precip9120ws,
  ssn.object = sfsr_ssn,
  type       = c("flowcon", "flowuncon", "euclid")
)
plot(tg)

# ----------------------------------------
# FIT EXAMPLE SPATIAL STREAM NETWORK MODEL

# Spatial covariance components:
#   tail-up   = covariance among flow-connected locations
#   tail-down = covariance among flow-connected AND flow-unconnected locations
#   Euclidean = spatial covariance not explicitly tied to network topology
#   nugget    = spatially independent/local variation
?ssn_lm
ssn_mod = ssn_lm(
  formula       = meanDailyMax ~ elevws + precip9120ws,
  ssn.object    = sfsr_ssn,
  tailup_type   = "exponential",
  taildown_type = "exponential",
  euclid_type   = "gaussian",
  additive      = "afvArea"
)

# ---------------------
# EXAMINE MODEL RESULTS
summary(ssn_mod)
coef(ssn_mod)
varcomp(ssn_mod)
tidy(ssn_mod, conf.int = TRUE)
glance(ssn_mod)

# -----------------------------
# COMPARE COVARIANCE STRUCTURES

# fit a simpler model containing only a tail-down covariance component
ssn_mod2 = ssn_lm(
  formula = meanDailyMax ~ elevws + precip9120ws,
  ssn.object = sfsr_ssn,
  taildown_type = "exponential"
)
glances(ssn_mod, ssn_mod2) # compare model-fit statistics

# NOTE: ssn_lm() uses REML by default. AIC/AICc comparisons among REML models are appropriate when the fixed-effects structure is the same and only the
# covariance structure differs. If fixed effects differ among models, fit with ML for likelihood-based comparisons such as AIC/AICc.

# --------------------------------
# COMPARE FIXED-EFFECT STRUCTURES
# use maximum likelihood (ML) when comparing models with different fixed-effect structures

# simplest
ml_mod1 = ssn_lm(
  formula       = meanDailyMax ~ 1,
  ssn.object    = sfsr_ssn,
  tailup_type   = "exponential",
  taildown_type = "exponential",
  euclid_type   = "gaussian",
  additive      = "afvArea",
  estmethod     = "ml"
)

# elevation only
ml_mod2 = ssn_lm(
  formula       = meanDailyMax ~ elevws,
  ssn.object    = sfsr_ssn,
  tailup_type   = "exponential",
  taildown_type = "exponential",
  euclid_type   = "gaussian",
  additive      = "afvArea",
  estmethod     = "ml"
)

# precip only
ml_mod3 = ssn_lm(
  formula       = meanDailyMax ~ precip9120ws,
  ssn.object    = sfsr_ssn,
  tailup_type   = "exponential",
  taildown_type = "exponential",
  euclid_type   = "gaussian",
  additive      = "afvArea",
  estmethod     = "ml"
)

# both
ml_mod4 = ssn_lm(
  formula       = meanDailyMax ~ elevws + precip9120ws,
  ssn.object    = sfsr_ssn,
  tailup_type   = "exponential",
  taildown_type = "exponential",
  euclid_type   = "gaussian",
  additive      = "afvArea",
  estmethod     = "ml"
)
glances(ml_mod1, ml_mod2, ml_mod3, ml_mod4)

# --------------------------
# LEAVE-ONE-OUT CROSS VALIDATION

# cross-validation provides another way to compare predictive performance
# lower RMSPE indicates better out-of-sample predictive performance
loocv_mod1 = loocv(ml_mod1)
loocv_mod2 = loocv(ml_mod2)
loocv_mod3 = loocv(ml_mod3)
loocv_mod4 = loocv(ml_mod4)
loocv_mod1$RMSPE
loocv_mod2$RMSPE
loocv_mod3$RMSPE
loocv_mod4$RMSPE

# mod2 had lowest AIC, but 4 has lowest RMSPE; I'll move forward with 2

# -----------------
# MODEL DIAGNOSTICS

# add fitted values, residuals, leverage, Cook's distance, etc.
aug_ssn_mod = augment(ml_mod2)
aug_ssn_mod
names(aug_ssn_mod)

# optional: save diagnostics as a geopackage
# st_write(
#   aug_ssn_mod,
#   here("output", "aug_ssn_mod.gpkg"),
#   delete_dsn = TRUE
# )

# standard diagnostic plot: fitted values vs standardized residuals
plot(ml_mod2, which = 1)

# SSN2 provides six standard diagnostic plots
# plot(ssn_mod, which = 1:6)

# --------------------------
# PREDICTION (KRIGING)

# predict temperature at the 1-km prediction sites
pred_1km = predict(ml_mod2, newdata = "pred_sfsr")
head(pred_1km)

# augment prediction points with model predictions
pred_1km = augment(ml_mod2, newdata = "pred_sfsr")
head(pred_1km)
pred_1km[, ".fitted"]

# map predicted temperatures
ggplot() +
  geom_sf(data = sfsr_ssn$edges, linewidth = 0.3, color = "grey60") +
  geom_sf(data = pred_1km, aes(color = .fitted), size = 1) +
  scale_color_viridis_c(option = "H") +
  labs(
    title = "Predicted Mean Daily Max Stream Temperature",
    color = "Temperature (°C)"
  ) +
  theme_bw()

# I should also get and plot prediction errors

# block prediction estimates the average response over a set of prediction locations rather than separate point predictions
predict(
  ml_mod2,
  newdata  = "pred_sfsr",
  block    = TRUE,
  interval = "prediction"
)

### END SCRIPT
