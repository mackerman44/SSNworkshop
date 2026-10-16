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
names(sfsr_ssn$preds)
nrow(sfsr_ssn$obs)
nrow(sfsr_ssn$preds$pred_sfsr)
summary(sfsr_ssn$obs$meanDailyMean)

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
  geom_sf(data = sfsr_ssn$obs, aes(color = meanDailyMean), size = 1.5) +
  scale_color_viridis_c(na.value = "grey70") +
  theme_bw() +
  labs(title = "Observed August Mean Daily Temperature",
       color = "Temperature (°C)")

# ----------------------------------
# CREATE STREAM DISTANCE MATRICES
# Distances are needed between observations and between observations and predictions.
ssn_create_distmat(
  ssn.object    = sfsr_ssn,
  predpts       = "pred_sfsr",
  among_predpts = FALSE,
  overwrite     = TRUE
)

# ----------------------------------------
# FIT PRELIMINARY SPATIAL STREAM MODEL
# Intercept-only: estimate overall mean and spatial covariance, without
# landscape predictors. Tail-up covariance uses watershed-area AFVs.
ssn_mod = ssn_lm(
  formula       = meanDailyMean ~ 1,
  ssn.object    = sfsr_ssn,
  tailup_type   = "exponential",
  taildown_type = "exponential",
  additive      = "afvArea"
)

# ----------------------
# EXAMINE MODEL RESULTS
summary(ssn_mod)
coef(ssn_mod)
varcomp(ssn_mod)

# --------------------------
# PREDICTION (KRIGING)
# Predict August mean daily temperature at prepared prediction points.
pred_1km = augment(ssn_mod, newdata = "pred_sfsr")
head(pred_1km)
summary(pred_1km$.fitted)

# map predictions
ggplot() +
  geom_sf(data = sfsr_ssn$edges, linewidth = 0.3, color = "grey60") +
  geom_sf(data = pred_1km, aes(color = .fitted), size = 1) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(title = "Predicted August Mean Daily Temperature",
       color = "Temperature (°C)")

### END SCRIPT
