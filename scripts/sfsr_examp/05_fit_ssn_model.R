# FIT SOUTH FORK SALMON RIVER SSN MODEL
#
# Author(s): Mike Ackerman
#
# Purpose: Evaluate fixed effects and spatial covariance structures for
#   August mean daily maximum stream temperature, and predict at two
#   sets of South Fork Salmon River prediction points.
#
# Source: https://usepa.github.io/SSN2/articles/introduction.html
#
# Created: October 15, 2026
#   Last modified: October 21, 2026
#
# NOTES: Compare fixed effects using ML and covariance structures using
#   REML with identical observations. Refit the selected model using
#   REML and predict at NSI and custom 1-km prediction points.
#
#   Candidate predictors include catchment and watershed summaries of
#   elevation, precipitation, and long-term air temperature.

# clear environment
rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(sf)
library(SSN2)

# -------------------
# IMPORT ASSEMBLED SSN
ssn_path = here("output", "sfsr_examp", "ssn", "sfsr_examp.ssn")
stopifnot(dir.exists(ssn_path))

# import observations and both prediction datasets
sfsr_ssn = ssn_import(
  path        = ssn_path,
  include_obs = TRUE,
  predpts     = c("pred_nsi_sfsr", "pred_1km_sfsr"),
  overwrite   = TRUE
)

# ----------------------
# EXPLORE THE SSN OBJECT
summary(sfsr_ssn)
names(sfsr_ssn)
ssn_names(sfsr_ssn)

# observation sites
nrow(sfsr_ssn$obs)
names(sfsr_ssn$obs)

# prediction datasets
names(sfsr_ssn$preds)
lapply(sfsr_ssn$preds, dim)

# stream network
nrow(sfsr_ssn$edges)
names(sfsr_ssn$edges)

# ---------------------------
# IDENTIFY RESPONSE AND PREDICTORS

response = "meanDailyMax"

sc_covars = c(
  "elevcat", "elevws",
  "precip9120cat", "precip9120ws",
  "tmean9120cat", "tmean9120ws"
)

# verify required fields
stopifnot(
  response %in% names(sfsr_ssn$obs),
  all(sc_covars %in% names(sfsr_ssn$obs)),
  "afvArea" %in% names(sfsr_ssn$obs),
  all(vapply(
    sfsr_ssn$preds,
    \(x) all(sc_covars %in% names(x)) &&
      "afvArea" %in% names(x),
    logical(1)
  ))
)

# -----------------------
# EXPLORE OBSERVED TEMPERATURE

summary(sfsr_ssn$obs[[response]])

ggplot(sfsr_ssn$obs, aes(x = meanDailyMax)) +
  geom_histogram(bins = 20, color = "white") +
  theme_bw() +
  labs(
    title = "Observed August Stream Temperature",
    x = "Mean Daily Maximum Temperature (°C)",
    y = "Number of Observations"
  )

# map observations and both prediction datasets
ggplot() +
  geom_sf(data = sfsr_ssn$edges, color = "grey70", linewidth = 0.3) +
  geom_sf(data = sfsr_ssn$preds$pred_nsi_sfsr, color = "orange", size = 0.5) +
  geom_sf(data = sfsr_ssn$preds$pred_1km_sfsr, color = "purple", size = 0.5) +
  geom_sf(data = sfsr_ssn$obs, aes(color = meanDailyMax), size = 1.5) +
  scale_color_viridis_c(na.value = "grey70") +
  theme_bw() +
  labs(
    title = "SFSR SSN Observations and Prediction Sites",
    color = "August\nTemperature (°C)"
  )

# --------------------------------
# EXPLORE PREDICTOR CORRELATIONS

cor_mat = sfsr_ssn$obs |>
  st_drop_geometry() |>
  select(all_of(sc_covars)) |>
  cor(use = "pairwise.complete.obs")

round(cor_mat, 2)

# plot correlations
cor_df = as.data.frame(as.table(cor_mat)) |>
  rename(
    variable_1 = Var1,
    variable_2 = Var2,
    correlation = Freq
  )

ggplot(cor_df, aes(variable_1, variable_2, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(correlation, 2)), size = 3) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "StreamCat Predictor Correlations — Observations",
    x = NULL,
    y = NULL,
    fill = "Pearson r"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ----------------------------------------
# PREPARE A CONSISTENT MODELING DATASET

# use the same complete observations for all candidate models this avoids comparing AIC values across different sample sizes
obs_complete = sfsr_ssn$obs |>
  st_drop_geometry() |>
  filter(
    if_all(
      all_of(c(response, sc_covars, "afvArea")),
      is.finite
    )
  )

cat("Original observations:", nrow(sfsr_ssn$obs), "\n")
cat("Complete observations:", nrow(obs_complete), "\n")
stopifnot(nrow(obs_complete) >= 10)

# retain only complete observation records in the modeling SSN prediction points remain unchanged
sfsr_ssn$obs = sfsr_ssn$obs |> filter(OBSPRED_ID %in% obs_complete$OBSPRED_ID)
stopifnot(nrow(sfsr_ssn$obs) == nrow(obs_complete))

# ----------------------------------
# CREATE STREAM DISTANCE MATRICES

# calculate distances for both prediction datasets
# among_predpts = FALSE avoids potentially large prediction-by-prediction
# distance matrices; these are unnecessary for individual predictions
ssn_create_distmat(
  ssn.object    = sfsr_ssn,
  predpts       = c("pred_nsi_sfsr", "pred_1km_sfsr"),
  among_predpts = FALSE,
  overwrite     = TRUE
)

# -------------------------------
# EXPLORE SPATIAL AUTOCORRELATION

# preliminary Torgegram with watershed elevation and precipitation
tg = Torgegram(
  formula = meanDailyMax ~ elevws + precip9120ws,
  ssn.object = sfsr_ssn,
  type = c("flowcon", "flowuncon", "euclid")
)

plot(tg)

# ------------------------------------
# COMPARE FIXED-EFFECT STRUCTURES (ML)

# begin with interpretable models and avoid including catchment and watershed versions of the same metric in a single candidate model
fixed_formulas = list(
  intercept       = meanDailyMax ~ 1,
  elevation_ws    = meanDailyMax ~ elevws,
  #elevation_cat   = meanDailyMax ~ elevcat,
  precipitation   = meanDailyMax ~ precip9120ws,
  air_temperature = meanDailyMax ~ tmean9120ws,
  elev_precip     = meanDailyMax ~ elevws + precip9120ws,
  elev_airtemp    = meanDailyMax ~ elevws + tmean9120ws,
  precip_airtemp  = meanDailyMax ~ precip9120ws + tmean9120ws,
  three_predictor = meanDailyMax ~ elevws + precip9120ws + tmean9120ws
)

# use a common covariance structure for fixed-effect comparisons; ML is required for likelihood-based comparison of different formulas
fixed_models = lapply(fixed_formulas, function(f) {
  ssn_lm(
    formula       = f,
    ssn.object    = sfsr_ssn,
    tailup_type   = "exponential",
    taildown_type = "exponential",
    euclid_type   = "gaussian",
    additive      = "afvArea",
    estmethod     = "ml"
  )
})

# compare model fit
fixed_results = purrr::imap_dfr(
  fixed_models,
  \(mod, nm) {
    glance(mod) |>
      mutate(model = nm, .before = 1)
  }
) |>
  arrange(AICc)

print(fixed_results)

# select lowest-AICc fixed-effect structure
best_fixed_name = fixed_results$model[1]
best_formula    = fixed_formulas[[best_fixed_name]]

cat("Selected fixed-effects model:", best_fixed_name, "\n")
print(best_formula)

# --------------------------------
# DEFINE CANDIDATE COVARIANCE STRUCTURES

# compare covariance structures while holding fixed effects constant
cov_structures = tribble(
  ~model,            ~tailup,       ~taildown,     ~euclid,
  "tailup",          "exponential", "none",        "none",
  "taildown",        "none",        "exponential", "none",
  "tailup_taildown", "exponential", "exponential", "none",
  "full",            "exponential", "exponential", "gaussian"
)

# ------------------------------------
# COMPARE COVARIANCE STRUCTURES (REML)

# fit candidate models with identical fixed effects and observations
cov_models = lapply(seq_len(nrow(cov_structures)), function(i) {
  ssn_lm(
    formula       = best_formula,
    ssn.object    = sfsr_ssn,
    tailup_type   = cov_structures$tailup[i],
    taildown_type = cov_structures$taildown[i],
    euclid_type   = cov_structures$euclid[i],
    additive      = "afvArea",
    estmethod     = "reml"
  )
})

names(cov_models) = cov_structures$model

# compare model fit using AICc
cov_results = purrr::imap_dfr(
  cov_models,
  \(mod, nm) glance(mod) |> mutate(model = nm, .before = 1)
) |>
  arrange(AICc)

print(cov_results)

# identify best covariance structure
best_cov_name = cov_results$model[1]

best_cov = cov_structures |>
  filter(model == best_cov_name)

print(best_cov)

# --------------------------
# LEAVE-ONE-OUT CROSS VALIDATION

# evaluate predictive performance for covariance candidates
cv_results = purrr::imap_dfr(
  cov_models,
  \(mod, nm) {
    loocv(mod) |>
      mutate(model = nm, .before = 1)
  }
) |>
  arrange(RMSPE)

print(cv_results)

# compare AICc and cross-validation rankings
model_comparison = cov_results |>
  select(model, AIC, AICc, BIC) |>
  left_join(
    cv_results |>
      select(model, RMSPE, bias, cor2),
    by = "model"
  ) |>
  arrange(AICc)

print(model_comparison)

# --------------------------------
# FIT FINAL SELECTED MODEL (REML)

# retrieve selected model
final_mod = cov_models[[best_cov_name]]

# examine model results
summary(final_mod)
coef(final_mod)
varcomp(final_mod)
tidy(final_mod, conf.int = TRUE)
glance(final_mod)

# -----------------
# MODEL DIAGNOSTICS

# augment observed sites with fitted values and residuals
aug_obs = augment(final_mod)

head(aug_obs)
names(aug_obs)

# fitted versus residuals
ggplot(aug_obs, aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Model Residuals",
    x = "Fitted Temperature (°C)",
    y = "Residual (°C)"
  )

# standard SSN2 diagnostic plots
plot(final_mod, which = 1:6)

# map model residuals
ggplot() +
  geom_sf(data = sfsr_ssn$edges, color = "grey70", linewidth = 0.3) +
  geom_sf(data = aug_obs, aes(color = .resid), size = 1.5) +
  scale_color_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0
  ) +
  theme_bw() +
  labs(
    title = "Spatial Distribution of Model Residuals",
    color = "Residual (°C)"
  )

# --------------------------
# PREDICTION (KRIGING)

# predict separately at both prediction datasets
# include standard errors and 95% prediction intervals
pred_results = lapply(c("pred_nsi_sfsr", "pred_1km_sfsr"), function(nm) {
  augment(
    final_mod,
    newdata = nm,
    se_fit = TRUE,
    interval = "prediction",
    level = 0.95
  )
})

names(pred_results) = c("pred_nsi_sfsr", "pred_1km_sfsr")

# inspect predictions
lapply(pred_results, head)
lapply(pred_results, names)

# --------------------------------
# MAP PREDICTED TEMPERATURES

# original NSI prediction points
ggplot() +
  geom_sf(data = sfsr_ssn$edges, color = "grey70", linewidth = 0.3) +
  geom_sf(data = pred_results$pred_nsi_sfsr, aes(color = .fitted), size = 0.8) +
  scale_color_viridis_c(option = "H") +
  theme_bw() +
  labs(
    title = "Predicted August Temperature — NSI Points",
    color = "Temperature (°C)"
  )

# custom 1-km prediction points
ggplot() +
  geom_sf(data = sfsr_ssn$edges, color = "grey70", linewidth = 0.3) +
  geom_sf(data = pred_results$pred_1km_sfsr, aes(color = .fitted), size = 0.8) +
  scale_color_viridis_c(option = "H") +
  theme_bw() +
  labs(
    title = "Predicted August Temperature — 1-km Points",
    color = "Temperature (°C)"
  )

# --------------------------------
# MAP PREDICTION UNCERTAINTY

# prediction standard errors at 1-km points
ggplot() +
  geom_sf(data = sfsr_ssn$edges, color = "grey70", linewidth = 0.3) +
  geom_sf(data = pred_results$pred_1km_sfsr, aes(color = .se.fit), size = 0.8) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(
    title = "Prediction Standard Errors — 1-km Points",
    color = "SE (°C)"
  )

# --------------------------------
# SAVE MODEL OUTPUTS

# create output directory
model_dir = here("output", "sfsr_examp", "models")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

# save final model and comparison tables
saveRDS(final_mod, file.path(model_dir, "sfsr_final_model.rds"))
write_csv(fixed_results, file.path(model_dir, "fixed_effects_comparison.csv"))
write_csv(model_comparison, file.path(model_dir, "covariance_comparison.csv"))

# save observation diagnostics
st_write(
  aug_obs,
  file.path(model_dir, "sfsr_model_diagnostics.gpkg"),
  delete_dsn = TRUE
)

# save both prediction datasets
for (nm in names(pred_results)) {
  st_write(
    pred_results[[nm]],
    file.path(model_dir, paste0(nm, "_predictions.gpkg")),
    delete_dsn = TRUE
  )
}

### END SCRIPT
