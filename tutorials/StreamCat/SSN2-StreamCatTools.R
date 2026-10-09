# SSN2 StreamCatTools R Script

# Load packages ####
library(SSN2)
library(StreamCatTools)
library(ggplot2)


# Read in bugs Data ####
path <- "bugs.ssn"
bugs <- ssn_import(
  path = path,
  predpts = "pred",
  overwrite = TRUE
)

## print a summary
summary(bugs)

## visualize stream network, prediction data, observed data
ggplot() +
  # geom_sf(data = bugs$edges, linewidth = 0.1, alpha = 0.3, color = "lightgrey") +
  geom_sf(data = bugs$preds$pred, color = "black", size = 0.9) +
  geom_sf(data = bugs$obs, aes(color = Rich), size = 2) +
  scale_color_viridis_c(option = "H", limits = c(0, 59)) +
  theme_bw(base_size = 16)

## create distance matrices needed for spatial statistical modeling
ssn_create_distmat(bugs, predpts = "pred", overwrite = TRUE)

#:::::::::::::::::::::::::::::::::::::
# Using StreamCatTools ####
#:::::::::::::::::::::::::::::::::::::

# Read in StreamCat Data ####

## temperature and percent evergreen forest cover
sc_vars <- c("tmean8110", "pctconif2008")

## access StreamCat data based on COMID for observed sites
sc_data <- sc_get_data(
  metric = sc_vars,
  aoi = "cat",
  comid = bugs$obs$comid
)

## Look at the first few rows
head(sc_data)

## Add StreamCat variables to the observed data
bugs$obs <- merge(bugs$obs, sc_data, by = "comid")

## visualize stream network, StreamCat variables

### Percent evergreen cover (2008) at catchment level
ggplot() +
  # geom_sf(data = bugs$edges, linewidth = 0.1, alpha = 0.3, color = "lightgrey") +
  geom_sf(data = bugs$obs, aes(color = pctconif2008cat), size = 2) +
  scale_color_viridis_c(option = "A") +
  theme_bw(base_size = 16)

### Mean annual temperature (1981-2010) at catchment level
ggplot() +
  # geom_sf(data = bugs$edges, linewidth = 0.1, alpha = 0.3, color = "lightgrey") +
  geom_sf(data = bugs$obs, aes(color = tmean8110cat), size = 2) +
  scale_color_viridis_c(option = "D") +
  theme_bw(base_size = 16)

## access StreamCat data based on COMID for prediction sites
sc_pred_data <- sc_get_data(
  metric = sc_vars,
  aoi = "cat",
  comid = bugs$preds$pred$comid
)

## Add StreamCat variables to the observed data
bugs$preds$pred <- merge(bugs$preds$pred, sc_pred_data, by = "comid")

#:::::::::::::::::::::::::::::::::::::
# Fit an SSN Model Using StreamCat Data ####
#:::::::::::::::::::::::::::::::::::::

## Fit ssn_lm model with tailup and taildown components
ssn_mod <- ssn_lm(
  formula = Rich ~ ELEV_DEM + tmean8110cat + pctconif2008cat,
  ssn.object = bugs,
  tailup_type = "exponential",
  taildown_type = "exponential",
  additive = "afvArea"
)

## Summarize the model
summary(ssn_mod)

## augment prediction data with predictions
aug_preds <- augment(ssn_mod, newdata = "pred")

## visualize stream network, predictions, observed data
ggplot() +
  # geom_sf(data = bugs$edges, linewidth = 0.1, alpha = 0.3, color = "lightgrey") +
  geom_sf(data = aug_preds, aes(color = .fitted), size = 0.9) +
  geom_sf(data = bugs$obs, aes(color = Rich), size = 2) +
  scale_color_viridis_c(option = "H", name = "Rich", limits = c(0, 59)) +
  theme_bw(base_size = 16)

## block Kriging for an estimate of the average in the region
predict(ssn_mod, newdata = "pred", block = TRUE, interval = "prediction")

#:::::::::::::::::::::::::::::::::::::
# Other StreamCatTools Features ####
#:::::::::::::::::::::::::::::::::::::

## get variable information
var_info <- sc_get_params(param = "variable_info")
head(var_info)

## get names of metrics (and other information)
sc_get_metric_names(category = "Climate", aoi = c("Cat", "Ws"))

## get COMIDs from coordinates
sct_comid <- sc_get_comid(bugs$obs)

### format is comma separated string, so separate into standard R vector
sct_comid <- unlist(strsplit(sct_comid, ","))

### create a data frame to compare against COMIDs in data (they should match)
dat_comid <- data.frame(
  bugs_comid = bugs$obs$comid,
  sct_comid = sct_comid
)
head(dat_comid)
