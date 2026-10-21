
# ATTACH STREAMCAT COVARIATES TO NSI
#
# Author(s): Mike Ackerman
#
# Purpose: Explore, retrieve, and attach StreamCat landscape covariates
#   to South Fork Salmon River NSI hydrography.
#
# Source: https://usepa.github.io/StreamCatTools/
#
# Created: October 15, 2026
#   Last modified: October 20, 2026
#
# NOTES: Explore StreamCat metrics by category, source dataset, spatial
#   support, and temporal coverage. Retrieve four candidate predictors
#   for August stream temperature and evaluate their distributions,
#   missing values, spatial patterns, and correlations.
#
#   Selected covariates:
#     elevws          - mean watershed elevation
#     tmean9120ws     - mean watershed air temperature (1991-2020)
#     precip9120ws    - mean watershed precipitation (1991-2020)
#     pctconif2019cat - evergreen forest cover in local catchment (2019)
#
#   StreamCat API access is required for metadata and data retrieval.

# clear environment
rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(sf)
library(StreamCatTools)

# --------------
# LOAD NSI DATA
nsi_flow_sfsr = st_read(here("data", "NSI", "processed", "nsi_flow_sfsr.gpkg"))

# ------------------------------------
# EXPLORE AVAILABLE STREAMCAT METRICS

# retrieve full StreamCat variable catalog
sc_vars = sc_get_params(param = "variable_info")

# inspect metadata structure
glimpse(sc_vars)
nrow(sc_vars)


sc_vars |> count(category, sort = TRUE)          # count metrics by category
sc_vars |> count(category, dataset, sort = TRUE) # count metrics by source dataset
sc_vars |> count(aoi, sort = TRUE)               # examine available spatial supports
sc_vars |> count(year, sort = TRUE)              # examine temporal coverage

# organize metadata into a searchable catalog
sc_catalog = sc_vars |>
  select(
    category,
    metric,
    short_description,
    long_description,
    units,
    aoi,
    year,
    dataset,
    source_name,
    source_URL
  ) |>
  arrange(category, metric)

# write metric catalog to file
write_csv(sc_catalog, here("data", "metadata", "streamcat_metric_catalog.csv"))

# search for variables potentially relevant to stream temperature
sc_candidates = sc_catalog |>
  filter(
    str_detect(
      str_to_lower(
        paste(
          category,
          metric,
          short_description,
          long_description,
          dataset
        )
      ),
      "elev|precip|temp|forest|canopy|shade|solar|snow|baseflow|riparian"
    )
  )

# inspect climate-related metrics
sc_catalog |>
  filter(str_detect(str_to_lower(category), "climate")) |>
  select(metric, short_description, units, aoi, year, dataset)

# --------------------------------
# IDENTIFY VALID STREAMCAT METRICS

# retrieve valid metric names
sc_metrics = sc_get_params(param = "metric_names") |>
  unlist() |>
  as.character()

# search valid names for potential predictors
grep(
  "elev|tmean9120|precip9120|pctconif2019|bfi",
  sc_metrics,
  value = TRUE,
  ignore.case = TRUE
)

# --------------------------------
# SELECT STREAMCAT COVARIATES

# selected metrics for August stream temperature modeling
selected_metrics = c("elev", "tmean9120", "precip9120")

# verify selected metrics are available (should be TRUE)
all(selected_metrics %in% sc_metrics)

# inspect metadata for selected metrics
selected_info = sc_catalog |>
  filter(
    metric %in% selected_metrics |
      str_detect(metric, c(selected_metrics))
  )

# view
print(selected_info, width = Inf)

# --------------------------------
# RETRIEVE STREAMCAT DATA

# identify unique NSI COMIDs
nsi_comids = nsi_flow_sfsr |>
  st_drop_geometry() |>
  filter(!is.na(COMID)) |>
  distinct(COMID) |>
  pull(COMID) |>
  as.character()

# test retrieval using five COMIDs
sc_test = sc_get_data(
  comid  = head(nsi_comids, 5),
  metric = paste(selected_metrics, collapse = ","),
  aoi    = "cat,ws"
)

# inspect test results
glimpse(sc_test)
names(sc_test)

# retrieve selected metrics for all SFSR COMIDs
# NOTE: If the API fails, consider smaller batches.
sc_data = sc_get_data(
  comid  = nsi_comids,
  metric = paste(selected_metrics, collapse = ","),
  aoi    = "cat,ws"
)

# inspect returned streamcat covariate data
glimpse(sc_data)
names(sc_data)
nrow(sc_data)

# --------------------------------
# QUALITY CONTROL OF STREAMCAT DATA

# standardize COMID type
sc_data = sc_data |>
  rename(COMID = comid) |>
  mutate(COMID = as.character(COMID))

# identify returned covariates
sc_covars = setdiff(names(sc_data), "COMID")

# verify unique COMIDs
stopifnot(
  !anyNA(sc_data$COMID),
  !anyDuplicated(sc_data$COMID)
)

# verify selected covariates were returned
stopifnot(
  all(model_covars %in% names(sc_data))
)

# identify NSI COMIDs without matching StreamCat records
missing_comids = nsi_flow_sfsr |>
  st_drop_geometry() |>
  mutate(COMID = as.character(COMID)) |>
  distinct(COMID) |>
  anti_join(sc_data, by = "COMID")

# inspect unmatched COMIDs
head(missing_comids)

# summarize selected covariates
sc_summary = sc_data |>
  select(all_of(sc_covars)) |>
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) |>
  group_by(variable) |>
  summarise(
    n = n(),
    n_missing = sum(is.na(value)),
    pct_missing = 100 * mean(is.na(value)),
    minimum = if (all(is.na(value))) NA_real_ else min(value, na.rm = TRUE),
    median = if (all(is.na(value))) NA_real_ else median(value, na.rm = TRUE),
    mean = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE),
    maximum = if (all(is.na(value))) NA_real_ else max(value, na.rm = TRUE),
    .groups = "drop"
  )

print(sc_summary)

# --------------------------------
# ATTACH STREAMCAT COVARIATES TO NSI

# join StreamCat metrics to processed NSI flowlines
nsi_flow_sfsr_sc = nsi_flow_sfsr |>
  mutate(COMID = as.character(COMID)) |>
  left_join(
    sc_data,
    by = "COMID",
    relationship = "many-to-one"
  )

# verify the join preserved all flowlines
nrow(nsi_flow_sfsr_sc) == nrow(nsi_flow_sfsr)

# inspect covariate coverage
nsi_flow_sfsr_sc |>
  st_drop_geometry() |>
  summarise(
    n_segments = n(),
    across(
      all_of(sc_covars),
      ~ sum(!is.na(.x)),
      .names = "n_{.col}"
    )
  )

# --------------------------------
# EXPLORE COVARIATE DISTRIBUTIONS

# reshape selected covariates for plotting
sc_long = nsi_flow_sfsr_sc |>
  st_drop_geometry() |>
  select(all_of(sc_covars)) |>
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  )

# plot distributions
ggplot(sc_long, aes(x = value)) +
  geom_histogram(bins = 30, color = "white") +
  facet_wrap(~ variable, scales = "free") +
  theme_bw() +
  labs(
    title = "StreamCat Covariate Distributions — SFSR",
    x = "Covariate Value",
    y = "Number of Stream Segments"
  )

# --------------------------------
# MAP STREAMCAT COVARIATES

# map watershed elevation
ggplot() +
  geom_sf(data = nsi_flow_sfsr_sc, aes(color = elevws), linewidth = 0.7) +
  scale_color_viridis_c(na.value = "grey80") +
  theme_bw() +
  labs(
    title = "Mean Watershed Elevation (ws) — SFSR",
    color = "Elevation"
  )
ggplot() +
  geom_sf(data = nsi_flow_sfsr_sc, aes(color = elevcat), linewidth = 0.7) +
  scale_color_viridis_c(na.value = "grey80") +
  theme_bw() +
  labs(
    title = "Mean Watershed Elevation (cat) — SFSR",
    color = "Elevation"
  )

# map watershed air temperature
ggplot() +
  geom_sf(data = nsi_flow_sfsr_sc, aes(color = tmean9120ws), linewidth = 0.7) +
  scale_color_viridis_c(na.value = "grey80") +
  theme_bw() +
  labs(
    title = "Mean Watershed Air Temperature (ws) — SFSR",
    color = "Air temp."
  )
ggplot() +
  geom_sf(data = nsi_flow_sfsr_sc, aes(color = tmean9120cat), linewidth = 0.7) +
  scale_color_viridis_c(na.value = "grey80") +
  theme_bw() +
  labs(
    title = "Mean Watershed Air Temperature (cat) — SFSR",
    color = "Air temp."
  )

# map watershed precipitation
ggplot() +
  geom_sf(data = nsi_flow_sfsr_sc, aes(color = precip9120ws), linewidth = 0.7) +
  scale_color_viridis_c(na.value = "grey80") +
  theme_bw() +
  labs(
    title = "Mean Watershed Precipitation (ws) — SFSR",
    color = "Precip."
  )
ggplot() +
  geom_sf(data = nsi_flow_sfsr_sc, aes(color = precip9120cat), linewidth = 0.7) +
  scale_color_viridis_c(na.value = "grey80") +
  theme_bw() +
  labs(
    title = "Mean Watershed Precipitation (cat) — SFSR",
    color = "Precip."
  )

# --------------------------------
# EXAMINE COVARIATE CORRELATIONS

# calculate Pearson correlations among selected covariates
cor_mat = nsi_flow_sfsr_sc |>
  st_drop_geometry() |>
  select(all_of(sc_covars)) |>
  cor(
    use = "pairwise.complete.obs",
    method = "pearson"
  )

# inspect correlation matrix
round(cor_mat, 2)

# reshape correlation matrix for plotting
cor_df = as.data.frame(as.table(cor_mat)) |>
  rename(
    variable_1 = Var1,
    variable_2 = Var2,
    correlation = Freq
  )

# plot correlation matrix
ggplot(
  cor_df,
  aes(
    x = variable_1,
    y = variable_2,
    fill = correlation
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = round(correlation, 2)),
    size = 3
  ) +
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
    title = "StreamCat Covariate Correlations — SFSR",
    x = NULL,
    y = NULL,
    fill = "Pearson r"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# --------------------------------
# SAVE STREAMCAT-ENRICHED HYDROGRAPHY

st_write(
  nsi_flow_sfsr_sc,
  here("data", "NSI", "processed", "nsi_flow_sfsr_streamcat.gpkg"),
  delete_dsn = TRUE
)

### END SCRIPT
