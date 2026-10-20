# PREPARE NORWEST OBSERVATIONS
#
# Author(s): Mike Ackerman
#
# Purpose: Prepare NorWeST observed stream temperature data for the South Fork Salmon River.
#
# Created: October 15, 2026
#   Last modified: October 20, 2026
#
# NOTES: I'm going to summarize summer (e.g., August) stream temperatures by deployment (site-year) from 
#   legacy NorWeST data, join the summaries to observation locations, and trim to the South Fork Salmon River.

# clear environment
rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(sf)
library(readxl)
library(lubridate)

# ----------------------
# LOAD SPATIAL REFERENCE

# load processed NSI hydrography and study area
nsi_flow_sfsr = st_read(here("data", "NSI", "processed", "nsi_flow_sfsr.gpkg"))
study_area    = st_read(here("data", "study_areas", "study_area_sfsr.gpkg"))

# confirm the coordinate reference systems match
stopifnot(st_crs(study_area) == st_crs(nsi_flow_sfsr))

# -----------------------------
# PREP NORWEST OBSERVATION DATA
norwest_dir = "C:/Users/mikea/nezperce365/MFO Research - General/Riverscapes/Data/Temp/NorWeSTv2"

# load daily temperature summaries
norwest_daily = read_xlsx(
  path = file.path(
    norwest_dir,
    "NorWeST_v2_ObservedStreamTempDailySummaries_Salmon_AllDays",
    "NorWeST_v2_ObservedStreamTempDailySummaries_Salmon_AllDays.xlsx"
  ),
  sheet = "DailyMetricOutputForWeb"
)

# summarize August temperatures by deployment (site-year)
obs_df = norwest_daily |>
  mutate(
    year = year(SampleDate),
    across(c(DailyMean, DailyMax), ~ na_if(.x, -9999))
  ) |>
  # here, we'll just use August temperature data
  filter(month(SampleDate) == 8) |>
  group_by(OBSPRED_ID, year) |>
  summarise(
    nDaysMean     = sum(!is.na(DailyMean)),
    nDaysMax      = sum(!is.na(DailyMax)),
    meanDailyMean = mean(DailyMean, na.rm = TRUE),
    meanDailyMax  = mean(DailyMax, na.rm = TRUE),
    maxDailyMax   = max(DailyMax, na.rm = TRUE),
    nMean18       = sum(DailyMean > 18, na.rm = TRUE),
    nMax18        = sum(DailyMax > 18, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    across(
      c(meanDailyMean, meanDailyMax, maxDailyMax),
      ~ ifelse(is.nan(.x) | is.infinite(.x), NA_real_, .x)
    )
  )

# NOTE: OBSPRED_ID should be unique to a site-year deployment, not a site or locations

# verify whether OBSPRED_ID uniquely identifies site-year deployments
obs_df |>
  count(OBSPRED_ID) |>
  filter(n > 1)

# load observation locations
obs_pts = st_read(
  file.path(
    norwest_dir,
    "NorWeSTv2_Salmon_ObservedTempPoints",
    "NorWeSTv2_Salmon_ObservedTempPoints.shp"
  )
) |>
  # transform CRS to share same CRS as hydrography
  st_transform(st_crs(nsi_flow_sfsr))

# inspect
ggplot() +
  geom_sf(data = study_area, fill = NA, linewidth = 0.6) +
  geom_sf(data = obs_pts, size = 0.5, color = "red") +
  theme_bw()

# -----------------------------------
# SOME QC OF NORWEST OBSERVATION DATA

# identify site-year deployments WITHOUT August daily temperature data
unmatched_obs = obs_pts |>
  st_filter(study_area) |>
  # returns site-year deployments where OBSPRED_ID isn't in obs_df
  anti_join(obs_df, by = "OBSPRED_ID")

# by sample year, site-year deployments WITHOUT August daily temperature data
unmatched_obs |>
  st_drop_geometry() |>
  count(SampleYear) |>
  arrange(SampleYear)

# check whether unmatched observations have daily records in other months
unmatched_check = unmatched_obs |>
  st_drop_geometry() |>
  select(OBSPRED_ID, SampleYear) |>
  left_join(
    norwest_daily |>
      mutate(month = lubridate::month(SampleDate)) |>
      group_by(OBSPRED_ID) |>
      summarise(
        n_daily = n(),
        months = paste(sort(unique(month)), collapse = ", "),
        .groups = "drop"
      ),
    by = "OBSPRED_ID"
  )

# inspect results
unmatched_check |>
  count(n_daily, months, sort = TRUE)

# ----------------------------------------------
# SPATIAL JOIN DEPLOYMENT SUMMARIES TO LOCATIONS

# join temperature summary to spatial site-year observation locations
obs_sfsr = obs_pts |>
  select(OBSPRED_ID, PERMA_FID, COMID, GNIS_NAME, SampleYear) |>
  # inner_join excludes observations without a OBSPRED_ID in obs_df
  inner_join(
    obs_df, 
    by = "OBSPRED_ID",
    relationship = "one-to-one"
  ) |>
  # trim to study area
  st_filter(study_area)

# inspect the join
ggplot() +
  geom_sf(data = nsi_flow_sfsr, linewidth = 0.3, color = "steelblue") +
  geom_sf(data = obs_sfsr, size = 2, color = "red") +
  geom_sf(data = study_area, fill = NA, linewidth = 0.6) +
  theme_bw()

# more inspection
nrow(obs_sfsr)    # unique site-year deployments
sum(!is.na(obs_sfsr$meanDailyMean))
sum(is.na(obs_sfsr$meanDailyMean))

# check deployment year against date-derived year
obs_sfsr |>
  st_drop_geometry() |>
  filter(!is.na(year), SampleYear != year) |>
  select(OBSPRED_ID, PERMA_FID, SampleYear, year)

# examine temperature metrics
summary(obs_sfsr$meanDailyMean)
summary(obs_sfsr$meanDailyMax)
summary(obs_sfsr$maxDailyMax)
hist(obs_sfsr$nMean18)
hist(obs_sfsr$nMax18)

# check unique observation IDs before saving
stopifnot(!anyDuplicated(obs_sfsr$OBSPRED_ID))

# plot observation locations
ggplot() +
  geom_sf(data = study_area, fill = NA) +
  geom_sf(data = nsi_flow_sfsr, color = "steelblue", linewidth = 0.3) +
  geom_sf(
    data = filter(obs_sfsr, !is.na(meanDailyMean)),
    aes(color = meanDailyMean),
    size = 1.5
  ) +
  scale_color_viridis_c(option = "turbo", na.value = "grey70") +
  theme_bw() +
  labs(
    title = "NorWeST August Stream Observations — SFSR",
    color = "Mean Daily Mean\nTemperature (°C)"
  )

# -----------------
# SAVE OBSERVATIONS

# retain observations with missing August metrics for now; filter when modeling
st_write(
  obs_sfsr,
  here("data", "norwest", "processed", "sfsr_norwest_obs.gpkg"),
  delete_dsn = TRUE
)

### END SCRIPT
