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

# --------------
# LOAD NSI DATA
sfsr_nsi = st_read(here("data", "hydrography", "NSI", "processed", "nsi_flow_sfsr.gpkg"))

# --------------------------
# LOAD POPULATION BOUNDARIES
load(here("data", "steelhead_pops.rda"))
sth_pop = st_transform(sth_pop, st_crs(sfsr_nsi))
st_crs(sth_pop) == st_crs(sfsr_nsi)

# inspect
sort(unique(sth_pop$TRT_POPID))   # all available populations
pop_polys = sth_pop |>            # populations of interest
  filter(TRT_POPID %in% c("SFSEC-s", "SFMAI-s"))

# map
ggplot() +
  geom_sf(data = pop_polys, aes(fill = TRT_POPID)) +
  theme_bw()

# dissolve pop_polys into a single study area
study_area = pop_polys |>
  summarise()

rm(sth_pop, pop_polys)

# -----------------------------
# PREP NORWEST OBSERVATION DATA
nrwst_dir = "C:/Users/mikea/OneDrive - nezperce365/Workspace/norwest/"

# load daily temperature summaries
nrwst_daily = read_xlsx(
  path = file.path(
    nrwst_dir,
    "NorWeST_v2_ObservedStreamTempDailySummaries_Salmon_AllDays.xlsx"
  ),
  sheet = "DailyMetricOutputForWeb"
)

# summarize daily temperatures by deployment (site-year)
obs_df = nrwst_daily |>
  mutate(
    year = year(SampleDate),
    month = month(SampleDate),
    across(c(DailyMean, DailyMax), ~ na_if(.x, -9999))
  ) |>
  # let's just use August data
  filter(month == 8) |>
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

# check whether deployments are unique
obs_df |>
  count(OBSPRED_ID) |>
  filter(n > 1)

# load observation locations
salmon_obs_pts = st_read(
  file.path(
    nrwst_dir,
    "NorWeSTv2_Salmon_ObservedTempPoints",
    "NorWeSTv2_Salmon_ObservedTempPoints.shp"
  )
) |>
  st_transform(st_crs(sfsr_nsi))

# join temperature metrics to spatial observation locations
sfsr_obs = salmon_obs_pts |>
  select(OBSPRED_ID, PERMA_FID, COMID, GNIS_NAME, SampleYear) |>
  left_join(obs_df, by = "OBSPRED_ID") |>
  # trim to study area
  st_filter(study_area)

# inspect the join
nrow(sfsr_obs)
sum(!is.na(sfsr_obs$meanDailyMean))
sum(is.na(sfsr_obs$meanDailyMean))

# check deployment year against date-derived year
sfsr_obs |>
  st_drop_geometry() |>
  filter(!is.na(year), SampleYear != year) |>
  select(OBSPRED_ID, PERMA_FID, SampleYear, year)

# examine temperature metrics
summary(sfsr_obs$meanDailyMean)
summary(sfsr_obs$meanDailyMax)
summary(sfsr_obs$maxDailyMax)

# check unique observation IDs before saving
stopifnot(!anyDuplicated(sfsr_obs$OBSPRED_ID))

# plot observation locations
ggplot() +
  geom_sf(data = study_area, fill = NA) +
  geom_sf(data = sfsr_nsi, color = "steelblue", linewidth = 0.3) +
  geom_sf(
    data = filter(sfsr_obs, !is.na(meanDailyMean)),
    aes(color = meanDailyMean),
    size = 1.5
  ) +
  scale_color_viridis_c(na.value = "grey70") +
  theme_bw() +
  labs(color = "Mean Daily\nTemperature (°C)")

# write observations to file
# retain observations with missing August metrics for now; filter when modeling
st_write(
  sfsr_obs,
  here("data", "norwest", "processed", "sfsr_norwest_obs.gpkg"),
  delete_dsn = TRUE
)

### END SCRIPT
