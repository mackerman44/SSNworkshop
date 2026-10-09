##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# R script for SSNbler_vignette.pdf 
#
# Includes additional comments and R code not included in the
# SSNbler vignette
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# ------- Install SSNbler (if necessary) and load library ------

## Install package from CRAN repository 
## install.packages("SSNbler")

## Install package from Github
##install.packages("remotes") ## Install remotes if needed
##library(remotes)

## Install the latest version of SSNbler from github
##remotes::install_github("pet221/SSNbler", ref = "main")

## Load SSNbler library
library(SSNbler)

# ------------- Import and view the input data ----------------

## Copy the example dataset to a temporary directory
copy_streams_to_temp()
path <- paste0(tempdir(), "/streamsdata")

## List local input files
list.files(path)

## Load the sf package and import the streams, observation sites, and
## two prediction datasets
library(sf)
MF_streams <- st_read(paste0(path, "/MF_streams.gpkg"))
MF_obs <- st_read(paste0(path, "/MF_obs.gpkg"))
MF_pred1km <- st_read(paste0(path, "/MF_pred1km.gpkg"))
MF_CapeHorn <- st_read(paste0(path, "/MF_CapeHorn.gpkg"))

## Notice that the imported data are of class sf data.frame
class(MF_obs)

## Look at column names and their definitions
names(MF_obs)
help(MF_obs)

## Plot the data using ggplot2
library(ggplot2)
ggplot() +
  geom_sf(data = MF_streams) +
  geom_sf(data = MF_CapeHorn, color = "gold", size = 1.7) +
  geom_sf(data = MF_pred1km, colour = "purple", size = 1.7) +
  geom_sf(data = MF_obs, color = "blue", size = 2) +
  coord_sf(datum = st_crs(MF_streams))

# --------------- Build the LSN -------------------------------

## Set the lsn.path variable
lsn.path <- paste0(tempdir(), "/mf04")

## Check the minimum line feature length
## snap_tolerance must be < this value
min(st_length(MF_streams))

## Build the LSN and take note of the messages printed to the console
edges <- lines_to_lsn(
  streams = MF_streams,  ## sf object with LINESTRING geometry 
  lsn_path = lsn.path,   ## path to output directory
  check_topology = TRUE, ## should topology be checked?
  snap_tolerance = 0.05, ## snap endnodes within this distance
  topo_tolerance = 20,   ## flag topological errors within this distance
  overwrite = TRUE       ## overwrite existing files?
)

## list new files saved in lsn.path
list.files(lsn.path)

# ----------- Incorporate sites into LSN -----------------------

## Incorporate observations. Pay attention to output messages in R
## console to ensure all sites are snapped successfully.
obs <- sites_to_lsn(
  sites = MF_obs,        ## sf object with POINT geometry
  edges = edges,         ## sf object output from lines_to_lsn()
  lsn_path = lsn.path,
  file_name = "obs",     ## base filename for output sites
  snap_tolerance = 100,  ## snap sites within this distance of edge
  save_local = TRUE,     ## should outputs be saved locally?
  overwrite = TRUE
)

## Predictions at 1km intervals
preds <- sites_to_lsn(
  sites = MF_pred1km,
  edges = edges,
  save_local = TRUE,
  lsn_path = lsn.path,
  file_name = "pred1km.gpkg",
  snap_tolerance = 100,
  overwrite = TRUE
)

## Prediction sites in Cape Horn Creek
capehorn <- sites_to_lsn(
  sites = MF_CapeHorn,
  edges = edges,
  save_local = TRUE,
  lsn_path = lsn.path,
  file_name = "CapeHorn.gpkg",
  snap_tolerance = 100,
  overwrite = TRUE
)

## list files saved in lsn.path. Notice that geopackages 
## for the sites have been added locally
list.files(lsn.path)

## Notice the new columns rid, ratio, and snapdist
names(obs)

## Summarise snapdist, looking for sites that were snapped large
## distances
summary(obs$snapdist)

## Look at ratio values to ensure the range is 0<= ratio <= 1
summary(obs$ratio)

# --------- Calculate upstream distance for edges ------------

edges <- updist_edges(
  edges = edges,         ## sf object output from lines_to_lsn()
  save_local = TRUE,
  lsn_path = lsn.path,
  calc_length = TRUE     ## should a new Length column be 
)                        ## added to edges?


## View edges column names
names(edges)

## Summarise upDist
summary(edges$upDist)

# ---------- Calculate upstream distance for sites -----------

site.list <- updist_sites(
  sites = list(          ## A named list of sf objects created in 
    obs = obs,           ## sites_to_lsn()   
    pred1km = preds,
    CapeHorn = capehorn
  ),
  edges = edges,         ## sf object output from lines_to_lsn()
  length_col = "Length", ## column in edges containing length
  save_local = TRUE,
  lsn_path = lsn.path
)

## site.list is a named list of sf objects, with length = 3
length(site.list)

## View output site.list names
names(site.list)

## sf objects can be accessed within the list using the list names
## For example, get the class, dimensions, and column names of obs
class(site.list$obs)
dim(site.list$obs)  
names(site.list$obs)

## Plot the upstream distances for edges and obs
ggplot() +
  geom_sf(data = edges, aes(color = upDist)) +
  geom_sf(data = site.list$obs, aes(color = upDist)) +
  coord_sf(datum = st_crs(MF_streams)) +
  scale_color_viridis_c()

# ------- Calculate AFV for edges --------------------------------

## Summarize h2oAreaKm2 and check for zeros
summary(edges$h2oAreaKm2) 

edges <- afv_edges(
  edges = edges,
  infl_col = "h2oAreaKm2",   ## edges column used to calc segment PI
  segpi_col = "areaPI",      ## name of new edges segment PI column
  afv_col = "afvArea",       ## name of new edges AFV column
  lsn_path = lsn.path
)

## Look at edges column names. Notice the two new columns, 
## areaPI and afvArea
names(edges)

## Summarize the AFV column to ensure it ranges from 0 <= AFV <= 1
summary(edges$afvArea)

# ------- Calculate AFV for sites -------------------------------

site.list <- afv_sites(
  sites = site.list,    ## named list of sf objects
  edges = edges,        ## sf object output from lines_to_lsn()
  afv_col = "afvArea",  ## name of edges column containing AFVs
  save_local = TRUE,   
  lsn_path = lsn.path
)

## Output is the same named list of sf objects as input sites,
## with a new column added
names(site.list)

## View column names in pred1km. Notice the new afvArea column
names(site.list$pred1km)

## Summarise AFVs in pred1km to ensure they range from 0 <= AFV <= 1
## Remember, large numbers of zeros can be problematic
summary(site.list$pred1km$afvArea)

# ---------- Assemble the SSN Object -----------------------------

mf04_ssn <- ssn_assemble(
  edges = edges,                      ## sf object
  lsn_path = lsn.path,
  obs_sites = site.list$obs,          ## sf object of observations
  preds_list = site.list[c("pred1km", ## named list of sf objects
                           "CapeHorn")],
  ssn_path = paste0(path, "/MiddleFork04.ssn"), ## output directory
  import = TRUE,                      ## should an SSN object be returned?
  check = TRUE,                       ## check validity of SSN object?
  afv_col = c("afvArea"),             ## names of AFV columns to check
  overwrite = TRUE
)

## Look at files in MiddleFork04.ssn
list.files(paste0(path, "/MiddleFork04.ssn"))

## Check class for mf04_ssn
class(mf04_ssn)

## SSN object is a named list
is.list(mf04_ssn)

## Get SSN object (list) names
names(mf04_ssn)

## Print path to local .ssn
mf04_ssn$path

## Print names of prediction datasets in SSN object
names(mf04_ssn$preds)

## Look at the observations in the SSN. Notice the new columns
## added when the SSN was assembled (pid, locID, netID, netgeom)
class(mf04_ssn$obs)
names(mf04_ssn$obs)
View(mf04_ssn$obs)

## Plot mf04_ssn
ggplot() +
  geom_sf(
    data = mf04_ssn$edges,       ## plot edges sf object
    color = "medium blue",       ## in blue
    aes(linewidth = h2oAreaKm2)  ## linewidth proportional to h2oAreaKm2
  ) +
  scale_linewidth(range = c(0.1, 2.5)) + ## set linewidth aesthetics
  geom_sf(
    data = mf04_ssn$preds$pred1km, ## plot pred1km 
    size = 1.5,                    ## set point size
    shape = 21,                    ## set point style
    fill = "white",                ## set point fill colour
    color = "dark grey"            ## set point outline colour
  ) +
  geom_sf(
    data = mf04_ssn$obs,           ## plot observation points
    size = 1.7,                    ## set point size
    aes(color = Summer_mn)         ## colour by column name
  ) +
  coord_sf(datum = st_crs(MF_streams)) +  ## set datum for plot graticules
  scale_color_viridis_c() +               ## use viridis colour palette
  labs(color = "Temperature", linewidth = "WS Area") + ## labels for legend
  theme(
    legend.text = element_text(size = 8),  ## legend text size
    legend.title = element_text(size = 10) ## legend title text size
  )

# ---- Create Distance Matrices ----------------------------------
## load SSN2 library
library(SSN2)

## Generate hydrologic distance matrices for observations
## and prediction set pred1km
ssn_create_distmat(mf04_ssn, predpts = "pred1km")

## Look at distance matrix files created for obs
list.files(paste0(mf04_ssn$path, "/distance/obs"))

## Look at distance matrix files created for pred1km
list.files(paste0(mf04_ssn$path, "/distance/pred1km"))

## Get a named list of distance matrices for obs
obs.distmat<- ssn_get_stream_distmat(mf04_ssn,
                                     name = "obs")
names(obs.distmat)          ## look at list element names
colnames(obs.distmat[[2]])  ## column names for dist.net2
View(obs.distmat[[2]])      ## View dist.net2

## Create symmetric hydrologic distance matrix
obs.distmat2 <- obs.distmat[[2]] + t(obs.distmat[[2]])
View(obs.distmat2)

# ---- Fit spatial stream-network model --------------------------

## Fit a spatial linear model to Summer mean temperature with a
## mixture of TU/TD/EUC covariance models
ssn_mod <- ssn_lm(
  formula = Summer_mn ~ ELEV_DEM + AREAWTMAP,
  ssn.object = mf04_ssn,
  tailup_type = "exponential",
  taildown_type = "spherical",
  euclid_type = "gaussian",
  additive = "afvArea"
)

## Summarise model results
summary(ssn_mod)

## Get class for model output
class(ssn_mod)

## Notice that the ssn.object is stored in ssn_mod
## for reproducibility
summary(ssn_mod$ssn.object)

