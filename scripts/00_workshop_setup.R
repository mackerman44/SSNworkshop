# SSN WORKSHOP SETUP
#
# Author(s): Mike Ackerman
#
# Purpose: Set up R packages for the Spatial Stream Network (SSN) workshop, October 14-16, 2026, Boise
#
# Created: October 2, 2026
#   Last modified:

# clear environment
rm(list = ls())

# ----------------
# INSTALL PACKAGES

# install remotes, if needed
# install.packages("remotes")

# SSN2 and SSNbler are available from CRAN
# install.packages("SSN2")
# install.packages("SSNbler")

# workshop organizers recommend installing the latest versions from GitHub ("develop" branches are also available)
# remotes::install_github("USEPA/SSN2", ref = "main")
# remotes::install_github("pet221/SSNbler", ref = "main")

# -----------------
# PACKAGE CITATIONS
citation("SSN2")
citation("SSNbler")

# ----------------
# PACKAGE VERSIONS
packageVersion("SSN2")
packageVersion("SSNbler")

# ----------------
# LOAD PACKAGES
library(SSN2)
library(SSNbler)

# -------------------
# SESSION INFORMATION
R.version.string
sessionInfo()

### END SCRIPT