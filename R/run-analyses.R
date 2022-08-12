# Title: XXX
# Authors: Kilian Eichenseer & Lewis A. Jones
# Updated: 2022-08-11
#============================= RUN ANALYSES ==================================#

# Data preparation -------------------------------------------------------------
source("./R/subscripts/data_preparation/01_data_wrangle_hollis2019.R")
rm(list = ls())
source("./R/subscripts/data_preparation/02_prepare_hollis2019.R")
rm(list = ls())
source("./R/subscripts/data_preparation/03_palaeorotate_data.R")
rm(list = ls())
