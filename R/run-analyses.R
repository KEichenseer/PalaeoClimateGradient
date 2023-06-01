# Title: Bayesian multi-proxy reconstruction of early Eocene latitudinal temperature gradients
# Authors: Kilian Eichenseer & Lewis A. Jones
# Updated: 2022-08-11
#============================= RUN ANALYSES ==================================#
#
# check for missing packages and install if needed
list.of.packages <- c("cli", "doParallel","foreach","parallel","mvnfast","truncnorm",
                      "broom", "ggplot2", "ggtext", "ggthemes", "readxl", "mapproj", "magick",
                      "png", "ggpubr", "mcmcse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>=1) install.packages(new.packages)

# rlang needs to be version >= 1.1.0!
if (packageVersion("rlang") < "1.1.0") {
  # If the version is older than 1.1.0, update and reload
  install.packages("rlang")
}

# vctrs needs to be version >= 0.5.0 !
if (packageVersion("vctrs") < "0.5.0") {
  # If the version is older than 0.5.0 , update and reload
  install.packages("vctrs")
}

# ggplot2 needs to be version >= 3.4.0!
if (packageVersion("ggplot2") < "3.4.0") {
  # If the version is older than 3.4.0, update and reload
  install.packages("ggplot2")
}

# Data preparation -------------------------------------------------------------
source("./R/data_preparation/01_data_wrangle_hollis2019.R")
rm(list = ls())
source("./R/data_preparation/02_prepare_hollis2019.R")
rm(list = ls())
source("./R/data_preparation/03_palaeorotate_data.R")
rm(list = ls())
#
#
## Figures
#
# Fig. 1 - Geographical distribution of EECO proxy data and palaeotemperatures
source("./R/figures/fig_1.R", print.eval=TRUE)
#
# for Fig. 2 - Model reconstructions of simulated latitudinal temperature gradients
source("./R/figures/fig_2.R", print.eval=TRUE)
#
# Fig. 3 - Modern latitudinal temperature gradients - takes a minute
source("./R/figures/fig_3.R", print.eval=TRUE)
#
# Fig. 4 - Eocene latitudinal temperature gradient
source("./R/figures/fig_4.R", print.eval=TRUE)
#
#
## Analyses
#
# Model runs (careful, takes very long!) 
#
# Estimate modern latitudinal temperature gradient - run model
# !!! WARNING !!! This can take 30 minutes or longer!
source("./R/analyses/run_modern_climate_model.R")
#
# !!! WARNING !!! slow
# Estimate EECO latitudinal temperature gradient - run model
source("./R/analyses/run_EECO_climate_model.R")
#
# !!! WARNING !!! takes 10+ hours
# Estimate simulated latitudinal temperature gradients - run model
source("./R/analyses/run_simulated_gradients.R")
#
#
# processing model results
#
# global mean temperature
source("./R/analyses/processing/global_mean_temperatures.R")
# gradients and temperature differences
source("./R/analyses/processing/calculate_temperature_gradients_from_models.R")
# sdy
source("./R/analyses/processing/assess_sdy.R")
# R2 for simulated gradients
source("./R/analyses/processing/simulated_gradients_R_squared.R")

