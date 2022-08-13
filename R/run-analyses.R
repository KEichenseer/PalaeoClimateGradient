# Title: XXX
# Authors: Kilian Eichenseer & Lewis A. Jones
# Updated: 2022-08-11
#============================= RUN ANALYSES ==================================#
#
# check for missing packages and install if needed
list.of.packages <- c("cli", "doParallel","foreach","parallel","mvnfast","truncnorm",
                      "broom", "ggplot2", "ggthemes", "readxl", "mapproj", "magick",
                      "png", "ggpubr", "mcmcse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>=1) install.packages(new.packages)

# Data preparation -------------------------------------------------------------
source("./R/data_preparation/01_data_wrangle_hollis2019.R") ### Not yet working - Check!
rm(list = ls())
source("./R/data_preparation/02_prepare_hollis2019.R")
rm(list = ls())
source("./R/data_preparation/03_palaeorotate_data.R")
rm(list = ls())
#
# Analyses
# Fig. 1 - Geographical distribution of EECO proxy data and palaeotemperatures
source("./R/figures/Fig_1.R", print.eval=TRUE)
#
# for Fig. 2 - plot priors
source("./R/figures/prior_plot.R")
#
# Estimate modern latitudinal temperature gradient - run model
# !!! WARNING !!! This can take 30 minutes or longer!
source("./R/analyses/run_modern_climate_model.R")
#
# Estimate modern latitudinal temperature gradient from samples drawn
# with the EECO palaeolatitudinal sample distribution - run model
source("./R/analyses/run_modern_climate_models_EECO_p_lat_samples.R")
#
# Fig. 3 - plot modern gradients
source("./R/figures/fig_3.R", print.eval=TRUE)
#
# Estimate EECO latitudinal temperature gradient
source("./R/analyses/run_EECO_climate_model.R")
