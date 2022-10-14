source("./R/functions/model_processing/temp_from_gradient.R")
source("./R/functions/model_processing/combine_posterior.R")
# Load data -------------------------------------------------------------------
eeco_mod <- readRDS("./results/eeco/eeco_climate_model_output.rds")
ecol_data <- readRDS("./data/processed/distribution_matrix_for_EECO_model.rds")
geochem_data <- readRDS("./data/processed/observation_matrix_for_EECO_model.rds")


# Data preparation ------------------------------------------------------------
eeco_mod <- combine_posterior(mod = eeco_mod, burnin = 100000)

# add proxy type to ecol_data
ecol_data$proxy <- rep("coral reef", 11)
ecol_data$proxy[1:2] <- "Avicennia"
ecol_data$proxy[3:7] <- "A. and Rhizophoraceae"

# Calculate eeco temperatures
lat <- seq(from = 0, to = 90, by = 0.1)
eeco_temp <- temp_from_gradient(lat = lat, model_out = eeco_mod)

