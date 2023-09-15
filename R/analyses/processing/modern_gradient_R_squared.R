# Load libraries ---------------------------------------------------------
library(terra)
library(ggplot2)
library(ggpubr)
source("./R/functions/model_processing/temp_from_gradient.R")
source("./R/functions/model_processing/combine_posterior.R")
# Load data --------------------------------------------------------------
temp <- rast("./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")
modern_sample <- readRDS("./results/modern/modern_climate_model_output.rds")
eocene_sample <- readRDS("./results/modern/modern_sample_eocene_p_lat_gradient.RDS")
eocene_reps <- readRDS("./results/modern/modern_sample_summary_stats.RDS")
# Data preparation -------------------------------------------------------
# Reduce resolution of raster
r <- rast(res = 1)
temp <- resample(x = temp, y = r)
# Extract temperature data
temp <- terra::as.data.frame(x = temp, xy = TRUE)
colnames(temp) <- c("lng", "lat", "SST")
# Remove data from continents
temp <- na.omit(temp)
# Make lats absolute values
temp$lat <- abs(temp$lat)
eocene_reps$p_lat <- abs(eocene_reps$p_lat)
# Add 1 deg latitudinal bins
temp$bin <- ceiling(temp$lat) - 0.5
# Calculate median per lat bin
grad <- tapply(temp$SST, temp$bin, median)
grad <- data.frame(lat = as.numeric(names(grad)), SST = as.vector(grad))
# save for later use
saveRDS(grad,"results/modern/empirical_median_1-deg-lat.rds")
# Calculate modern temperatures from sample
lat <- seq(from = 0, to = 90, by = 0.1)
modern_gradient <- temp_from_gradient(lat = lat, model_out = modern_sample)
# Calculate modern temperatures from palaeodistribution of Eocene data
#eocene_sample <- combine_posterior(mod = eocene_sample, burnin = 5000)
#eocene_sample <- temp_from_gradient(lat = lat, model_out = eocene_sample)
# Plot ------------------------------------------------------------------------

#if(s == 4)     grad <- gradient(lats, four_gradients_model[[g]][[s+1]])

true_grad <- gradient_formulas_static[[g]](lats)

var_fit <- var(modern_gradient$median)
var_res <- var(temp$SST - gradient(temp$lat, apply(modern_sample,2,median)[1:4]))

R_squared_gelman <- var_fit/(var_fit+var_res)

saveRDS(R_squared_gelman, "results/modern/Rsquared_modern.rds")
