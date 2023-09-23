# Supplementary Figure S2 - gradient on Northern and Southern hemisphere

# load figures for plotting, distributions ect.
source("./R/functions/auxiliary_functions.R")
source("./R/functions/model_components/gradient.R")
source("./R/functions/model_processing/temp_from_gradient.R")

# read model output
mode <- readRDS("./results/eeco/eeco_climate_model_output.rds")
mode_all <- readRDS("./results/eeco/eeco_climate_model_output_combined.rds")

# read the input data
obsmat <- readRDS("./data/processed/observation_matrix_for_EECO_model.rds")
distrmat <- readRDS("./data/processed/distribution_matrix_for_EECO_model.rds")

# calculate location means of observation input data for easier plotting
ylist <- lapply(unique(obsmat$sample), function(x) obsmat$temperature[which(obsmat$sample == x)])
yobs_mean= c(sapply(ylist,mean))
n_obs <- length(yobs_mean)
# add distribution means
mean_t <- c(yobs_mean,distrmat$mu)
lats <- mode[[1]]$lat # latitudes of locations is saved in model output


# read model output
mode_north <- readRDS("./results/SM/eeco_climate_model_northern_hemisphere_output.rds")
mode_all_north <- readRDS("./results/SM/eeco_climate_model_northern_hemisphere_output_combined.rds")
mode_south <- readRDS("./results/SM/eeco_climate_model_southern_hemisphere_output.rds")
mode_all_south <- readRDS("./results/SM/eeco_climate_model_southern_hemisphere_output_combined.rds")

dat <- readRDS("data/processed/Hollis_processed_EECO_2022_07_19.rds")
#
obsmat_north <- data.frame(sample = (paste((dat$p_lat),dat$longitude, dat$proxy)),
                     p_lat = (dat$p_lat), temperature = dat$temperature,
                     sd = dat$temperature_sd,
                     proxy = dat$proxy)
# order observation matrix for easily keeping tracks of mean estimates
obsmat_north <- obsmat_north[which(obsmat$p_lat>=0),]
obsmat_north <- obsmat_north[with(obsmat_north, order(p_lat, sample)),]

obsmat_south<- data.frame(sample = (paste((dat$p_lat),dat$longitude, dat$proxy)),
                           p_lat = (dat$p_lat), temperature = dat$temperature,
                           sd = dat$temperature_sd,
                           proxy = dat$proxy)
# order observation matrix for easily keeping tracks of mean estimates
obsmat_south <- obsmat_south[which(obsmat_south$p_lat<=0),]
obsmat_south$p_lat <- abs(obsmat_south$p_lat)
obsmat_south <- obsmat_south[with(obsmat_south, order(p_lat, sample)),]

# calculate location means of observation input data for easier plotting
ylist_north <- lapply(unique(obsmat_north$sample), function(x) obsmat_north$temperature[which(obsmat_north$sample == x)])
yobs_mean_north= c(sapply(ylist_north,mean))
n_obs_north <- length(yobs_mean_north)
# add distribution means
mean_t_north <- c(yobs_mean_north,distrmat$mu)
lats_north <- mode_north[[1]]$lat # latitudes of locations is saved in model output


# calculate location means of observation input data for easier plotting
ylist_south <- lapply(unique(obsmat_south$sample), function(x) obsmat_south$temperature[which(obsmat_south$sample == x)])
yobs_mean_south= c(sapply(ylist_south,mean))
n_obs_south <- length(yobs_mean_south)
# add distribution means - none
mean_t_south <- c(yobs_mean_south)
lats_south <- mode_south[[1]]$lat # latitudes of locations is saved in model output

lat <- seq(0,90,1)

all_temp <- temp_from_gradient(lat = lat, model_out = mode_all)
north_temp <- temp_from_gradient(lat = lat, model_out = mode_all_north)
south_temp <- temp_from_gradient(lat = lat, model_out = mode_all_south)




### Modern
source("./R/functions/model_processing/temp_from_gradient.R")
source("./R/functions/model_processing/combine_posterior.R")
# Load data --------------------------------------------------------------
temp <- rast("./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")
modern_sample <- readRDS("./results/modern/modern_climate_model_output.rds")
# Data preparation -------------------------------------------------------
# Reduce resolution of raster
r <- rast(res = 1)
temp <- resample(x = temp, y = r)
# Extract temperature data
temp <- terra::as.data.frame(x = temp, xy = TRUE)
colnames(temp) <- c("lng", "lat", "SST")
# Remove data from continents
temp <- na.omit(temp)

# north south seperate
temp_ns <- temp
# Add 1 deg latitudinal bins
temp_ns$bin <- ceiling(temp_ns$lat) - 0.5
# Calculate median per lat bin
grad_ns <- tapply(temp_ns$SST, temp_ns$bin, median)
grad_ns <- data.frame(lat = as.numeric(names(grad_ns)), SST = as.vector(grad_ns))


# Make lats absolute values
temp$lat <- abs(temp$lat)
# Add 1 deg latitudinal bins
temp$bin <- ceiling(temp$lat) - 0.5
# Calculate median per lat bin
grad <- tapply(temp$SST, temp$bin, median)
grad <- data.frame(lat = as.numeric(names(grad)), SST = as.vector(grad))


# Plot north south
png("figures/SM/FigS3_north_south_eeco.png", width = 5.5, height = 3.5, units = "in", res = 400)

par(las = 1, mar = c(4.2,4.2,.5,.5), mgp = c(2.5,0.8,0))
plot(0,0,type = "n", xlim = c(-94,94), ylim = c(-3,40), yaxs = "i", xaxs = "i",
     xlab = expression("latitude ("*degree*")"), ylab = expression("sea surface temperature ("*degree*"C)"), xaxt = "n")

error_polygon(lat, all_temp$l_ci_95, all_temp$u_ci_95, col = rgb(0,0,0,0.2))
error_polygon(-lat, all_temp$l_ci_95, all_temp$u_ci_95, col = rgb(0,0,0,0.2))

error_polygon(lat, north_temp$l_ci_95, north_temp$u_ci_95, col = rgb(0,0,0.8,0.2))

error_polygon(-lat, south_temp$l_ci_95, south_temp$u_ci_95, col = rgb(0.8,0,0,0.2))


points(lat, all_temp$median, type = "l", lwd = 2, col = "black")
points(-lat, all_temp$median, type = "l", lwd = 2, col = "black")

points(lat, north_temp$median, type = "l", lwd = 2, col =  rgb(0,0,0.8))
points(-lat, south_temp$median, type = "l", lwd = 2, col =  rgb(0.8,0,0))

points(grad$lat, grad$SST, lwd = 2, type = "l", lty = 2)
points(-grad$lat, grad$SST, lwd = 2, type = "l", lty = 2)

points(grad_ns$lat[which(grad_ns$lat < 0)], grad_ns$SST[which(grad_ns$lat < 0)], lwd = 2, type = "l", lty = 2,
       col =  rgb(0.8,0,0))
points(grad_ns$lat[which(grad_ns$lat > 0)], grad_ns$SST[which(grad_ns$lat > 0)], lwd = 2, type = "l", lty = 2,
       col =  rgb(0,0,0.8))

# add median of the poster estimates of location means
plot_posterior(mode_north[[1]], col_obs = rgb(0,0,1), col_dist = rgb(0,.75,.72)) # select the first run, doesn't yet work for combined chains
plot_posterior(mode_south[[1]], col_obs = rgb(1,0,0), col_dist = rgb(1,0,0), negative = TRUE) # select the first run, doesn't yet work for combined chains

axis(1,seq(-90,90,30))

text(-87.5,38.25,"Southern Hemisphere", col = rgb(.8,0,0), adj = 0)
text(87.5,38.25,"Northern Hemisphere", col = rgb(0,0,0.8), adj = 1)

legend("bottom", legend = c("EECO", "modern"), lwd = 2, lty = c(1,2), bty = "n")

dev.off()

# 
# ## add the empirical locality means 
# # colour by type 
# cols = rep(rgb(.85,0,0,0.5),length(lats)) 
# cols[(n_obs+1):length(lats)] <- rgb(0,0,.85,0.5) # setb distribution colour to blue
# # plot empirical means
# points(lats,mean_t,pch = 1, col = cols, bg = NA, lwd = 2)
# 
# legend("topright",c("posterior mean (proxy obs.)", "posterior mean (distribution)",
#                     "empirical mean (proxy obs.)", "assigned distribution mean"),
#        pch = c(17,17,1,1), col = c(rgb(0.8,0.5,0,0.75),rgb(0,0.5,.8,0.75),
#                                    rgb(0.85,0,0,0.75),rgb(0,0,.85,0.75)))

