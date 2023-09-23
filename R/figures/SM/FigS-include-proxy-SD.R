# Supplementary Figure S2 - gradient with just geochem.
source("./R/options.R")

# load figures for plotting, distributions ect.
source("./R/functions/auxiliary_functions.R")
source("./R/functions/model_components/gradient.R")
source("./R/functions/model_processing/temp_from_gradient.R")

# read model output
mode <- readRDS("./results/eeco/eeco_climate_model_output.rds")
mode_all <- readRDS("./results/eeco/eeco_climate_model_output_combined.rds")

mode_sd <-readRDS("results/SM/eeco_climate_model__with_observation_sd_output.rds")
mode_all_sd <-readRDS("results/SM/eeco_climate_model__with_observation_sd_output_combined.rds")


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


cex1 <- c(1.05,1,1.15,0.92)
alpha2 <- 0.75
pch1 <- c(16,17,18,15)
proxind1 <- sapply(obsmat$proxy, function(x) which(x == c("d18O","MgCa","TEX86","d47")))
alpha3 <- 0.9
col_proxgrad <- rgb(0.9,0.3,0,0.9)
col_proxgrad_shade <- rgb(0.9,0.3,0,0.27)
xlab1 <- expression("absolute latitude ("*degree*")")
ylab1 <- expression("sea surface temperature ("*degree*"C)")
lat_obs <- mode[[1]]$lat
alpha2 <- 0.75
col1 <- c(rgb(.9,.8,.5,alpha2), rgb(.8,.7,.5,alpha2), rgb(.7,.6,.5,alpha2), rgb(.65,.6,.6,alpha2))


lat <- seq(0,90,1)

all_temp <- temp_from_gradient(lat = lat, model_out = mode_all)
all_sd_temp <- temp_from_gradient(lat = lat, model_out = mode_all_sd)



png("figures/SM/FigS_sd_included.png", width = 4.75, height = 3.5, units = "in", res = 400)
par(mar=c(4.2,4.2,0.5,0.5), mgp = c(2.5,0.8,0)
    ,las = 1)
plot(obsmat$p_lat,obsmat$temperature, pch = pch1[proxind1], col = col1[proxind1], xlim = c(-2,92), ylim = c(0,41.5),
     xlab = xlab1, ylab = ylab1, cex = cex1[proxind1], type = "n", xaxs = "i", yaxs = "i")

error_polygon(lat, all_temp$l_ci_95, all_temp$u_ci_95, col =  rgb(0,0,0,0.2))

error_polygon(lat, all_sd_temp$l_ci_95, all_sd_temp$u_ci_95, col = rgb(0,.6,1,0.15))


points(lat, all_temp$median, type = "l", lwd = 2, col = "black")
points(lat, all_sd_temp$median, type = "l", lwd = 2, col = rgb(0,.6,1,0.75))

plot_posterior(mode_sd[[1]], col_obs = rgb(0,0.6,1,1), col_dist = NA, add=T)# line_col = rgb(0,0,1,0.5), confint_col = rgb(0,0,1,0.1))

plot_posterior(mode[[1]], col_obs = rgb(0,0,0,0.5), col_dist = NA, add=T,
               pch = 19)# line_col = rgb(0,0,1,0.5), confint_col = rgb(0,0,1,0.1))

legend("bottomleft", legend = c("original analysis", "with proxy uncertainty"), fill = c(rgb(0,0,0,0.75), rgb(0,0.6,1,1)),
       bty = "n", col = NA)

dev.off()
