# Supplementary Figure S2 - gradient with just geochem.
source("./R/options.R")

# load figures for plotting, distributions ect.
source("./R/functions/auxiliary_functions.R")
source("./R/functions/model_components/gradient.R")
# read model output
mode <- readRDS("./results/eeco/eeco_climate_model_output.rds")
mode_all <- readRDS("./results/eeco/eeco_climate_model_output_combined.rds")
mode_ecowide <- readRDS("./results/SM/eeco_climate_model_wider_ecol.rds")
mode_all_ecowide <-  readRDS("./results/SM/eeco_climate_model_wider_ecol_combined.rds")

# read the input data
obsmat <- readRDS("./data/processed/observation_matrix_for_EECO_model.rds")
distrmat <- readRDS("./data/processed/distribution_matrix_for_EECO_model.rds")


wider_coral_mean <- mean(c(35.6,16.0))
wider_coral_sd <- round(1/(1.96*2) * (35.6-16.0),2)
proxy_distributions$mean[3] <- wider_coral_mean
proxy_distributions$sd[3] <- wider_coral_sd

wider_avi_rhi_mean <- wider_coral_mean
wider_avi_rhi_sd <- wider_coral_sd
proxy_distributions$mean[2] <- wider_avi_rhi_mean
proxy_distributions$sd[2] <- wider_avi_rhi_sd

distrmat2 = data.frame(p_lat = abs(bioprox$p_lat), 
                      mu = proxy_distributions$mean[proxy_index],
                      scale = proxy_distributions$sd[proxy_index],
                      shape = proxy_distributions$shape[proxy_index],
                      distribution = proxy_distributions$distribution[proxy_index])



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
xlab1 <- expression("|latitude| ("*degree*")")
ylab1 <- expression("temperature ("*degree*"C)")
lat_obs <- mode[[1]]$lat
alpha2 <- 0.75
col1 <- c(rgb(.9,.8,.5,alpha2), rgb(.8,.7,.5,alpha2), rgb(.7,.6,.5,alpha2), rgb(.65,.6,.6,alpha2))



png("figures/SM/FigS_wider-ecol-limits.png", width = 4.75, height = 3.5, units = "in", res = 400)
par(mar=c(4.2,4.2,0.5,0.5), mgp = c(2.5,0.8,0)
    ,las = 1)
plot(obsmat$p_lat,obsmat$temperature, pch = pch1[proxind1], col = col1[proxind1], xlim = c(-2,92), ylim = c(14,41.5),
     xlab = xlab1, ylab = ylab1, cex = cex1[proxind1], type = "n", xaxs = "i", yaxs = "i")

plot_gradient(mode_all,add=T, line_col = rgb(0,0,0,1), confint_col = rgb(0,0,0,0.2))

plot_gradient(mode_all_ecowide, add=T, line_col = col_proxgrad, confint_col = col_proxgrad_shade)



sapply(1:11, function(x) points(rep(distrmat2$p_lat[x],2), distrmat2$mu[x] + distrmat2$scale[x] * c(-2,2), type = "l", lwd = 2,
                                col = rgb(0.9,0.3,0,0.75)))

sapply(1:11, function(x) points(rep(distrmat$p_lat[x],2), distrmat$mu[x] + distrmat$scale[x] * c(-2,2), type = "l", lwd = 2,
                                col = rgb(0,0,0,0.5)))

points(distrmat$p_lat, distrmat$mu, col = rgb(0,0,0,0.75), pch = c(rep(17,2), rep(18,5), rep(19,4)),
       lwd = 1, cex = 1)

points(distrmat2$p_lat, distrmat2$mu, col = rgb(0.9,0.3,0,0.75), pch = c(rep(17,2), rep(18,5), rep(19,4)),
       lwd = 1, cex = 1)


legend("topright", legend = c("original analysis", "wider ecological limits"), fill = c(rgb(0,0,0,0.75), rgb(0.9,0.3,0,0.75)),
       bty = "n")

dev.off()
