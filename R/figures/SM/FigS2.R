# Supplementary Figure S2 - gradient on Northern and Southern hemisphere

# load figures for plotting, distributions ect.
source("./R/functions/auxiliary_functions.R")
source("./R/functions/model_components/gradient.R")
# read model output
mode <- readRDS("./results/eeco/eeco_climate_model_output.rds")
mode_all <- readRDS("./results/eeco/eeco_climate_model_output_combined.rds")
mode_prox <- readRDS("./results/eeco/eeco_climate_model_output_just_proxy.rds")
mode_all_prox <- combine_posterior(mode_prox, burnin = 100000)

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
col_proxgrad <- rgb(0.85,0.6,0,0.9)
col_proxgrad_shade <- rgb(0.85,0.6,0,0.27)
xlab1 <- expression("|latitude| ("*degree*")")
ylab1 <- expression("temperature ("*degree*"C)")
lat_obs <- mode[[1]]$lat
alpha2 <- 0.75
col1 <- c(rgb(.9,.8,.5,alpha2), rgb(.8,.7,.5,alpha2), rgb(.7,.6,.5,alpha2), rgb(.65,.6,.6,alpha2))
eecoprox <- readRDS("results/eeco/eeco_climate_model_output_just_proxy.rds")
eecoproxcomb <- combine_posterior(eecoprox, burnin = 100000)


mu_obs <- apply(eecoprox[[1]]$yestimate[1001:6000,],2,mean)
low_obs <- apply(eecoprox[[1]]$yestimate[1001:6000,],2,function(x) quantile(x,0.025))
upp_obs <- apply(eecoprox[[1]]$yestimate[1001:6000,],2,function(x) quantile(x,0.975))
ind1_mu <- sapply(unique(obsmat$sample), function(x) which(unique(obsmat$proxy[which(obsmat$sample==x)]) == c("d18O","MgCa","TEX86","d47")))




png("figures/SM/FigS2_geochem_eeco.png", width = 4.75, height = 3.5, units = "in", res = 400)
par(mar=c(4.2,4.2,0.5,0.5), mgp = c(2.5,0.8,0)
    ,las = 1)
plot(obsmat$p_lat,obsmat$temperature, pch = pch1[proxind1], col = col1[proxind1], xlim = c(-2,92), ylim = c(14,41.5),
     xlab = xlab1, ylab = ylab1, cex = cex1[proxind1], type = "n", xaxs = "i", yaxs = "i")

plot_gradient(mode_all,add=T, line_col = rgb(0,0.3,0.8,0.5), confint_col = rgb(0,0.3,0.8,0.15))

plot_gradient(eecoproxcomb,add=T, line_col = col_proxgrad, confint_col = col_proxgrad_shade)

points(lat_obs[1:23], mu_obs[1:23], col = col1[ind1_mu], pch = pch1[ind1_mu], lwd = 3, cex = 0.1+cex1[ind1_mu])
sapply(1:23, function(x) points(rep(lat_obs[x],2), c(low_obs[x],upp_obs[x]), type = "l", lwd = 2, col = col1[ind1_mu[x]]))

points(lat_obs[1:23], mu_obs[1:23], col = col1[ind1_mu], pch = pch1[ind1_mu], lwd = 3, cex = 0.1+cex1[ind1_mu])
sapply(1:23, function(x) points(rep(lat_obs[x],2), c(low_obs[x],upp_obs[x]), type = "l", lwd = 2, col = col1[ind1_mu[x]]))

text(15,39.5,"only geochemical proxies", col = col_proxgrad, adj = 0)
text(5,26,"all proxies", col = rgb(0,0.3,0.8,0.5), adj = 0)

dev.off()
