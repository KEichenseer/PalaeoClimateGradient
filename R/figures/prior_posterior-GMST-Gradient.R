# Supplementary Figure S2 - gradient on Northern and Southern hemisphere

# load figures for plotting, distributions ect.
source("./R/functions/auxiliary_functions.R")
source("./R/functions/model_components/gradient.R")
# read model output
# Eocene
mode_all <- readRDS("./results/eeco/eeco_climate_model_output_combined.rds")
# prior
modp_all <- readRDS("./results/eeco/just_prior_model_output.rds_combined.rds")
# modern with eocene sampling
modem_list <- readRDS("./results/modern/modern_sample_eocene_p_lat_gradient.rds")
modem_all <- combine_posterior(modem_list,5000)



# calculate temperatures at lat 0 and lat 90 for every iteration
eo_t0 <- gradient(0,mode_all[,1:4])
prior_t0 <- gradient(0,modp_all[,1:4])
moeo_t0 <- gradient(0,modem_all[,1:4])

eo_t90 <- gradient(90,mode_all[,1:4])
eo_t67 <- gradient(66.6,mode_all[,1:4])

prior_t90 <- gradient(90,modp_all[,1:4])
moeo_t90 <- gradient(90,modem_all[,1:4])

# gradient
eo_grad <- eo_t0 - eo_t90
eo_grad67 <- eo_t0 - eo_t67

prior_grad <- prior_t0 - prior_t90
moeo_grad <- moeo_t0 - moeo_t90

# densities
d_eo <- density(eo_grad, n=length(eo_grad))
d_eo67 <- density(eo_grad67, n=length(eo_grad67))

d_prior <- density(prior_grad, n=length(prior_grad))
d_moeo <- density(moeo_grad, n=length(moeo_grad))

plot(0,0,type = "n", xlim = c(0,50), ylim = c(0,0.26), xaxs = "i", yaxs = "i", ylab = "posterior density")
error_polygon(d_prior$x,d_prior$y,rep(0,length(d_prior$y)))
points(d_prior$x, d_prior$y, type = "l", col = rgb(0,0,0,0.4), lwd = 2)
error_polygon(d_eo67$x,d_eo67$y,rep(0,length(d_eo67$y)), col = rgb(.8,.6,0,0.2))
points(d_eo67, type = "l", col = rgb(.8,.6,0,0.4), lwd = 2,lty=2)
error_polygon(d_eo$x,d_eo$y,rep(0,length(d_eo$y)), col = rgb(1,0,0,0.2))
points(d_eo, type = "l", col = rgb(1,0,0,0.4), lwd = 2)
error_polygon(d_moeo$x,d_moeo$y,rep(0,length(d_moeo$y)), col = rgb(0,0.5,1,0.2))
points(d_moeo, type = "l", col = rgb(0,0.5,1,0.4), lwd = 2)


### Global mean sea surface temperatures

# latitudinal weights in 1 deg lat bands
alpha1 <- seq(1,90,1)
alpha2 <- seq(0,89,1)
latweight <- sin(pi*alpha1/180) - sin(pi*alpha2/180) # sums to 1 - careful when using lat. subsets


# eocene temperature gradient
eotemp <- gradient(seq(0.5,89.5,1), mode_all[,1:4],0)
# prior temperature gradient
priortemp <- gradient(seq(0.5,89.5,1), modp_all[,1:4],0)
# modern eocene sample temperature gradient
motemp <- gradient(seq(0.5,89.5,1), modem_all[,1:4],0)

# global T
eotemp_global_all <- apply(eotemp,1, function(x) (sum(x*latweight)), simplify = T)
priortemp_global_all <- apply(priortemp,1, function(x) (sum(x*latweight)), simplify = T)
motemp_global_all <- apply(motemp,1, function(x) (sum(x*latweight)), simplify = T)

# densities
d_eo <- density(eotemp_global_all, n=length(eotemp_global_all))
d_prior <- density(priortemp_global_all, n=length(priortemp_global_all))
d_moeo <- density(motemp_global_all, n=length(motemp_global_all))

plot(0,0,type = "n", xlim = c(10,36), ylim = c(0,0.81), xaxs = "i", yaxs = "i", ylab = "posterior density")
error_polygon(d_prior$x,d_prior$y,rep(0,length(d_prior$y)))
points(d_prior, type = "l", col = rgb(0,0,0,0.4), lwd = 2)
error_polygon(d_eo$x,d_eo$y,rep(0,length(d_eo$y)), col = rgb(1,0,0,0.2))
points(d_eo, type = "l", col = rgb(1,0,0,0.4), lwd = 2)
error_polygon(d_moeo$x,d_moeo$y,rep(0,length(d_moeo$y)), col = rgb(0,0.5,1,0.2))
points(d_moeo, type = "l", col = rgb(0,0.5,1,0.4), lwd = 2)


