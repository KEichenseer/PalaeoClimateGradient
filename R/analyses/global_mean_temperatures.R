# global mean temperatures for the modern and the Eocene

source("R/functions/model_components/gradient.R")
source("R/functions/auxiliary_functions.R")

mode_all <- readRDS("results/eeco/eeco_climate_model_output_combined.rds")
modm_all <- readRDS("results/modern/modern_climate_model_output.rds")
# modern with eocene sampling
modem_list <- readRDS("./results/modern/modern_sample_eocene_p_lat_gradient.rds")
modem_all <- combine_posterior(modem_list,5000)
# just proxy
eeco_prox <- readRDS("results/eeco/eeco_climate_model_output_just_proxy.rds")
eeco_prox <- combine_posterior(eeco_prox, burnin = 100000)

# eeco northern and southern
eeco_north <- readRDS("results/SM/eeco_climate_model_northern_hemisphere_output_combined.rds")
eeco_south <- readRDS("results/SM/eeco_climate_model_southern_hemisphere_output_combined.rds")

# latitudinal weights in 1 deg lat bands
alpha1 <- seq(1,90,1)
alpha2 <- seq(0,89,1)
latweight <- sin(pi*alpha1/180) - sin(pi*alpha2/180) # sums to 1 - careful when using lat. subsets
# eocene temperature gradient
eotemp <- gradient(seq(0.5,89.5,1), mode_all[,1:4],0)
# eocene temperature gradient quantiles
eotempq <- apply(eotemp,2,function(x) quantile(x, probs = c(0.025,0.5,0.975)))
# eocene global average temperature with 95% CI
eotemp_global <- as.list(format(round(quantile(apply(eotemp,1, function(x) sum(x*latweight), simplify = TRUE),probs = c(0.025,0.5,0.975)),1),nsmall = 1))
# modern temperature gradient
motemp <- gradient(seq(0.5,89.5,1), modm_all[,1:4],0)
# modern temperature gradient quantiles
motempq <- apply(motemp,2,function(x) quantile(x, probs = c(0.025,0.5,0.975)))
# modern global average temperature with 95% CI
motemp_global <- as.list(format(round(quantile(apply(motemp,1, function(x) sum(x*latweight), simplify = TRUE),probs = c(0.025,0.5,0.975)),1),nsmall = 1))
# modern eocene sampled temperature gradient
moeotemp <- gradient(seq(0.5,89.5,1), modem_all[,1:4],0)
# modern eocene sampled temperature gradient quantiles
moeotempq <- apply(moeotemp,2,function(x) quantile(x, probs = c(0.025,0.5,0.975)))
# modern eocene sampled global average temperature with 95% CI
moeotemp_global <- as.list(format(round(quantile(apply(moeotemp,1, function(x) sum(x*latweight), simplify = TRUE),probs = c(0.025,0.5,0.975)),1),nsmall = 1))


# eocene proxy only temperature gradient
eotempprox <- gradient(seq(0.5,89.5,1), eeco_prox[,1:4],0)
# eocene proxy only global average temperature with 95% CI
eotempprox_global <- as.list(format(round(quantile(apply(eotempprox,1, function(x) sum(x*latweight), simplify = TRUE),probs = c(0.025,0.5,0.975)),1),nsmall = 1))



# eocene northern temperature gradient
eeconorthtemp <- gradient(seq(0.5,89.5,1), eeco_north[,1:4],0)
# eocene northern temperature gradient quantiles
eeconorthtempq <- apply(eeconorthtemp,2,function(x) quantile(x, probs = c(0.025,0.5,0.975)))
# eocene northern global average temperature with 95% CI
eeconorthtemp_global <- as.list(format(round(quantile(apply(eeconorthtemp,1, function(x) sum(x*latweight), simplify = TRUE),probs = c(0.025,0.5,0.975)),1),nsmall = 1))
#
# eocene southern temperature gradient
eecosouthtemp <- gradient(seq(0.5,89.5,1), eeco_south[,1:4],0)
# eocene southern temperature gradient quantiles
eecosouthtempq <- apply(eecosouthtemp,2,function(x) quantile(x, probs = c(0.025,0.5,0.975)))
# eocene southern global average temperature with 95% CI
eecosouthtemp_global <- as.list(format(round(quantile(apply(eecosouthtemp,1, function(x) sum(x*latweight), simplify = TRUE),probs = c(0.025,0.5,0.975)),1),nsmall = 1))
#

moeotemp_global_all <- apply(moeotemp,1, function(x) (sum(x*latweight)), simplify = T)
motemp_global_all <- apply(motemp,1, function(x) (sum(x*latweight)), simplify = T)
eotemp_global_all <- apply(eotemp,1, function(x) (sum(x*latweight)), simplify = T)
eotempprox_global_all <- apply(eotempprox,1, function(x) (sum(x*latweight)), simplify = T)
eeconorthtemp_global_all <- apply(eeconorthtemp,1, function(x) (sum(x*latweight)), simplify = T)
eecosouthtemp_global_all <- apply(eecosouthtemp,1, function(x) (sum(x*latweight)), simplify = T)

# eeco south north difference 
eeco_south_north_diff_global <- as.list(format(round(quantile(eecosouthtemp_global_all-eeconorthtemp_global_all, probs = c(0.025,0.5,0.975)),1),nsmall = 1))

# modern eocene sample deviation of global average
modern_eocene_gradient_deviation <- as.list(format(round(quantile(moeotemp_global_all-motemp_global_all, probs = c(0.025,0.5,0.975)),1),nsmall = 1))

# eocene proxy deviation of all eocene average
eocene_eoceneproxy_gradient_deviation <- as.list(format(round(quantile(eotempprox_global_all-eotemp_global_all, probs = c(0.025,0.5,0.975)),1),nsmall = 1))

# save results
saveRDS(motemp_global,"results/modern/global_mean_modern.rds")
saveRDS(eotemp_global,"results/eeco/global_mean_eeco.rds")
saveRDS(moeotemp_global,"results/eeco/global_mean_modern_eocene_sampling_eeco.rds")
saveRDS(eotempprox_global,"results/eeco/global_mean_eeco_proxy_only.rds")

saveRDS(modern_eocene_gradient_deviation,"results/eeco/global_mean_modern_eocene_sampling_difference.rds")
saveRDS(eocene_eoceneproxy_gradient_deviation,"results/eeco/global_mean_eocene_eoceneproxy_difference.rds")
saveRDS(eeconorthtemp_global,"results/SM/global_mean_northern_eeco.rds")
saveRDS(eecosouthtemp_global,"results/SM/global_mean_southern_eeco.rds")
saveRDS(eeco_south_north_diff_global,"results/SM/global_mean_diff_south_north_eeco.rds")
