# global mean temperatures for the modern and the Eocene

source("R/functions/model_components/gradient.R")
mode_all <- readRDS("results/eeco/eeco_climate_model_output_combined.rds")
modm_all <- readRDS("results/modern/modern_climate_model_output.rds")
# modern with eocene sampling
modem_list <- readRDS("./results/modern/modern_sample_eocene_p_lat_gradient.rds")
modem_all <- combine_posterior(modem_list,5000)

# latitudinal weights in 1 deg lat bands
alpha1 <- seq(1,90,1)
alpha2 <- seq(0,89,1)
latweight <- sin(pi*alpha1/180) - sin(pi*alpha2/180) # sums to 1 - careful when using lat. subsets
# eocene temperature gradient
eotemp <- gradient(seq(0.5,89.5,1), mode_all[,1:4],0)
# eocene temperature gradient quantiles
eotempq <- apply(eotemp,2,function(x) quantile(x, probs = c(0.025,0.5,0.975)))
# eocene global average temperature with 95% CI
eotemp_global <- apply(eotempq,1, function(x) round(sum(x*latweight),1), simplify = FALSE)
# modern temperature gradient
motemp <- gradient(seq(0.5,89.5,1), modm_all[,1:4],0)
# modern temperature gradient quantiles
motempq <- apply(motemp,2,function(x) quantile(x, probs = c(0.025,0.5,0.975)))
# modern global average temperature with 95% CI
motemp_global <- apply(motempq,1, function(x) round(sum(x*latweight),1), simplify = FALSE)
# modern eocene sampled temperature gradient
moeotemp <- gradient(seq(0.5,89.5,1), modem_all[,1:4],0)
# modern eocene sampled temperature gradient quantiles
moeotempq <- apply(moeotemp,2,function(x) quantile(x, probs = c(0.025,0.5,0.975)))
# modern eocene sampled global average temperature with 95% CI
moeotemp_global <- apply(moeotempq,1, function(x) round(sum(x*latweight),1), simplify = FALSE)
#

moeotemp_global_all <- apply(moeotemp,1, function(x) (sum(x*latweight)), simplify = T)
motemp_global_all <- apply(motemp,1, function(x) (sum(x*latweight)), simplify = T)

# modern eocene sample deviation of global average
modern_eocene_gradient_deviation <- as.list(format(round(quantile(moeotemp_global_all-motemp_global_all, probs = c(0.025,0.5,0.975)),1),nsmall = 1))

# save results
saveRDS(motemp_global,"results/modern/global_mean_modern.rds")
saveRDS(eotemp_global,"results/eeco/global_mean_eeco.rds")
saveRDS(moeotemp_global,"results/eeco/global_mean_modern_eocene_sampling_eeco.rds")
saveRDS(modern_eocene_gradient_deviation,"results/eeco/global_mean_modern_eocene_sampling_difference.rds")
