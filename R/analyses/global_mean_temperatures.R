# global mean temperatures for the modern and the Eocene

source("R/functions/model_components/gradient.R")
mode_all <- readRDS("results/eeco/eeco_climate_model_output_params.rds")
modm_all <- readRDS("results/modern/modern_climate_model_output.rds")

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
#
# save results
saveRDS(motemp_global,"results/modern/global_mean_modern.rds")
saveRDS(eotemp_global,"results/eeco/global_mean_eeco.rds")
