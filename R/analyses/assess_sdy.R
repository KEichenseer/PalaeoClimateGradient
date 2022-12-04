# compare gradient sdy
source("./R/functions/model_processing/combine_posterior.R")
modern_gradient_eeco_lat <- readRDS("./results/modern/modern_sample_eocene_p_lat_gradient.RDS")
modern_gradient_eeco_lat <- combine_posterior(mod = modern_gradient_eeco_lat, burnin = 5000)
eeco_gradient <- readRDS("./results/eeco/eeco_climate_model_output_combined.rds")

modern_eeco_sample_sdy <- as.list(format(round(quantile(eocene_sample$sdy,probs = c(.025,.5,.975)),1),nsmall = 1))

eeco_sdy <- as.list(format(round(quantile(mode_all$sdy,probs = c(.025,.5,.975)),1),nsmall = 1))

saveRDS(modern_eeco_sample_sdy, "results/modern/modern_eeco_sample_sdy.rds")
saveRDS(eeco_sdy, "results/eeco/eeco_sdy.rds")

# just proxy sdy
eeco_prox <- readRDS("results/eeco/eeco_climate_model_output_just_proxy.rds")
eeco_prox <- combine_posterior(eeco_prox, burnin = 100000)
eeco_prox_sdy <- as.list(format(round(quantile(eeco_prox$sdy,probs = c(.025,.5,.975)),1),nsmall = 1))
