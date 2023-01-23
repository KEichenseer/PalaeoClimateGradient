### calculating results
source("./R/functions/model_processing/temp_from_gradient.R")
source("./R/functions/model_components/gradient.R")

# modern
modm <- readRDS("./results/modern/modern_climate_model_output.rds")

# eocene
mode <- readRDS("./results/eeco/eeco_climate_model_output_combined.rds")

# modern with eocene sampling
modem_list <- readRDS("./results/modern/modern_sample_eocene_p_lat_gradient.rds")
modem <- combine_posterior(modem_list,5000)
  
# calculate temperatures at lat 0 and lat 90 for every iteration
modern_t0 <- gradient(0,modm[,1:4])
modern_eocene_t0 <- gradient(0,modem[,1:4])
eocene_t0 <- gradient(0,mode[,1:4])

modern_t90 <- gradient(90,modm[,1:4])
modern_eocene_t90 <- gradient(90,modm[,1:4])
eocene_t90 <- gradient(90,mode[,1:4])

n <- length(modern_t0)

# randomly recombine iterations of modern and eocene gradient and take difference
set.seed(1)
tdiff_lat0 <- as.list(format(round(quantile(eocene_t0[sample(1:n,n)] - modern_t0[sample(1:n,n)], probs = c(0.025,0.5,0.975)),1),nsmall = 1))
tdiff_lat90 <- as.list(format(round(quantile(eocene_t90[sample(1:n,n)] - modern_t90[sample(1:n,n)], probs = c(0.025,0.5,0.975)),1),nsmall = 1))

# calculate overall gradient for modern and eocene (polar - equatorial difference)
eocene_gradient <- as.list(format(round(quantile(eocene_t0 - eocene_t90, probs = c(0.025,0.5,0.975)),1),nsmall = 1))
modern_gradient <- as.list(format(round(quantile(modern_t0 - modern_t90, probs = c(0.025,0.5,0.975)),1),nsmall = 1))
modern_eocene_gradient <- as.list(format(round(quantile(modern_eocene_t0 - modern_eocene_t90, probs = c(0.025,0.5,0.975)),1),nsmall = 1))

# eocene sample deviation from modern gradient
modern_eocene_gradient_diff <- as.list(format(round(quantile((modern_eocene_t0 - modern_eocene_t90)-(modern_t0 - modern_t90), probs = c(0.025,0.5,0.975)),1),nsmall = 1))

# save results
saveRDS(tdiff_lat0, "./results/eeco/eocene_modern_difference_equator.rds")
saveRDS(tdiff_lat90, "./results/eeco/eocene_modern_difference_poles.rds")
saveRDS(eocene_gradient, "./results/eeco/eocene_overall_gradient.rds")
saveRDS(modern_gradient, "./results/modern/modern_overall_gradient.rds")
saveRDS(modern_eocene_gradient, "./results/modern/modern_eocene_sampling_overall_gradient.rds")
saveRDS(modern_eocene_gradient_diff, "./results/modern/modern_eocene_sampling_overall_gradient_diff.rds")
