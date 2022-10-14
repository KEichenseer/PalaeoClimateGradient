### calculating results
source("./R/functions/model_processing/temp_from_gradient.R")

# modern
modm <- readRDS("./results/modern/modern_climate_model_output.rds")

# eocene
mode <- readRDS("./results/eeco/eeco_climate_model_output_combined.rds")


# calculate temperatures at lat 0 and lat 90 for every iteration
modern_t0 <- gradient(0,modm[,1:4])
eocene_t0 <- gradient(0,mode[,1:4])

modern_t90 <- gradient(90,modm[,1:4])
eocene_t90 <- gradient(90,mode[,1:4])

n <- length(modern_t0)

# randomly recombine iterations of modern and eocene gradient and take difference
set.seed(1)
tdiff_lat0 <- as.list(quantile(eocene_t0[sample(1:n,n)] - modern_t0[sample(1:n,n)], probs = c(0.025,0.5,0.975)))
tdiff_lat90 <- as.list(quantile(eocene_t90[sample(1:n,n)] - modern_t90[sample(1:n,n)], probs = c(0.025,0.5,0.975)))

# calculate overall gradient for modern and eocene (polar - equatorial difference)
eocene_gradient <- as.list(quantile(eocene_t0 - eocene_t90, probs = c(0.025,0.5,0.975)))
modern_gradient <- as.list(quantile(modern_t0 - modern_t90, probs = c(0.025,0.5,0.975)))

# save results
