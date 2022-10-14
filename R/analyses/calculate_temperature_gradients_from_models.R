### calculating results
source("./R/functions/model_processing/temp_from_gradient.R")

# modern
modm <- readRDS("./results/modern/modern_climate_model_output.rds")

# eocene
mode <- readRDS("./results/eeco/eeco_climate_model_output_combined.rds")
# combine chains, discard burnin
mode <- combine_posterior(mode,100000)

# calculate temperatures at lat 0 and lat 90 for every iteration
modern_t0 <- gradient(0,modm[,1:4])
eocene_t0 <- gradient(0,mode[,1:4])

modern_t90 <- gradient(90,modm[,1:4])
eocene_t90 <- gradient(90,mode[,1:4])

n <- length(modern_t0)

# randomly recombine iterations of modern and eocene gradient and take difference
set.seed(1)
quantile(eocene_t0[sample(1:n,n)] - modern_t0[sample(1:n,n)], probs = c(0.025,0.5,0.975))
quantile(eocene_t90[sample(1:n,n)] - modern_t90[sample(1:n,n)], probs = c(0.025,0.5,0.975))

# calculate overall gradient for modern and eocene (polar - equatorial difference)
quantile(eocene_t0 - eocene_t90, probs = c(0.025,0.5,0.975))
quantile(modern_t0 - modern_t90, probs = c(0.025,0.5,0.975))

# calculate global average temperatures

(2*pi*1^2)*(cos(90-40)-cos(90-0))
