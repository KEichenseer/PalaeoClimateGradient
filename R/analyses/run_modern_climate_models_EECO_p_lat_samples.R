### Run 100 analyses with 100 modern samples, based on the Eocene palaeolatitudes
## -- non-hierarchical climate model on the modern data
#
# source options for analysis
source("R/options.R")
#
# read modern temperatures, sampled at the palaeolatitudes of the Eocene record
modt_samples <- readRDS("results/modern/modern_sample.RDS")

# source model script
source("R/functions/model_components/climate_model_modern.R")

# list to store result
library(doParallel)
nClust <- 4
cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)

nChains = 4
### Run model
mods <- foreach(s = 1:nChains) %dopar% {
  # load model functions
  source("R/functions/model_components/climate_model_modern.R")
  # set random seed
  set.seed(s)
  p_lat = abs(modt_samples[[s]][,1])
  temp = modt_samples[[s]][,2]
  
  coeff_inits <- c(rnorm(3,c(5,30,45), c(3,5,7)),exp(rnorm(1,log(0.1),0.6)))
  sdy_init <- exp(rnorm(1,log(2),0.6))
  run_MCMC_simple(n_iter = 25000, n_thin = 10,
                               x = p_lat, y = temp, 
                               coeff_inits = coeff_inits, sdy_init = sdy_init, 
                               prior_input  = priors, adapt_sd = 2500,
                               adapt_sd_decay = 100, start_adapt = 101, quiet = FALSE)
}

### stop cluster
parallel::stopCluster(cl)

### save results
saveRDS(mods,"results/modern/modern_sample_eocene_p_lat_gradient.rds")
