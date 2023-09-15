#### Run the climate model for the Early Eocene Climate Optimum
####
# source options for analysis
source("./R/options.R")
#
### Read data
## Hollis data
dat <- readRDS("data/processed/Hollis_processed_EECO_2022_07_19.rds")
#
# prepare for use in the model
obsmat <- data.frame(sample = (paste(abs(dat$p_lat),dat$longitude, dat$proxy)),
                     p_lat = abs(dat$p_lat), temperature = dat$temperature,
                     sd = dat$temperature_sd,
                     proxy = dat$proxy)
# order observation matrix for easily keeping tracks of mean estimates
obsmat <- obsmat[with(obsmat, order(p_lat, sample)),]


### Prepare for model run
# set up clusters for parallel computing
library(doParallel)
nClust <- n_chains # from the options script
cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)

### Run model
# source script
mode <- foreach(pc = 1:n_chains) %dopar% { # run 1 chain per cluster
  # call model functions
  source("R/functions/model_components/climate_model_EECO.R")
  # set random seed
  set.seed(pc)
  hierarchical_model(n_iter = n_iter, n_thin = n_thin,
                     obsmat = obsmat, distrmat = NULL, 
                     prior_input = priors,adapt_sd = adapt_sd)
}

### stop cluster
parallel::stopCluster(cl)

### save output
saveRDS(mode, "results/eeco/eeco_climate_model_output_just_proxy.rds")

# combine chains and discard burn-in
source("R/functions/model_processing/combine_posterior.R")
mode_all <- combine_posterior(mode,100000)

# save combined chains
saveRDS(mode_all, "results/eeco/eeco_climate_model_output_just_proxy_combined.rds")

### To check results use the "processing/assess_EECO-climate_model_output_script.R
