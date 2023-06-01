### Run a model on just the priors
library(doParallel)
# source options for analysis
source("./R/options.R")
#
nClust <- 4
cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)
nChains = nClust # run 1 chain per cluster

# source script
modp <- foreach(pc = 1:nChains) %dopar% {
  # call model functions
  source("R/functions/model_components/climate_model_EECO.R")
  # set random seed
  set.seed(pc)
  hierarchical_model(n_iter = 60000, n_thin = 10,
                     obsmat = NULL, distrmat = NULL, 
                     prior_input = priors,adapt_sd = 5000)
}

### stop cluster
parallel::stopCluster(cl)

### save output
saveRDS(modp, "results/eeco/just_prior_model_output.rds")

# combine chains and discard burn-in
modp_all <- combine_posterior(modp,10000)
# save combined chains
saveRDS(modp_all, "results/eeco/just_prior_model_output.rds_combined.rds")
