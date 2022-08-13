#### Run a non-hierarchical climate model on the modern data
#
# source options for analysis
source("R/options.R")
#
### load climate data
# Mean annual sea surface temperatures (Bio-Oracle)
set.seed(1)

sst <- raster::raster("./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")
r <- raster::raster(res = 1)
sst <- raster::resample(x = sst, y = r)

sst <- raster::as.data.frame(x = sst, xy = TRUE, centroids = TRUE)
sst <- na.omit(sst)

temp <- sst$Present.Surface.Temperature.Mean

# use absolute latitudes
lat <- abs(sst$y)

# set up clusters for parallel computing
library(doParallel)
nClust <- 4
cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)

nChains = nClust
### Run model
# source script
modm <- foreach(pc = 1:nChains) %dopar% {
  # call model functions
  source("R/functions/model_components/climate_model_modern.R")
  # set random seed
  set.seed(pc)
  run_MCMC_simple(n_iter = 1000,
                  n_thin = 10,
                  x = lat,
                  y = temp,
                  prior_input = priors,
                  adapt_sd = NULL)
}
                
### stop cluster
parallel::stopCluster(cl)


# discard 10k iterations as burnin each, and combine the results of the 4 chains
mod_all <- combine_posterior(modm,burnin = 10000)

# save output
saveRDS(mod_all, "R/subscripts/models/modern_climate_model_output.rds")