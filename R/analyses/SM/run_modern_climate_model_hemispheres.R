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

sst_north <- subset(sst,y > 0)
sst_south <- subset(sst,y < 0)


# set up clusters for parallel computing
library(doParallel)
nClust <- 4
cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)

nChains = nClust
### Run model
# source script
mod_north <- foreach(pc = 1:nChains) %dopar% {
  # call model functions
  source("R/functions/model_components/climate_model_modern.R")
  # set random seed
  set.seed(pc)
  run_MCMC_simple(n_iter = 60000,
                  n_thin = 10,
                  x = sst_north$y,
                  y = sst_north$Present.Surface.Temperature.Mean,
                  prior_input = priors,
                  adapt_sd = 5000)
}
                
### stop cluster
parallel::stopCluster(cl)

cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)

nChains = nClust
### Run model
# source script
mod_south <- foreach(pc = 1:nChains) %dopar% {
  # call model functions
  source("R/functions/model_components/climate_model_modern.R")
  # set random seed
  set.seed(pc)
  run_MCMC_simple(n_iter = 60000,
                  n_thin = 10,
                  x = abs(sst_south$y),
                  y = sst_south$Present.Surface.Temperature.Mean,
                  prior_input = priors,
                  adapt_sd = 5000)
}

### stop cluster
parallel::stopCluster(cl)



# discard 10k iterations as burnin each, and combine the results of the 4 chains
source("R/functions/model_processing/combine_posterior.R")
mod_all_north <- combine_posterior(mod_north,burnin = 1000)
mod_all_south <- combine_posterior(mod_south,burnin = 1000)

source("R/functions/model_components/gradient.R")
source("R/functions/auxiliary_functions.R")
plot_gradient(mod_all_north, ylim = c(-4,30))
plot_gradient(mod_all_south, line_col = "red", confint_col = rgb(1,0,0,0.2), add = T)


# save output
saveRDS(mod_all, "results/modern/modern_climate_model_output.rds")
