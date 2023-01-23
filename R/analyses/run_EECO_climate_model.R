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

saveRDS(obsmat, "./data/processed/observation_matrix_for_EECO_model.rds")

### Mangrove and Coral data
# read and assign proxy type
bioprox <- readRDS("data/processed/bio_proxies_2022_08_08.RDS")
bioprox$proxy <- rep("Reef",nrow(bioprox))
bioprox$proxy[which(bioprox$taxa=="Avicennia")] <- "Avicennia"
bioprox$proxy[which(bioprox$type=="mangrove" & bioprox$taxa!="Avicennia")] <- 
  "Avicennia-Rhizophoraceae"

# create distribution matrix for use with model (distributions defined in options)
proxy_index <- sapply(bioprox$proxy, function(f) which(proxy_distributions$name==f))

distrmat = data.frame(p_lat = abs(bioprox$p_lat), 
                      mu = proxy_distributions$mean[proxy_index],
                      scale = proxy_distributions$sd[proxy_index],
                      shape = proxy_distributions$shape[proxy_index],
                      distribution = proxy_distributions$distribution[proxy_index])
saveRDS(distrmat, "./data/processed/distribution_matrix_for_EECO_model.rds")
### Prepare for model run
# set up clusters for parallel computing
library(doParallel)
nClust <- 4
cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)

nChains = nClust # run 1 chain per cluster
### Run model
# source script
mode <- foreach(pc = 1:nChains) %dopar% {
  # call model functions
  source("R/functions/model_components/climate_model_EECO.R")
  # set random seed
  set.seed(pc)
  hierarchical_model(n_iter = 600000, n_thin = 100,
                  obsmat = obsmat, distrmat = distrmat, 
                  prior_input = priors,adapt_sd = 5000)
}

### stop cluster
parallel::stopCluster(cl)

### save output
saveRDS(mode, "results/eeco/eeco_climate_model_output.rds")

# combine chains and discard burn-in
mode_all <- combine_posterior(mode,100000)
# save combined chains
saveRDS(mode_all, "results/eeco/eeco_climate_model_output_combined.rds")

### To check results use the "assess_EECO-climate_model_output_script.R
