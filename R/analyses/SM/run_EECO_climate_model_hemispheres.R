
###
### seperate gradients- northern hemisphere
###

### Read data
## Hollis data
dat <- readRDS("data/processed/Hollis_processed_EECO_2022_07_19.rds")
#
# prepare for use in the model
obsmat <- data.frame(sample = (paste((dat$p_lat),dat$longitude, dat$proxy)),
                     p_lat = (dat$p_lat), temperature = dat$temperature,
                     sd = dat$temperature_sd,
                     proxy = dat$proxy)
# order observation matrix for easily keeping tracks of mean estimates
obsmat <- obsmat[with(obsmat, order(p_lat, sample)),]
obsmat <- obsmat[which(obsmat$p_lat>=0),]

### Mangrove and Coral data
# read and assign proxy type
bioprox <- readRDS("data/processed/bio_proxies_2022_08_08.RDS")
bioprox$proxy <- rep("Reef",nrow(bioprox))
bioprox$proxy[which(bioprox$taxa=="Avicennia")] <- "Avicennia"
bioprox$proxy[which(bioprox$type=="mangrove" & bioprox$taxa!="Avicennia")] <- 
  "Avicennia-Rhizophoraceae"

# create distribution matrix for use with model (distributions defined in options)
proxy_index <- sapply(bioprox$proxy, function(f) which(proxy_distributions$name==f))

distrmat = data.frame(p_lat = (bioprox$p_lat), 
                      mu = proxy_distributions$mean[proxy_index],
                      scale = proxy_distributions$sd[proxy_index],
                      shape = proxy_distributions$shape[proxy_index],
                      distribution = proxy_distributions$distribution[proxy_index])
distrmat <- distrmat[which(distrmat$p_lat>=0),]

### Prepare for model run
# set up clusters for parallel computing
library(doParallel)
nClust <- 4
cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)

nChains = nClust # run 1 chain per cluster
### Run model
# source script
mode_north <- foreach(pc = 1:nChains) %dopar% {
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
saveRDS(mode_north, "results/SM/eeco_climate_model_northern_hemisphere_output.rds")

# combine chains and discard burn-in
mode_all_north <- combine_posterior(mode_north,100000)
# save combined chains
saveRDS(mode_all_north, "results/SM/eeco_climate_model_northern_hemisphere_output_combined.rds")


###
### seperate gradients- southern hemisphere
###

### Read data
## Hollis data
dat <- readRDS("data/processed/Hollis_processed_EECO_2022_07_19.rds")
#
# prepare for use in the model
obsmat <- data.frame(sample = (paste((dat$p_lat),dat$longitude, dat$proxy)),
                     p_lat = (dat$p_lat), temperature = dat$temperature,
                     sd = dat$temperature_sd,
                     proxy = dat$proxy)
obsmat <- obsmat[which(obsmat$p_lat<=0),]
obsmat$p_lat <- abs(obsmat$p_lat)
# order observation matrix for easily keeping tracks of mean estimates
obsmat <- obsmat[with(obsmat, order(p_lat, sample)),]

### Mangrove and Coral data
# read and assign proxy type
bioprox <- readRDS("data/processed/bio_proxies_2022_08_08.RDS")
bioprox$proxy <- rep("Reef",nrow(bioprox))
bioprox$proxy[which(bioprox$taxa=="Avicennia")] <- "Avicennia"
bioprox$proxy[which(bioprox$type=="mangrove" & bioprox$taxa!="Avicennia")] <- 
  "Avicennia-Rhizophoraceae"

# create distribution matrix for use with model (distributions defined in options)
proxy_index <- sapply(bioprox$proxy, function(f) which(proxy_distributions$name==f))

distrmat = data.frame(p_lat = (bioprox$p_lat), 
                      mu = proxy_distributions$mean[proxy_index],
                      scale = proxy_distributions$sd[proxy_index],
                      shape = proxy_distributions$shape[proxy_index],
                      distribution = proxy_distributions$distribution[proxy_index])
distrmat <- distrmat[which(distrmat$p_lat<=0),]
# no entries

### Prepare for model run
# set up clusters for parallel computing
library(doParallel)
nClust <- 4
cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)

nChains = nClust # run 1 chain per cluster
### Run model
# source script
mode_south <- foreach(pc = 1:nChains) %dopar% {
  # call model functions
  source("R/functions/model_components/climate_model_EECO.R")
  # set random seed
  set.seed(pc)
  hierarchical_model(n_iter = 600000, n_thin = 100,
                     obsmat = obsmat, distrmat = NULL, 
                     prior_input = priors,adapt_sd = 5000)
}

### stop cluster
parallel::stopCluster(cl)

### save output
saveRDS(mode_south, "results/SM/eeco_climate_model_southern_hemisphere_output.rds")

# combine chains and discard burn-in
mode_all_south<- combine_posterior(mode_south,100000)
# save combined chains
saveRDS(mode_all_south, "results/SM/eeco_climate_model_southern_hemisphere_output_combined.rds")
