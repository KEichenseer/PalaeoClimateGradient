#### Run the climate model for the Early Eocene Climate Optimum
####
# source options for analysis
source("./R/options.R")
#
# change prior on M
priors$f3 <- function(x,log) dunif(x, min = 0, max = 90, log = log) # dnorm(x, mean = 42, sd = 25, log = log) #
# change prior on K
priors$f2 <- function(x,lower,log) dtnorm(x, lower, upper = Inf, mean = 28, sd = 15, log = log)

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

### Prepare for model run
# set up clusters for parallel computing
library(doParallel)
nClust <- n_chains # from the options script
cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)

### Run model
# source script
mode3 <- foreach(pc = 1:n_chains) %dopar% { # run 1 chain per cluster
  # call model functions
  source("R/functions/model_components/climate_model_EECO.R")
  # set random seed
  set.seed(pc)
  hierarchical_model(n_iter = 0.1*n_iter, n_thin = 0.1*n_thin,
                  obsmat = obsmat, distrmat = distrmat, 
                  prior_input = priors,adapt_sd = adapt_sd)
}

### stop cluster
parallel::stopCluster(cl)

### save output
saveRDS(mode, "results/SM/eeco_climate_model_output_different_M_prior.rds")

# combine chains and discard burn-in
source("R/functions/model_processing/combine_posterior.R")
mode_all <- combine_posterior(mode,100000*0.1)
mode_all2 <- combine_posterior(mode2,100000*0.1)
mode_all3 <- combine_posterior(mode3,100000*0.1)

plot_gradient(mode_all, ylim = c(-5,37))
plot_gradient(mode_all3, add = T, line_col = rgb(0,1,1,.7), confint_col = rgb(0,0,1,0.2))
plot_gradient(mode_all2, add = T, line_col = rgb(1,0,0,.7), confint_col = rgb(1,0,0,0.2))

plot_chains(mode)

# save combined chains
saveRDS(mode_all, "results/eeco/eeco_climate_model_output_combined.rds")

### To check results use the "processing/assess_EECO-climate_model_output_script.R


eeco_temp_ori <- temp_from_gradient(lat = 0:90, model_out = mode_all)
eeco_temp_3 <- temp_from_gradient(lat = 0:90, model_out = mode_all3)

points(eeco_temp_ori$lat,eeco_temp_ori$median, col = "black", type = "l", ylim = c(-3,36))
error_polygon(eeco_temp_ori$lat,eeco_temp_ori$l_ci_95, eeco_temp_ori$u_ci_95, col = rgb(0,0,0,0.2))

points(eeco_temp_3$lat,eeco_temp_3$median, col = "blue", type = "l", ylim = c(-3,36))
error_polygon(eeco_temp_3$lat,eeco_temp_3$l_ci_95, eeco_temp_3$u_ci_95, col = rgb(0,1,1,0.2))


legend("bottomleft", legend = c("original", "relaxed priors on M (and K)"), col = c("black", "red"),
       lwd = 1)
