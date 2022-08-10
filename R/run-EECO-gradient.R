
### Hollis data

dat <- readRDS("data/processed/Hollis_processed_EECO_2022_07_19.rds")



# prepare for use in the model
obsmat <- data.frame(sample = (paste(abs(data_sub$p_lat),data_sub$longitude, data_sub$proxy)),
                     latitude = abs(data_sub$p_lat), temperature = data_sub$temperature,
                     sd = data_sub$temperature_sd,
                     proxy = data_sub$proxy)

### Mangrove and Coral data
# read and assign proxy type
bioprox <- readRDS("data/processed/bio_proxies_2022_08_08.RDS")
bioprox$proxy <- rep("Reef",nrow(bioprox))
bioprox$proxy[which(bioprox$taxa=="Avicennia")] <- "Avicennia"
bioprox$proxy[which(bioprox$type=="mangrove" & bioprox$taxa!="Avicennia")] <- "Avicennia-Rhizophoraceae"

# define proxy distributions
proxy_distributions <- data.frame(name = c("Avicennia", "Avicennia-Rhizophoraceae", "Reef"),
               distribution = c("normal", "normal", "normal"),
                mean = c(mean(c(15.6,22.5)), mean(c(20.7,29.5)), 27.6),
                sd = c((22.5-15.6)/4, c(29.5-20.7)/4, (29.5-21)/4),
                shape = rep(NA,3))

# create distribution matrix for use with model
proxy_index <- sapply(bioprox$proxy, function(f) which(proxy_distributions$name==f))

distrmat = data.frame(latitude = abs(bioprox$p_lat), 
                      location = proxy_distributions$mean[proxy_index],
                      scale = proxy_distributions$sd[proxy_index],
                      shape = proxy_distributions$shape[proxy_index],
                      distribution = proxy_distributions$distribution[proxy_index])

### Define the priors
prior_fun <- readRDS("results/modern/prior_from_modern_gradient.RDS")

xval <- list(seq(-5,35,0.1),
             seq(-5,60,0.1),
             seq(0,90,0.1),
             seq(0,0.6,0.01))
source("R/subscripts/AuxiliaryFunctions.R")
plot_prior(prior_fun,xval)
### Prepare for model run
# Source model script
source("R/subscripts/ClimateParallelSD.R")

# set up clusters for parallel computing
library(doParallel)
nClust <- 4
cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)

### Run model
mode <- climate_parallel_sd(nChains = nClust, nIter = 25000, nThin = 5, obsmat = obsmat, distrmat = distrmat, 
                            coeff_inits = NULL, sdy_init = NULL, 
                            yest_inits = NULL, sdyest_inits = NULL, prior_fun = prior_fun,
                            proposal_var_inits = c(3,3,3,0.2), adapt_sd = 2500,
                            adapt_sd_decay = 100, start_adapt = 101, quiet = FALSE) 
# Run model without bio proxies
mode_nobio <- climate_parallel_sd(nChains = nClust, nIter = 25000, nThin = 5, obsmat = obsmat, distrmat = NULL, 
                            coeff_inits = NULL, sdy_init = NULL, 
                            yest_inits = NULL, sdyest_inits = NULL, prior_fun = prior_fun,
                            proposal_var_inits = c(3,3,3,0.2), adapt_sd = 2500,
                            adapt_sd_decay = 100, start_adapt = 101, quiet = FALSE) 

# Run model with no sd on observations
obsmat$sd <- NA
mode_nosd <- climate_parallel_sd(nChains = nClust, nIter = 25000, nThin = 5, obsmat = obsmat, distrmat = distrmat, 
                                  coeff_inits = NULL, sdy_init = NULL, 
                                  yest_inits = NULL, sdyest_inits = NULL, prior_fun = prior_fun,
                                  proposal_var_inits = c(3,3,3,0.2), adapt_sd = 2500,
                                  adapt_sd_decay = 100, start_adapt = 101, quiet = FALSE) 


### stop cluster
parallel::stopCluster(cl)

#######################################
#
# Assess output

plot_chains(mode_nosd)

mode_all <- combine_posterior(mode,4000)
mode_all_nobio <- combine_posterior(mode_nobio,4000)
mode_all_nosd <- combine_posterior(mode_nosd,4000)


mcmcse::multiESS(mode_all[,1:4])
plot_gradient(mode_all,ylim = c(10,37))
plot_gradient(mode_all_nobio,ylim = c(10,37), confint_col = rgb(0,.7,0.5,0.2),line_col = rgb(0,.7,0.5,1), add = T)
plot_gradient(mode_all_nosd,ylim = c(10,37), confint_col = rgb(0.7,0,0.5,0.2),line_col = rgb(0.7,0,0.5,1), add = T)

plot_posterior(mode[[3]])

