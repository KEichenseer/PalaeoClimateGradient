
### Hollis data

dat <- readRDS("data/processed/Hollis_processed_2022_07_19.rds")

# select data from one stage to test, exclude NA data
data_sub <- subset(dat,EECO == 1 & !(is.na(temperature)) & !(is.na(paleolat_Meredith)) &
                     (is.na(depth_habitat) | depth_habitat %in% c("Mixed-layer",     "Mixed layer")))

# order data
data_sub <- data_sub[with(data_sub, order(abs(paleolat_Meredith), longitude)),]

# prepare for use in the model
obsmat <- data.frame(sample = (paste(abs(data_sub$paleolat_Meredith),data_sub$longitude, data_sub$proxy)),
                     latitude = abs(data_sub$paleolat_Meredith), temperature = data_sub$temperature,
                     sd = data_sub$temperature_sd,
                     proxy = data_sub$proxy)

### Mangrove and Coral data
# for now use this placeholder

bioprox <- readRDS("data/processed/bio_proxies_2022_08_08.RDS")

distrmat = data.frame(latitude = 79.4 , ## palaeorotated from Faddeevsky Island: palaeorotate(data.frame(lat = 75.5, lng = 144, age = 52))
                      location = mean(c(15.6,20.8)),
                      scale = 1.33,
                      shape = NA,
                      distribution = "normal")

### Define the priors
prior_fun <- list(  
  f1 = function(x,log) dsnorm(x,location = -2.66, scale = 20, alpha = 20, log = log), # prior on A (lower asymptote)
  f2 = function(x,lower,log) dtnorm(x, lower, upper = Inf, mean = 30, sd = 10, log = log), # prior on upper asymptote
  f3 = function(x,log) dnorm(x, mean = 45, sd = 12, log = log), # prior on M
  f4 = function(x,log) dlnorm(x, mean = -2.4, sd = 0.6, log = log)# prior on Q
)

xval <- list(seq(-5,35,0.1),
             seq(-5,60,0.1),
             seq(0,90,0.1),
             seq(0,0.5,0.01))
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
mod2 <- climate_parallel_sd(nChains = nClust, nIter = 10000, nThin = 5, obsmat = obsmat, distrmat = distrmat, 
                            coeff_inits = NULL, sdy_init = NULL, 
                            yest_inits = NULL, sdyest_inits = NULL, prior_fun = prior_fun,
                            proposal_var_inits = c(2,2,2,0.2), adapt_sd = FALSE,
                            adapt_sd_decay = 100, start_adapt = 101, quiet = FALSE) 


### stop cluster
parallel::stopCluster(cl)

#######################################
#
# Assess output
source("R/subscripts/AuxiliaryFunctions.R")

plot_chains(modm)
plot_gradient(mod1[[1]])
plot_gradient(mod2[[1]], add = T, line_col = "red",confint_col = rgb(1,0,0,0.2))

