### load climate data
# Mean annual sea surface temperatures (Bio-Oracle)

sst <- raster("./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")

# sample points randomly from the raster (within specified latitudinal limits) 
source("R/subscripts/AuxiliaryFunctions.R")
coords <- fastRandomPoints_lat(sst$Present.Surface.Temperature.Mean,2000,-90,90)

# extract temperatures from the sampled points
sstr <- raster::extract(sst,coords)

lat <- abs(coords[,2])

plot(lat,sstr, pch = 19
     , col = rgb(0,0,0,0.1))

### Define the priors
prior_fun <- list(  
  f1 = function(x,log) dsnorm(x,location = -2.66, scale = 20, alpha = 20, log = log), # prior on A (lower asymptote)
  f2 = function(x,lower,log) dtnorm(x, lower, upper = Inf, mean = 30, sd = 10, log = log), # prior on upper asymptote
  f3 = function(x,log) dnorm(x, mean = 45, sd = 12, log = log), # prior on M
  f4 = function(x,log) dlnorm(x, mean = -2.4, sd = 0.6, log = log)# prior on Q
)

# set up clusters for parallel computing
library(doParallel)
nClust <- 4
cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)

### Run model
# source script
source("R/subscripts/ClimateParallelSimple.R")
#run
system.time({modm <- climate_simple_parallel(nChains = nClust, nIter = 100000, nThin = 10,
                                x = lat, y = sstr, 
                                coeff_inits = NULL, sdy_init = NULL, 
                                prior_fun = prior_fun,
                                proposal_var_inits = c(2,2,2,0.2), adapt_sd = 2000,
                                adapt_sd_decay = 100, start_adapt = 101, quiet = FALSE)})

### stop cluster
parallel::stopCluster(cl)

### assess posterior

plot_chains(modm)
plot_gradient(modm[[1]],burnin = 1000, ylim = c(-2,30))
plot_gradient(modm[[2]],burnin = 1000, add = T, confint_col = rgb(0,0.5,1,0.2), line_col = rgb(0,0.5,1,1))
points(lat,sstr, pch = 21, bg = rgb(0,0,0,0.2), col = NA)
mod_all <- combine_posterior(modm,10000)

# mean coefficient values:
apply(mod_all,2,mean)

### We will base the prior on these, but add very large uncertainties to not make the prior 
#   dominate the analysis

# parameter A - skew normal distribution (avoiding estimating A values below freezing point of sea water)
x1 <- seq(-6,36,0.01)
dens <- dsnorm(x1,location = -3.13, scale = 12, alpha = 30,log = FALSE)
plot(x1,dens,type = "l",lwd = 2)
abline(v = c(-4,-2,0),lty = 2)
x1[which.max(dens)]
# we put the maximum density on the modern estimate for A (-1.9)

# parameter K - normal distribution
x2 <- seq(-5,65,0.01)
dens <- dnorm(x2,mean = 30, sd = 10)
plot(x2,dens,type = "l",lwd = 2)
x2[which.max(dens)]
# we put the maximum density on the modern estimate for K (30.0)

# parameter M - normal distribution
x3 <- seq(0,90,0.1)
dens <- dnorm(x3,mean = 42.1, sd = 12)
plot(x,dens,type = "l",lwd = 2)
x3[which.max(dens)]
# we put the maximum density on the modern estimate for M (42.1)

# parameter B - gamma distribution
x4 <- seq(0,0.62,0.001)
#dens <- dlnorm(x4, mean = -1.85, sd = 0.6, log = FALSE)
dens <- dgamma(x4, 3.2,20, log = FALSE)
plot(x4,dens,type = "l",lwd = 2)
#points(x4,dens,type = "l", col = rgb(.8,0,0,0.7),lwd = 2)
x4[which.max(dens)]
abline(h=0)
# we put the maximum density on the modern estimate for B (0.11)

source()

### Define the priors
prior_fun <- list(  
  f1 = function(x,log) dsnorm(x,location = -3.13, scale = 12, alpha = 30, log = log), # prior on A (lower asymptote)
  f2 = function(x,lower,log) dtnorm(x, lower, upper = Inf, mean = 30, sd = 10, log = log), # prior on upper asymptote
  f3 = function(x,log) dnorm(x, mean = 42.1, sd = 12, log = log), # prior on M
  f4 = function(x,log) dgamma(x, shape = 3.2, rate = 20, log = log)# prior on Q
)

# plot all priors
plot_prior(prior_fun,xval=list(x1,x2,x3,x4))

### Loop to run 100 time with 100 samples