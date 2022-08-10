### load climate data
# Mean annual sea surface temperatures (Bio-Oracle)
set.seed(1)

sst <- raster::raster("./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")
r <- raster::raster(res = 1)
sst <- raster::resample(x = sst, y = r)

sst <- raster::as.data.frame(x = sst, xy = TRUE, centroids = TRUE)
sst <- na.omit(sst)

# sample points randomly from the raster (within specified latitudinal limits) 
source("R/subscripts/AuxiliaryFunctions.R")

temp <- sst$Present.Surface.Temperature.Mean

lat <- abs(sst$y)

plot(lat,temp, pch = 21,col=NA, bg = rgb(0,0,0,0.1))

### Define the priors
prior_fun <- list(  
  f1 = function(x,log) dsnorm(x,location = -3.03, scale = 12, alpha = 30, log = log), # prior on A (lower asymptote)
  f2 = function(x,lower,log) dtnorm(x, lower, upper = Inf, mean = 28.3, sd = 10, log = log), # prior on upper asymptote
  f3 = function(x,log) dnorm(x, mean = 42.1, sd = 12, log = log), # prior on M
  f4 = function(x,log) dgamma(x, shape = 3.2, rate = 20, log = log)# prior on Q
)

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
  source("R/subscripts/models/modern_climate_model.R")
  # set random seed
  set.seed(pc)
  run_MCMC_simple(n_iter = 100000, n_thin = 10,
                               x = lat, y = temp, 
                               coeff_inits = NULL, sdy_init = NULL, 
                               logprior_input = prior_fun,
                               proposal_var_inits = c(2,2,2,0.2), adapt_sd = 2500,
                               adapt_sd_decay = 100, start_adapt = 101, quiet = FALSE)
}
                
### stop cluster
parallel::stopCluster(cl)

### assess posterior

plot_chains(modm)
plot_gradient(modm[[1]],burnin = 1000, ylim = c(-2,30))
plot_gradient(modm[[2]],burnin = 1000, add = T, confint_col = rgb(0,0.5,1,0.2), line_col = rgb(0,0.5,1,1))
points(lat,ssts, pch = 21, bg = rgb(0,0,0,0.2), col = NA)

# discard 10k iterations as burnin each, and combine the results of the 4 chains
mod_all <- combine_posterior(modm,burnin = 10000)
mcmcse::multiESS(mod_all[,1:5])
saveRDS(mod_all,"results/modern/modern_sample_gradient.RDS")

lat <- seq(0,90,0.1)
med_grad <- gradient(lat,apply(mod_all[,1:4],2,median),0)
plot(lat,med_grad)

plot_sample_gradient(modm[[1]])

# mean coefficient values:
apply(mod_all,2,median)

### We will base the prior on these, but add very large uncertainties to not make the prior 
#   dominate the analysis

# parameter A - skew normal distribution (avoiding estimating A values below freezing point of sea water)
x1 <- seq(-6,36,0.01)
dens <- dsnorm(x1,location = -3.03, scale = 12, alpha = 30,log = FALSE)
plot(x1,dens,type = "l",lwd = 2)
abline(v = c(-4,-2,0),lty = 2)
x1[which.max(dens)]
# we put the maximum density on the modern estimate for A (-1.7)

# parameter K - normal distribution
x2 <- seq(-5,65,0.01)
dens <- dnorm(x2,mean = 28.3, sd = 10)
plot(x2,dens,type = "l",lwd = 2)
x2[which.max(dens)]
# we put the maximum density on the modern estimate for K+A (28.3)

# parameter M - normal distribution
x3 <- seq(0,90,0.1)
dens <- dnorm(x3,mean = 42.1, sd = 12)
plot(x3,dens,type = "l",lwd = 2)
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


### Run 100 analyses with 100 modern samples, based on the Eocene palaeolatitudes


### Define the priors
prior_fun <- list(  
  f1 = function(x,log) dsnorm(x,location = -3.03, scale = 12, alpha = 30, log = log), # prior on A (lower asymptote)
  f2 = function(x,lower,log) dtnorm(x, lower, upper = Inf, mean = 28.3, sd = 10, log = log), # prior on upper asymptote
  f3 = function(x,log) dnorm(x, mean = 42.1, sd = 12, log = log), # prior on M
  f4 = function(x,log) dgamma(x, shape = 3.2, rate = 20, log = log)# prior on Q
)

saveRDS(prior_fun,"results/modern/prior.RDS")

# plot all priors
xlist = list(seq(-6,36,0.01),
             seq(-5,65,0.01),
             seq(0,90,0.1),
             seq(0,0.62,0.001)
             )
plot_prior(prior_fun,xval=xlist)


# read modern samples
modt_samples <- readRDS("results/modern/modern_sample.RDS")

### Run model
# source script
source("R/subscripts/models/modern_climate_model.R")


mods <- list(NULL)

set.seed(1)

system.time({
for(s in 1:100) {
  
  print(s)
  p_lat = abs(modt_samples[[s]][,1])
  temp = modt_samples[[s]][,2]

  coeff_inits <- c(rnorm(3,c(5,30,45), c(3,5,7)),exp(rnorm(1,log(0.1),0.6)))
  sdy_init <- exp(rnorm(1,log(2),0.6))
  mods[[s]] <- run_MCMC_simple(n_iter = 25000, n_thin = 10,
                                             x = p_lat, y = temp, 
                                             coeff_inits = coeff_inits, sdy_init = sdy_init, 
                                             logprior_input  = prior_fun,
                                            proposal_var_inits = c(3,3,3,0.2), adapt_sd = 2500,
                                             adapt_sd_decay = 100, start_adapt = 101, quiet = FALSE)
}
})

saveRDS(mods,"results/modern/modern_sample_eocene_p_lat_gradient.RDS")

# assess posterior of some chains 
plot_chains(mods[1:4])
mods_all <- combine_posterior(mods,5000)
plot_gradient(mods_all)
#saveRDS(mods_all,"results/modern/combined_100_modern_gradients_with_modern_T_and_Eocene_palaeolats.RDS")

mcmcse::multiESS(mods[[100]]$params[,1:4])
n_iter = 25000
burnin = 5000
n_thin = 10
# plot results
latx <- seq(0,90,0.1)

med_grad <- array(NA_real_,dim = c(100,length(latx)))
for(i in 1:100) {
  print(i)
  grads <- gradient(latx,mods[[i]]$params[(burnin/n_thin):(n_iter/n_thin),1:4],0)
  med_grad[i,] <- apply(grads,2,median)
}
### Plot the gradients
ngrad <- 100
s <- 2
plot(0,0,type = "n", xlim = c(0,90), ylim = c(-5,35), xlab = "|latitude|", ylab = "temperature (deg C)")
#points(abs(coords[,2]),ssts, pch = 21, bg = rgb(0,0,0,0.1), col = NA)

for(i in 1:100) points(latx,med_grad[i,],type = "l", col = rgb(i/ngrad,0.25,((ngrad-i+1)/ngrad),0.5),lwd = 2)
points(modt_samples[[i]][,1],modt_samples[[i]][,2],pch = 21, bg = rgb(0,0,0.8,0.33), col = NA)

plot_gradient(modm[[1]],burnin = 1000, add = T, confint_col = rgb(0,0,0,0.2), line_col = rgb(0,0,0,1), lwd = 4)

mean(modm[[1]]$params$sdy)

mods_all <- combine_posterior(mods,5000)
plot_gradient(mods_all,add=F)
