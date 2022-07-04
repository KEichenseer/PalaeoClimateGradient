
### A test with StabisoDB and coral reef data
### 02/07/2022
###

mystage <- "Turonian"

# load isotope data
iso <- readRDS("data/processed/StabisoDB_processed_26_06_2022.rds")

# select data from one stage to test, exclude NA data
iso_sub <- subset(iso,stage_2020 == mystage & !(is.na(paleolat)) & !(is.na(temperature)))#table(iso$stage_2020)
iso_sub <- iso_sub[with(iso_sub, order(abs(paleolat), longitude)),]
# prepare for use in the model
iso_mod <- data.frame(sample = (paste(abs(iso_sub$paleolat),iso_sub$longitude)),
                      latitude = abs(iso_sub$paleolat), temperature = iso_sub$temperature)
samples1 <- unique(iso_mod$sample)
for(i in 1:length(samples1)) iso_mod$sample[which(iso_mod$sample==samples1[i])] <- i

nsamples <- length(samples1)
samples <-  1:(max(1,nsamples))
sample_lats <- sapply(samples, function(x) unique(iso_mod$lat[which(iso_mod$sample==x)]))


# load coral reef data
coral <- readRDS("data/processed/PARED_coralreefs_processed_02_07_2022.rds")

coral_sub <- subset(coral,early_stage == mystage & !(is.na(pal_lat_scotese)))#table(iso$stage_2020)
coral_sub <- coral_sub[with(coral_sub, order(abs(pal_lat_scotese), longit)),]
# prepare for use in the model
if(nrow(coral_sub) >= 1) {coral_distrmat <- data.frame(latitude = abs(coral_sub$pal_lat_scotese),
                                                       location = 22.8,
                                                       scale = 10,
                                                       shape = 4,
                                                       distribution = "skew-normal")
} else coral_distrmat <- data.frame(NULL)

# visualise data
plot(abs(iso_mod$latitude),iso_mod$temperature, xlim = range(c(iso_mod$latitude, coral_distrmat$latitude)),
     ylim = c(0,50))

z1 <- list(NULL)
y1 <- list(NULL)
yskew_rho <-  as.numeric(coral_distrmat[,"shape"])/sqrt(1+as.numeric(coral_distrmat[,"shape"])^2)
N = 5000
i1 <- 0
for(d in 1:nrow(coral_distrmat)) {
  i1=i1+1
  z1[[i1]] <- truncnorm::rtruncnorm(n = N, a = 0, b = Inf, mean = 0, sd = 1)
  #y1 <- epsilon1 + a1*z1 + rnorm(N,0,omega1)
  y1[[i1]] <- as.numeric(coral_distrmat[d,"location"]) + as.numeric(coral_distrmat[d,"scale"])*yskew_rho[d]*z1[[i1]] +
    as.numeric(coral_distrmat[d,"scale"])*sqrt(1-(yskew_rho[d]^2))*rnorm(N,0,1)
  
  beanplot::beanplot(y1[[i1]], add = T,
                     at = as.numeric(coral_distrmat[d,"latitude"]), maxwidth = 0.5, side = "second",
                     what = c(0,1,0,0), col = rgb(0,0,0,0.2), border = NA, xpd = T)
}

# wait a second - how do we deal with situationsn with 1 value per locality?
## ok, this is taken care of by the prior on yest ... although it is counterintuitive, as yest deviates less 
## with larger priors on sdyest. Definetly some prior sensitivity tests are in order. 
### resolved, it is an inverse gamma distribution, hence larger scale leads to smaller deviations, all good.

# call model functions
source("R/subscripts/ClimateGradientModel.R")

# initial temperature values for the sample locations
yest_inits <- c(unlist(sapply(unique(iso_mod$sample), function(x) mean(iso_mod$temperature[which(iso_mod$sample == x)]))),
                       coral_distrmat$location + coral_distrmat$scale * sqrt(2/pi) * 
                         coral_distrmat$shape/sqrt(1+coral_distrmat$shape^2) )
                
nIter <- 20000


# run the model
mod <- run_MCMC(nIter = nIter, obsmat = iso_mod, distrmat = coral_distrmat, coeff_inits = c(10,20,45,0.1), sdy_init = 1, 
                yest_inits = yest_inits, sdyest_inits = rep(2,length(unique(iso_mod$sample))),
                proposal_var_inits = c(2,2,2,0.2))


burnin <- 4000+1 # burnin - to be discarded for the plots

# function for plotting the 95 % CI shading
source("R/subscripts/AuxiliaryFunctions.R")

## visualise data and output
plot(iso_mod$latitude,iso_mod$temperature, pch = 19, col = rgb(0,0,0,0.33), xaxs = "i", yaxs="i",
     xlim = c(0,90), ylim = c(-3,45),
     xlab = expression ("absolute latitude ("*degree*")"), ylab = expression("temperature ("*degree~"C)"))


z1 <- list(NULL)
y1 <- list(NULL)
yskew_rho <-  as.numeric(coral_distrmat[,"shape"])/sqrt(1+as.numeric(coral_distrmat[,"shape"])^2)
N = 5000
i1 <- 0
for(d in 1:nrow(coral_distrmat)) {
  i1=i1+1
  z1[[i1]] <- truncnorm::rtruncnorm(n = N, a = 0, b = Inf, mean = 0, sd = 1)
  #y1 <- epsilon1 + a1*z1 + rnorm(N,0,omega1)
  y1[[i1]] <- as.numeric(coral_distrmat[d,"location"]) + as.numeric(coral_distrmat[d,"scale"])*yskew_rho[d]*z1[[i1]] +
    as.numeric(coral_distrmat[d,"scale"])*sqrt(1-(yskew_rho[d]^2))*rnorm(N,0,1)
  
  beanplot::beanplot(y1[[i1]], add = T,
                     at = as.numeric(coral_distrmat[d,"latitude"]), maxwidth = 5, side = "second",
                     what = c(0,1,0,0), col = rgb(0,0,0,0.2), border = NA, xpd = T)
}

latitude <- seq(0,90,by=0.5)
temperature <- gradient(x = latitude, coeff = apply(mod$params[burnin:nIter,1:4],2,median), sdy = 0)

points(latitude, temperature, type = "l", lwd = 2)
sample_it <- sample(burnin:nIter,1000)
grad_025 <- sapply(1:length(latitude), function(f) quantile(apply(mod$params[sample_it,1:4],1,function(a) gradient(
  x=latitude[f], coeff =  a, sdy = 0)), probs = 0.025))
grad_975 <- sapply(1:length(latitude), function(f) quantile(apply(mod$params[sample_it,1:4],1,function(a) gradient(
  x=latitude[f], coeff =  a, sdy = 0)), probs = 0.975))

error_polygon(latitude,grad_025,grad_975,rgb(0,0,0,0.25))


sapply(samples,function(x) points(sample_lats[x], mean(mod$yestimate[burnin:nIter,x]), pch = 17, col = rgb(1,0,0,0.75)))
sapply(samples,function(x) points(rep(sample_lats[x],2), quantile(mod$yestimate[burnin:nIter,x], probs = c(0.05,0.95)), 
                                  type = "l", col = rgb(1,0,0,0.75)))

distr_ind <- (nsamples+1):(nsamples+nrow(coral_distrmat))
sapply(1:length(distr_ind),function(x) points(coral_distrmat$latitude[x], median(mod$yestimate[burnin:nIter,distr_ind[x]]), pch = 17, col = rgb(0,0.5,1,0.75)))
sapply(1:length(distr_ind),function(x) points(rep(coral_distrmat$latitude[x],2), quantile(mod$yestimate[burnin:nIter,distr_ind[x]], probs = c(0.05,0.95)), 
                                  type = "l", col = rgb(0,0.5,1,0.75)))


replicate(8, points(latitude, gradient(x = latitude, coeff = unlist(mod$params[sample(burnin:nIter,1),1:4]), sdy = 0), 
                    type = "l", col = rgb(0,0.25,0.5,0.33), lwd = 2))
#axis(2,seq(-5,30,5),c(NA,0,NA,10,NA,20,NA,30))

### Ghzelian: Investigate why the lat 20 point estimate is so offset
### solved: was a problem with assigning sample numbers. 

############################################
####
#### Test climate_parallel
####
source("R/subscripts/ClimateParallel.R")

library(foreach)
library(doParallel)

cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)

mod3 <- climate_parallel(nChains = 3, nIter = nIter, obsmat = iso_mod, distrmat = coral_distrmat)
stopImplicitCluster()

plot(mod3[[1]]$params$DKA[seq(1,20000,10)],type = "l", col = rgb(0.2,0.8,0,0.5))
points(mod3[[2]]$params$DKA[seq(1,20000,10)],type = "l", col = rgb(.8,0,.8,0.5))
points(mod3[[3]]$params$DKA[seq(1,20000,10)],type = "l", col = rgb(0,0.6,.8,0.5))

### Effective sample size 
mcmcse::multiESS(mod3$params[(0.5*nIter):nIter,1:4])

### tests plot function
plot_gradient(mod, ylim = c(-5,35))
plot_data(iso_mod,add = T)
plot_distr(coral_distrmat)
plot_posterior(mod)

