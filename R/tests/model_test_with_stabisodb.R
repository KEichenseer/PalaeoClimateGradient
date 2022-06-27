### A test with StabisoDB data
### 26/06/2022
###

# load data
dat <- readRDS("data/processed/StabisoDB_processed_26_06_2022.rds")

# select data from one stage to test, exclude NA data
dat_sub <- subset(dat,stage_2020 == "Changhsingian" & !(is.na(paleolat)) & !(is.na(temperature)))#table(dat$stage_2020)
dat_sub <- dat_sub[with(dat_sub, order(abs(paleolat), longitude)),]
# prepare for use in the model
dat_mod <- data.frame(sample = (paste(abs(dat_sub$paleolat),dat_sub$longitude)),
                       latitude = abs(dat_sub$paleolat), temperature = dat_sub$temperature)
samples1 <- unique(dat_mod$sample)
for(i in 1:length(samples1)) dat_mod$sample[which(dat_mod$sample==samples1[i])] <- i

samples <-  1:length(samples1)
sample_lats <- sapply(samples, function(x) unique(dat_mod$lat[which(dat_mod$sample==x)]))
  
# visualise data
plot(abs(dat_mod$latitude),dat_mod$temperature)

# wait a second - how do we deal with situationsn with 1 value per locality?
## ok, this is taken care of by the prior on yest ... although it is counterintuitive, as yest deviates less 
## with larger priors on sdyest. Definetly some prior sensitivity tests are in order. 
### resolved, it is an inverse gamma distribution, hence larger scale leads to smaller deviations, all good.

# call model functions
source("R/subscripts/ClimateGradientModel.R")

# initial temperature values for the sample locations
yest_inits <- sapply(unique(dat_mod$sample), function(x) mean(dat_mod$temperature[which(dat_mod$sample == x)]))

# run the model
mod <- run_MCMC(nIter = 20000, obsmat = dat_mod, distrmat = NULL, coeff_inits = c(10,20,45,0.1), sdy_init = 1, 
         yest_inits = yest_inits, sdyest_inits = rep(2,length(yest_inits)),
                     proposal_var_inits = c(2,2,2,0.2))


nIter <- 20000
burnin <- 4000+1 # burnin - to be discarded for the plots

# function for plotting the 95 % CI shading
source("R/subscripts/AuxiliaryFunctions.R")

## visualise data and output
plot(dat_mod$latitude,dat_mod$temperature, pch = 19, col = rgb(0,0,0,0.33), xaxs = "i", yaxs="i",
     xlim = c(0,90), ylim = c(-3,45),
     xlab = expression ("absolute latitude ("*degree*")"), ylab = expression("temperature ("*degree~"C)"))


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


points(sample_data$x, sample_data$y,  pch = 19, cex = 1.1, col = rgb(0,0,0,0.5), xpd = T)
replicate(8, points(latitude, gradient(x = latitude, coeff = unlist(mod$params[sample(burnin:nIter,1),1:4]), sdy = 0), 
                    type = "l", col = rgb(0,0.25,0.5,0.33), lwd = 2))
#axis(2,seq(-5,30,5),c(NA,0,NA,10,NA,20,NA,30))

### Ghzelian: Investigate why the lat 20 point estimate is so offset
### solved: was a problem with assigning sample numbers. 