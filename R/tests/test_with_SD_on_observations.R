### test ClimategradientModelwithSDonOBS
source("R/subscripts/ClimateGradientModel.R")

source("R/subscripts/AuxiliaryFunctions.R")

obsmat = data.frame(
  sample = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,6,7,7),
  latitude = c(9,9,9,11,11,11,13,13,13,30,30,30,60,60,40,70,70),
  temperature = c(30,33,30,36,39,36,29,30,29.5,22,24,23,8,9,19,6,5),
  sd = c(NA,NA,NA,2,3,2,1,1,1,1,1,NA,1,1,2,1,NA))

distrmat = NULL

nIter = 20000
obsmat = obsmat
proposal_var_inits = c(2,2,2,0.2)
adapt_sd = floor(0.1 * nIter)
adapt_sd_decay = max(floor(0.01*nIter),1)
start_adapt = 101
quiet = FALSE
yest_inits = NULL
sdyest_inits = NULL
# set random seed
set.seed(1)

# random setting of initial values for the regression parameters
coeff_inits = rep(NA,4)
coeff_inits[1] = rnorm(1,10,3) # 20,45),c(1,2,4.5)
coeff_inits[2] = coeff_inits[1] + truncnorm::rtruncnorm(1,0,Inf,12,6)
coeff_inits[3] = rnorm(1,45,7.5)
coeff_inits[4] = exp(rnorm(1,-2.3,0.25))
sdy_init = exp(rnorm(1,0.7,0.25))

# deterministic setting of initial values for the temperature and temperature sd estimates
if((!(is.null(distrmat)) | !(is.null(obsmat))) & is.null(yest_inits)) {
  yest_inits <- c(unlist(sapply(unique(obsmat$sample), function(x) mean(obsmat$temperature[which(obsmat$sample == x)]))),
                  distrmat$location + distrmat$scale * sqrt(2/pi) * 
                    distrmat$shape/sqrt(1+distrmat$shape^2) )
}
if(!(is.null(obsmat)) & is.null(sdyest_inits)) sdyest_inits <- rep(2,length(unique(obsmat$sample)))

source("R/subscripts/ClimateGradientModelwithSDonObs.R")
test4.4sd <- run_MCMC_sd_obs(nIter = 20000, obsmat = obsmat, distrmat = NULL, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                               proposal_var_inits = c(2,2,2,0.2), adapt_sd = floor(0.2 * nIter),
                               adapt_sd_decay = max(floor(0.005*nIter),1),start_adapt = start_adapt, quiet = FALSE)


source("R/subscripts/ClimateGradientModel.R")
test1.2 <- run_MCMC(nIter = 20000, obsmat = obsmat, distrmat = NULL, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                         proposal_var_inits = c(2,2,2,0.2), adapt_sd = floor(0.2 * nIter),
                         adapt_sd_decay = max(floor(0.005*nIter),1), quiet = FALSE)

plot(test4.4sd$params$A[1:400],col=rgb(1,0,0,0.4))
points(test4.3nosd$params$A,col=rgb(1,0,1,0.4))

plot_gradient(test4.3nosd, ylim = c(-1,43), add = F, burnin = 5000, line_col = rgb(0,0,0.75,.7),confint_col = rgb(0,0,0.75,0.2))
plot_gradient(test4.4sd, ylim = c(-1,43), add = F, burnin = 5000, line_col = rgb(0.75,0,0,.7),confint_col = rgb(0.75,0,0,0.2))
plot_gradient(test4.3sdlf, ylim = c(-1,43), add = T, burnin = 5000, line_col = rgb(0,0.75,0,.7),confint_col = rgb(0,0.75,0,0.2))



#plot_gradient(test7, ylim = c(-1,43), add = T, burnin = 22000, line_col = rgb(0.75,.5,0,.7),confint_col = rgb(0.75,.5,0,0.2))
plot_gradient(test3, ylim = c(-1,43), add = T, burnin = 25000, line_col = rgb(0,.75,0,.7),confint_col = rgb(0,.75,0,0.2))
plot_gradient(test4, ylim = c(-1,43), add = T, burnin = 25000, line_col = rgb(0.75,0,0,.7),confint_col = rgb(0.75,0,0,0.2))
plot_gradient(test13, ylim = c(-1,43), add = T, burnin = 15000, line_col = rgb(0,0,0.75,.7),confint_col = rgb(0,0,0.75,0.2))

plot_sample_gradient(test8, line_col  = rgb(0,0.7,0.3,0.33))
plot_sample_gradient(test3, line_col  = rgb(0,0.3,0.7,0.33))

plot_gradient(test3, ylim = c(-1,43), add = T, line_col = rgb(0,.7,1,.7),confint_col = rgb(0,.7,1,0.2))
plot_posterior(test3, col_obs = rgb(0,0.3,1,0.75))
plot_posterior(test4.3sdlf, col_obs = rgb(0,1,0,0.75))
plot_posterior(test4.4sd, col_obs = rgb(1,0,0,0.75))
plot_posterior(test4.3nosd, col_obs = rgb(0,0,1,0.75))

plot(test9$params$A[21000:30000], test9$params$DKA[21000:30000], type = "o", col = rgb(0,0,0,0.2))


plot(test4.3nosd$params$A,test4.3nosd$params$DKA,type = "l", col = rgb(0,0,0,0.2))
points(test4.2$params$A,test4.2$params$DKA,type = "l", col = rgb(1,0,0,0.2))

points(obsmat$latitude,obsmat$temperature, pch = 19, col = rgb(0,0,0,0.33))
plot(test$yestimate[,3])
plot(test$obs_yestimate[[2]][,1])

logposterior_norm_sd(x = x[n_p_ind], yest = yestimate[i,n_p_ind], ymean = yobs_mean_old,
                     sdyest = sdyest[i,], coeff = coefficients[i-1,],
                     sdy = sdy[i],new_y =new_y_old)

loglik_norm_sd <- function(x, yest, ymean, sdyest, coeff, sdy, new_y) {
  # extract regression coefficients
  coeff = unlist(coeff)
  A = coeff[1]
  DKA = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  yest0 <- new_y[sd_obs_ind]
  ymean0 <- obsmat$temperature[sd_obs_ind]
  ysd0 <- obsmat$sd[sd_obs_ind]
  
  ll0 <- sum(dnorm(yest0, ymean0, ysd0,log=TRUE),na.rm=T)
  
  ll1 <- sum(dnorm(yest, ymean, sdyest,log=TRUE),na.rm=T)
  
  pred = A + DKA/((1+(exp(Q*(x-M)))))
  ll2 <- sum(dnorm(yest, mean = pred, sd = sdy, log = TRUE))
  return(c(ll2+ll1+ll0))
}