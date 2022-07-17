### test ClimategradienModelwithSDonOBS

obsmat = data.frame(
  sample = c(1,1,1,2,2,2,3,3,4,4,4,5),
  latitude = c(10,10,10,11,11,11,12,12,30,30,30,60),
  temperature = c(30,33,30,40,43,40,29,30,22,24,23,8),
  sd = c(NA,NA,NA,2,4,2,1,1,1,1,NA,1))

distrmat = NULL

nIter = 1000
obsmat = obsmat
proposal_var_inits = c(2,2,2,0.2)
adapt_sd = floor(0.1 * nIter)
adapt_sd_decay = max(floor(0.01*nIter),1)
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


test2 <- run_MCMC_sd_obs(nIter = 20000, obsmat = obsmat, distrmat = NULL, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                            proposal_var_inits = c(2,2,2,0.2), adapt_sd = floor(0.1 * nIter),
                            adapt_sd_decay = max(floor(0.01*nIter),1), quiet = FALSE)
plot_gradient(test2)
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