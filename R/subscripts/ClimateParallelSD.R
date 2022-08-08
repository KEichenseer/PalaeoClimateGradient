# parallel chains of the runMCMC function
climate_parallel_sd <- function(nChains = 3, nIter = 1000, nThin = 1, obsmat = NULL, distrmat = NULL, coeff_inits = NULL, sdy_init = NULL, 
                             yest_inits = NULL, sdyest_inits = NULL, prior_fun = prior_fun,
                             proposal_var_inits = c(2,2,2,0.2), adapt_sd = floor(0.1 * nIter),
                             adapt_sd_decay = max(floor(0.01*nIter),1), start_adapt = 101, quiet = FALSE) {

  
  out <- foreach(pc = 1:nChains) %dopar% {
    # call model functions
    source("R/subscripts/ClimateGradientModelwithSDonObs.R") 
    
    # set random seed
    set.seed(pc)
    
    # random setting of initial values for the regression parameters
    coeff_inits = rep(NA,4)
    coeff_inits[1] = truncnorm::rtruncnorm(1,-5,Inf,10,5) # 20,45),c(1,2,4.5)
    coeff_inits[2] = coeff_inits[1] + truncnorm::rtruncnorm(1,0,Inf,15,7)
    coeff_inits[3] = rnorm(1,45,10)
    coeff_inits[4] = exp(rnorm(1,-2.3,0.35))
    sdy_init = exp(rnorm(1,0.7,0.25))
    
    # deterministic setting of initial values for the temperature and temperature sd estimates
    if((!(is.null(distrmat)) | !(is.null(obsmat))) & is.null(yest_inits)) {
      yest_inits <- c(unlist(sapply(unique(obsmat$sample), function(x) mean(obsmat$temperature[which(obsmat$sample == x)]))),
                      distrmat$location + distrmat$scale * sqrt(2/pi) * 
                        distrmat$shape/sqrt(1+distrmat$shape^2) )
    }
    if(!(is.null(obsmat)) & is.null(sdyest_inits)) sdyest_inits <- rep(2,length(unique(obsmat$sample)))
    
    logprior <- write_logprior(prior_fun,log = TRUE)
    
    run_MCMC_sd_obs(nIter = nIter, nThin = nThin, obsmat = obsmat, distrmat = distrmat, coeff_inits = coeff_inits, sdy_init = sdy_init, 
             yest_inits = yest_inits, sdyest_inits = sdyest_inits,
             logprior = logprior,
             proposal_var_inits = proposal_var_inits, adapt_sd = adapt_sd,
             adapt_sd_decay = adapt_sd_decay, start_adapt = 101, quiet
    )
  }
  return(out)
}
