# parallel chains of the runMCMC function
climate_parallel <- function(nChains = 3, nIter = 1000, obsmat = NULL, distrmat = NULL, coeff_inits = NULL, sdy_init = NULL, 
                            yest_inits = NULL, sdyest_inits = NULL,
                            proposal_var_inits = c(2,2,2,0.2), adapt_sd = floor(0.1 * nIter),
                            adapt_sd_decay = max(floor(0.005*nIter),1), quiet = FALSE) {
  foreach(pc = 1:nChains) %dopar% {
    # call model functions
    source("R/subscripts/ClimateGradientModel.R") 
    
    # set random seed
    set.seed(pc)
    
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
    if(!(is.null(distrmat)) & is.null(sdyest_inits)) sdyest_inits <- rep(2,length(unique(obsmat$sample)))
    
    run_MCMC(nIter = nIter, obsmat = obsmat, distrmat = distrmat, coeff_inits = coeff_inits, sdy_init = sdy_init, 
             yest_inits = yest_inits, sdyest_inits = sdyest_inits,
             proposal_var_inits = proposal_var_inits, adapt_sd = adapt_sd,
             adapt_sd_decay = adapt_sd_decay, quiet
    )
}
}
