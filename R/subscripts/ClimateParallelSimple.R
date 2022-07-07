# Wrapper function for parallel chains 
# for ClimateGradientModelSimple.R
climate_simple_parallel <- function(nChains = 3, nIter = 1000, latitude = NULL, temperature = NULL, 
                                    priorvec = NULL,
                                    coeff_inits = NULL, sdy_init = NULL) {
  foreach(pc = 1:nChains) %dopar% {
    # call model functions
    source("R/subscripts/ClimateGradientModelSimple.R") 
    make_prior(priorvec)
    
    # set random seed
    set.seed(pc)
    
    # random setting of initial values for the regression parameters
    coeff_inits = rep(NA,4)
    coeff_inits[1] = rnorm(1,10,3) # 20,45),c(1,2,4.5)
    coeff_inits[2] = coeff_inits[1] + truncnorm::rtruncnorm(1,0,Inf,12,6)
    coeff_inits[3] = rnorm(1,45,7.5)
    coeff_inits[4] = exp(rnorm(1,-2.3,0.25))
    sdy_init = exp(rnorm(1,0.7,0.25))
    
    
    run_MCMC_simple(x = latitude, y = temperature, nIter = nIter,
                    coeff_inits = coeff_inits, sdy_init = sdy_init
    )
  }
}
