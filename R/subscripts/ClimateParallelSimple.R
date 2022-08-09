# Wrapper function for parallel chains 
# for ClimateGradientModelSimple.R
climate_simple_parallel <- function(nChains = 3, nIter = 1000, nThin = 1, x = NULL, y = NULL, 
                                    coeff_inits = NULL, sdy_init = NULL,
                                    proposal_var_inits = c(2,2,2,0.2),
                                    prior_fun = NULL,
                                    adapt_sd = floor(0.1 * nIter),
                                    adapt_sd_decay = max(floor(0.01*nIter),1),
                                    start_adapt = min(101,nIter),
                                    quiet = FALSE
                                    ) {
  

  
  foreach(pc = 1:nChains) %dopar% {
    # call model functions
    source("R/subscripts/ClimateGradientModelSimple.R") 
    #make_prior(priorvec)
    
    # set random seed
    set.seed(pc)
    
    # random setting of initial values for the regression parameters
    coeff_inits = rep(NA,4)
    coeff_inits[1] = rnorm(1,10,3) # 20,45),c(1,2,4.5)
    coeff_inits[2] = coeff_inits[1] + truncnorm::rtruncnorm(1,0,Inf,12,6)
    coeff_inits[3] = rnorm(1,45,7.5)
    coeff_inits[4] = exp(rnorm(1,-2.3,0.25))
    sdy_init = exp(rnorm(1,0.7,0.25))
    
    logprior <- write_logprior(prior_fun,log = TRUE)
    
    run_MCMC_simple(x = x, y = y, nIter = nIter, nThin = nThin,
                    coeff_inits = coeff_inits, sdy_init = sdy_init,
                    proposal_var_inits = proposal_var_inits,
                    logprior = logprior,
                    adapt_sd = adapt_sd, 
                    adapt_sd_decay = adapt_sd_decay,
                    start_adapt = start_adapt,
                    quiet = TRUE
                    
    )
  }
}
