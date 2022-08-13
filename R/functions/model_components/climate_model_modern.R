######################################################
######################################################
### A Metropolis-Hastings algorithm for Latitudinal
### Temperature Gradients - Simple version
###
### first developed: February 2022
### revised: June - August 2022
######################################################
######################################################

# Main MCMCM function
run_MCMC_simple <- function(x, y, n_iter, n_thin = 1,
                            coeff_inits = NULL, sdy_init = NULL,
                            proposal_var_inits = NULL,
                            prior_input = NULL,
                            adapt_sd = NULL, 
                            adapt_sd_decay = NULL,
                            start_adapt = NULL,
                            quiet = FALSE){
  
  ### Load functions
  source("R/functions/model_components/dsnorm.R")  
  source("R/functions/model_components/dtnorm.R")  
  source("R/functions/model_components/gradient.R")  
  source("R/functions/model_components/write_logprior.R")  
  source("R/functions/model_components/loglik.R")  
  source("R/functions/model_components/logposterior.R")  
  source("R/functions/model_components/MH_propose_multi.R")  
  source("R/functions/model_components/weighted_cov.R")  
  
  
  ### Initialisation
  
  # random setting of initial values for the regression parameters
  if(is.null(coeff_inits)) {coeff_inits = rep(NA,4)
  coeff_inits[1] = rnorm(1,10,3) # 20,45),c(1,2,4.5)
  coeff_inits[2] = coeff_inits[1] + truncnorm::rtruncnorm(1,0,Inf,12,6)
  coeff_inits[3] = rnorm(1,45,7.5)
  coeff_inits[4] = exp(rnorm(1,-2.3,0.25))
  }
  if(is.null(sdy_init)) sdy_init = exp(rnorm(1,0.7,0.25))
  
  
  logprior <- write_logprior(prior_fun = prior_input, log = TRUE) # create prior
  
  if(is.null(proposal_var_inits)) proposal_var_inits <- c(2^2,2^2,2^2,0.2^2)
  if(is.null(adapt_sd)) adapt_sd <- floor(0.1 * n_iter)
  if(is.null(adapt_sd_decay)) adapt_sd_decay <- min(c(150,max(floor(0.1*adapt_sd),1))) # cap decay at 150 
  if(is.null(start_adapt)) start_adapt <- min(c(100,max(floor(0.05*adapt_sd),1))) # cap start adapt at 100
  
  
  save_it <- seq(1,n_iter,n_thin) # iterations to save
  
  coefficients = array(dim = c(n_iter,4)) # set up array to store coefficients
  coefficients[1,] = coeff_inits # initialise coefficients
  sdy = rep(NA_real_,n_iter) # set up vector to store sdy
  sdy[1] = sdy_init # intialise sdy
  A_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  n <- length(y)
  shape_sdy <- A_sdy+n/2 # shape parameter for the inverse gamma
  
  # for proposals / adaptation
  proposal_cov <- matrix(0,4,4)
  diag(proposal_cov) <- proposal_var_inits
  proposal_var_inits <- proposal_cov
  
  if(is.numeric(adapt_sd)){
    if(adapt_sd<10) stop("adapt_sd needs to be >=10")
    all_weights <- exp((-(adapt_sd-1)):0/adapt_sd_decay)
    adapt_it <- seq(start_adapt,adapt_sd,10) # adapt covariance only at every 10th iteration
  }
  
  proposal_factor <- 1 # to adjust acceptance rate
  
  accept = rep(NA,n_iter)
  
  # setup progress bar
  if (!quiet) cli::cli_progress_bar('Sampling', total = n_iter)
  
  ### The MCMC loop
  for (i in 2:n_iter){
    
    # update progress bar
    if (!quiet) cli::cli_progress_update(set = i, status = paste0('iteration ', i))
    
    ## 1. Gibbs step to estimate sdy
    sdy[i] = sqrt(1/rgamma(
      1,shape_sdy,B_sdy+0.5*sum((y-gradient(x,coefficients[i-1,],0))^2)))
    
    ## 2. Metropolis-Hastings step to estimate the regression coefficients
    # create matrix of proposal innovations as this is much faster than doing it anew at every it
    if(i == adapt_sd+1 | (i==2 & adapt_sd < 2)) proposal_innovation <-   mvnfast::rmvn(
      n = n_iter-adapt_sd, mu = rep(0,4),sigma = 2.4/sqrt(4)*proposal_cov)+
        rnorm(n = 4*(n_iter-adapt_sd), mean = rep(0,4),sd = 0.001)
    # create proposals
    if(i <= adapt_sd) proposal_coeff = MH_propose_multi(1,coefficients[i-1,],proposal_cov =  proposal_cov) # new proposed values
    if(i > adapt_sd) proposal_coeff = c(coefficients[i-1,1:3],log(coefficients[i-1,4])) + proposal_innovation[i-adapt_sd,]
    proposal_coeff[4] <- exp(proposal_coeff[4])
    
    if(i <= adapt_sd) proposal_coeff = MH_propose_multi(1,coefficients[i-1,],proposal_cov =  proposal_cov) # new proposed values
    if(i > adapt_sd) proposal_coeff = c(coefficients[i-1,1:3],log(coefficients[i-1,4])) + proposal_innovation[i-adapt_sd,]
    proposal_coeff[4] <- exp(proposal_coeff[4])
    
    #if(any(proposal[4] <= 0)) HR = 0 else # B needs to be >0
    # Hastings ratio of the proposal
    HR = exp(logposterior(x = x, y = y, coeff = proposal_coeff, sdy = sdy[i], logprior = logprior) -
               logposterior(x = x, y = y, coeff = coefficients[i-1,], sdy = sdy[i], logprior = logprior) +
               (-log(coefficients[i-1,4])) -
               (-log(proposal_coeff[4])))
    # accept proposal with probability = min(HR,1)
    if (runif(1) < HR){
      accept[i] = 1
      coefficients[i,] = proposal_coeff
      # if proposal is rejected, keep the values from the previous iteration
    }else{
      accept[i] = 0
      coefficients[i,] = coefficients[i-1,]

    }
    
    ###
    ###
    ### Adaptation step
    if(i <= adapt_sd && (i>=start_adapt & i %in% adapt_it)) {
      weights = all_weights[(adapt_sd-i+1):adapt_sd]
      proposal_cov <- weighted_cov(cbind(coefficients[1:i,1:3],log(coefficients[1:i,4])),weights = weights)
      if(any(diag(proposal_cov)==0)) proposal_cov <- proposal_var_inits/i
      #if(i %in% seq(2*floor(adapt_sd/10), adapt_sd, floor(adapt_sd/10))) {
      #
      #  if(mean(accept[(i-floor(adapt_sd/10)):i]) < 0.23) proposal_factor <- proposal_factor - 0.1*proposal_factor
      #  if(mean(accept[(i-floor(adapt_sd/10)):i]) > 0.45) proposal_factor <- proposal_factor + 0.1*proposal_factor
      # proposal_cov <- proposal_cov * proposal_factor
      #  while(any(eigen(proposal_cov)$values <= 0.000001)) {
      #   diag(proposal_cov) <- 1.25 * diag(proposal_cov)
      #    # print(paste("stuck in while loop at it",i))}
      # }
      #
      #}
    }
    
    # # create matrix of proposal innovations as this is much faster than doing it anew at every it
    # if(i == adapt_sd) proposal_innovation <-   mvnfast::rmvn(n = n_iter-adapt_sd, mu = rep(0,4),
    #                                                          sigma = 2.4/sqrt(4)*proposal_cov)+
    #   rnorm(n = 4*(n_iter-adapt_sd), mean = rep(0,4),
    #         sd = 0.001)
    
    
  } # end of the MCMC loop
  
  ###  Function output
  output = list(params = data.frame(A = coefficients[save_it,1],
                      dKA = coefficients[save_it,2],
                      M = coefficients[save_it,3],
                      B = coefficients[save_it,4],
                      sdy = sdy[save_it]))
  output$call = mget(names(formals()),sys.frame(sys.nframe()))

  return(output)
}
