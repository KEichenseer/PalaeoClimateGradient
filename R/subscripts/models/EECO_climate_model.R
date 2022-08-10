############## Climate gradient model
############## 26/06/2022 # updated July 2022 to add SD on observations

### !! Careful still need to make the obsmat sd checks conditional on obsmat not being NULL !!







#
# Main MCMCM function
hierarchical_model <- function(n_iter = 1000, n_thin = 1, obsmat = NULL, 
                               distrmat = NULL, logprior_input = NULL,
                               coeff_inits = NULL, sdy_init = NULL,
                               yest_inits = NULL, sdyest_inits = NULL,
                     proposal_var_inits = NULL, adapt_sd = NULL,
                     adapt_sd_decay = NULL,  start_adapt = NULL, 
                     A_sdy = 10, B_sdy = 10,
                     quiet = FALSE){
  
  ### Load functions
  source("R/functions/model_components/dsnorm.R")  
  source("R/functions/model_components/dtnorm.R")  
  source("R/functions/model_components/gradient.R")  
  source("R/functions/model_components/write_logprior.R")  
  source("R/functions/model_components/loglik.R")  
  source("R/functions/model_components/MH_propose_multi.R")  
  source("R/functions/model_components/weighted_cov.R")  
  source("R/functions/model_components/skew_mu.R")  
  
  ### Initialisation
  
  # random setting of initial values for the regression parameters
  if(is.null(coeff_inits)) {coeff_inits = rep(NA,4)
  coeff_inits[1] = rnorm(1,10,3) # 20,45),c(1,2,4.5)
  coeff_inits[2] = coeff_inits[1] + truncnorm::rtruncnorm(1,0,Inf,12,6)
  coeff_inits[3] = rnorm(1,45,7.5)
  coeff_inits[4] = exp(rnorm(1,-2.3,0.25))
  }
  if(is.null(sdy_init)) sdy_init = exp(rnorm(1,0.7,0.25))
  
  # deterministic setting of initial values for the temperature and temperature sd estimates
  if((!(is.null(distrmat)) | !(is.null(obsmat))) & is.null(yest_inits)) {
    yest_inits <- c(unlist(sapply(unique(obsmat$sample), function(x) mean(obsmat$temperature[which(obsmat$sample == x)]))),
                    distrmat$mu + distrmat$scale * sqrt(2/pi) * 
                      distrmat$shape/sqrt(1+distrmat$shape^2) )
  }
  if(!(is.null(obsmat)) & is.null(sdyest_inits)) sdyest_inits <- rep(2,length(unique(obsmat$sample)))
  
  
  logprior <- write_logprior(prior_fun = logprior_input, log = TRUE) # create prior
  
  if(is.null(proposal_var_inits)) proposal_var_inits <- c(2^2,2^2,2^2,0.2^2)
  if(is.null(adapt_sd)) adapt_sd <- floor(0.1 * n_iter)
  if(is.null(adapt_sd_decay)) adapt_sd_decay <- max(floor(0.1*adapt_sd),1)
  if(is.null(start_adapt)) start_adapt <- max(floor(0.05*adapt_sd),1)
  
  save_it <- seq(1,n_iter,n_thin) # iterations to save
  
  coefficients = array(dim = c(n_iter,4)) # set up array to store coefficients
  coefficients[1,] = coeff_inits # initialise coefficients
  sdy = rep(NA_real_,n_iter) # set up vector to store sdy
  sdy[1] = sdy_init # intialise sdy
  
  logpostold = NA
  accept = rep(NA,n_iter)
  x = NULL
  if(!is.null(obsmat)) {
    if(is.list(sapply(unique(obsmat$sample), function(x) unique(obsmat$p_lat[which(obsmat$sample == x)])))) warning(
      "something is wrong with the observation matrix. Do all obs from one sample have the same p_lat?")
    
    ylist <- lapply(unique(obsmat$sample), function(x) obsmat$temperature[which(obsmat$sample == x)])
    x <- as.numeric(sapply(unique(obsmat$sample), function(x) unique(obsmat$p_lat[which(obsmat$sample == x)])))
    n_p <- length(ylist)
    n_p_ind <- 1:n_p
    
    ## extra for observations with SD
    if(any(!is.na(obsmat$sd))) {
      
      new_y <- ylist
      
      ylist_sd <- lapply(unique(obsmat$sample), function(x) obsmat$sd[which(obsmat$sample == x)])
      
      ylist_sd_ind <- lapply(ylist_sd,function(x) which(!(is.na(as.numeric(x)))))
    
      sd_obs_ind <- which(!is.na((obsmat$sd)))
    #sd_obs <- obsmat$sd[sd_obs_ind]
    #mean_obs <- obsmat$temperature[sd_obs_ind]
    #n_sd_obs <- length(sd_obs_ind)
    obs_yestimate <-  list(NULL) # set up list of arrays to store estimates
    for(k in 1: length(ylist_sd)) {
      obs_yestimate[[k]] <- array(dim = c(n_iter,length(ylist_sd_ind[[k]]))) 
      obs_yestimate[[k]][1,] <- ylist[[k]][ylist_sd_ind[[k]]] # initialise coefficients with mean (or median)
    }

    }
  } else n_p <- 0
  
  if(!is.null(distrmat)) {
    
    x <- c(x,as.numeric(distrmat$p_lat))
    
    ynorm_mu <- as.numeric(distrmat$mu[which(distrmat$distribution == "normal")])
    ynorm_sd <- as.numeric(distrmat$scale[which(distrmat$distribution == "normal")])
    
    yskew_mu <- as.numeric(distrmat$mu[which(distrmat$distribution == "skew-normal")])
    yskew_sigma <- as.numeric(distrmat$scale[which(distrmat$distribution == "skew-normal")])
    yskew_lambda <- as.numeric(distrmat$shape[which(distrmat$distribution == "skew-normal")])
    
    n_norm <- length(ynorm_mu)
    n_norm_ind <- (n_p+1):(n_p+n_norm)
    
    n_skew <- length(yskew_mu)
    
    n_skew_ind <- (n_p+n_norm+1):(n_p+n_norm+n_skew)
    
    
  } else {
    n_norm <- 0
    n_skew <- 0
  }
  
  nbin <- n_p + n_norm + n_skew
  
  yestimate = array(dim = c(n_iter,nbin)) # set up array to store coefficients
  yestimate[1,] = yest_inits # initialise coefficients
  
  A_sdy = A_sdy # parameter for the prior on the inverse gamma distribution of sdy
  B_sdy = B_sdy # parameter for the prior on the inverse gamma distribution of sdy
  shape_sdy <- A_sdy+nbin/2 # shape parameter for the inverse gamma
  
  sdyest = array(dim = c(n_iter,n_p)) # set up vector to store sdy
  sdyest[1,] = sdyest_inits # intialise sdy
  
  ymean = NULL
  
  # for proposals / adaptation
  proposal_cov <- matrix(0,4,4)
  diag(proposal_cov) <- proposal_var_inits
  proposal_var_inits <- proposal_cov
  
  proposal_factor <- 1 # to adjust acceptance rate
  
  if(is.numeric(adapt_sd)){
  if(adapt_sd<10) stop("adapt_sd needs to be >=10")
  all_weights <- exp((-(adapt_sd-1)):0/adapt_sd_decay)
  adapt_it <- seq(start_adapt,adapt_sd,10) # adapt covariance only at every 10th iteration
  }
  
  if(n_p != 0) {
    #### Investigate these: Need to be broad for single obser
    A_sdyest = 1 # parameter for the prior on the inverse gamma distribution of sdyest
    B_sdyest = 1 # parameter for the prior on the inverse gamma distribution of sdyest
    ####
    yn = sapply(ylist,length)
    #yn[which(is.na(yobs))] = NA
    shape_sdyest =  A_sdyest+yn/2 # shape parameter for the inverse gamma
    ### n-1?!
    yobs_mean= c(sapply(ylist,mean))
    yobs_var = sapply(ylist,var)
    sumobs <- sapply(ylist,sum)
    
  }
  
  if(n_skew != 0) yskew_rho <- -yskew_lambda/sqrt(1+yskew_lambda^2) # not sure why but this needs to be negative
  
  logpost = rep(NA,n_iter)
  # start progress bar
  if (!quiet) cli::cli_progress_bar('Sampling', total = n_iter)
  
  ### The MCMC loop
  for (i in 2:n_iter){
    # update progress bar
    if (!quiet) cli::cli_progress_update(set = i, status = paste0('iteration ', i))
    
    pred = gradient(x,coefficients[i-1,],0)
    
    ### 1. Gibbs steps for data that have multiple points (estimate global mean and sd)
    
    ## 1.0 Gibbs step for data that have SD on point estimates
    ## extra for observations with SD
     
  if(n_p != 0) {  
    if(any(!(is.na(obsmat$sd)))) {
      for(k in 1:length(ylist_sd)){
        if(any(!(is.na(ylist_sd[[k]])))){
          
          n_obs0 <- length(ylist_sd[[k]][ylist_sd_ind[[k]]])
          sdyest1 <-  sdyest[i-1,k]
          sdyobs0 <-  ylist_sd[[k]][ylist_sd_ind[[k]]]
          musdyobs0 <- ylist[[k]][ylist_sd_ind[[k]]]
          muyest1 <- yestimate[i-1,k]
          

      obs_yestimate[[k]][i,] =   rnorm(n_obs0,
                                       sdyest1^2/(sdyobs0^2+sdyest1^2)*musdyobs0 + sdyobs0^2/(sdyobs0^2+sdyest1^2)*muyest1,
                                       sqrt(1/(1/sdyest1^2 + 1/sdyobs0^2))) 
        }
      }
    
      if(i >= 3) new_y_old <- new_y else new_y_old <- unlist(ylist)
      yobs_mean_old <- yobs_mean
      
      new_y <- ylist
      for(k in 1:length(ylist)) new_y[[k]][ylist_sd_ind[[k]]] <- obs_yestimate[[k]][i,]
      
      yobs_mean= c(sapply(new_y,mean))
      sumobs <- sapply(new_y,sum)
      new_y <- unlist(new_y)

      
    } 
  }
      
    ## 1.1.a Gibbs step to estimate yestimate (from data points)
    if(n_p != 0) {
      yestimate[i,1:n_p] = rnorm(n_p,
                                 1/(1/sdy[i-1]^2 + yn/sdyest[i-1,1:n_p]^2)*(pred[1:n_p]/sdy[i-1]^2 + sumobs/sdyest[i-1,1:n_p]^2),
                                 sqrt((1/sdy[i-1]^2 + yn/sdyest[i-1,1:n_p]^2)^(-1)))
      
      ## 1.1.b Gibbs step to estimate sdyest
      for(j in 1:n_p) sdyest[i,j] = sqrt(1/rgamma(1,
                                                  shape_sdyest[j],
                                                  (B_sdyest+0.5*sum((ylist[[j]]-yestimate[i,j])^2))))
    }
    
    ### 1.2. Gibbs steps for data that have sample mean and sd given (estimate global mean only)
    ## 1.2.a Gibbs step to estimate yestimate
    ### distribution from here: https://people.eecs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture5.pdf
    if(n_norm != 0) {
      yestimate[i,n_norm_ind] = rnorm(n_norm,
                                      sdy[i-1]^2/(ynorm_sd^2+sdy[i-1]^2)*ynorm_mu + ynorm_sd^2/(ynorm_sd^2+sdy[i-1]^2)*pred[n_norm_ind],
                                      sqrt(1/(1/sdy[i-1]^2 + 1/ynorm_sd^2))) ## careful check this line!!!!!!!!
    }
    ### 1.3. Gibbs steps for data that have location, scale and shape parameter given (skew-normal). Estimate global mean only
    ## 1.3.a Gibbs step to estimate yestimate
    ### distribution from here: http://koreascience.or.kr/article/JAKO200504840590864.pdf
    if(n_skew != 0) {
      z <- (yskew_mu - yestimate[i-1,n_skew_ind])/yskew_sigma
      y <- truncnorm::rtruncnorm(n_skew,0,Inf,yskew_rho*z,sqrt(1-yskew_rho^2))
      yestimate[i,n_skew_ind] = skew_mu(x=yskew_mu, y=y, sigma=yskew_sigma, rho=yskew_rho,
                                        mu_prior=pred[n_skew_ind], sigma_prior=sdy[i-1])
    }
    
    ### 3. Gibbs step to estimate sdy
    sdy[i] = sqrt(1/rgamma(1,
                           shape_sdy,
                           (B_sdy+0.5*sum((yestimate[i,]-pred)^2))))
    
    ## 4. Metropolis-Hastings step to estimate the regression coefficients
    # 4.0: create proposal innovations
    
    # create matrix of proposal innovations as this is much faster than doing it anew at every it
    if(i == adapt_sd+1 | (i==2 & adapt_sd < 2)) proposal_innovation <-   mvnfast::rmvn(
      n = n_iter-adapt_sd, mu = rep(0,4),sigma = 2.4/sqrt(4)*proposal_cov)+
      rnorm(n = 4*(n_iter-adapt_sd), mean = rep(0,4),sd = 0.001)
    
    # 4.1: create proposals
    if(i <= adapt_sd) proposal_coeff = MH_propose_multi(1,coefficients[i-1,],proposal_cov =  proposal_cov) # new proposed values
    if(i > adapt_sd) proposal_coeff = c(coefficients[i-1,1:3],log(coefficients[i-1,4])) + proposal_innovation[i-adapt_sd,]
    proposal_coeff[4] <- exp(proposal_coeff[4])
    
    if(any(proposal_coeff[4] <= 0)) { # | i == 2
      HR = 0
    } else {# B needs to be >0
      # Hastings ratio of the proposal

      logpostold <- loglik(x = x, y = yestimate[i,], 
                           coeff = coefficients[i-1,],sdy = sdy[i])
      logpostold = logpostold + logprior(coefficients[i-1,])
      

      logpostnew <- loglik(x = x, y = yestimate[i,], 
                          coeff = proposal_coeff,sdy = sdy[i])
      logpostnew = logpostnew + logprior(proposal_coeff)
      
      HR = exp(logpostnew -
                 logpostold -
                 log(coefficients[i-1,4]) +
                 log(proposal_coeff[4]))
    }
    # accept proposal with probability = min(HR,1)
    if (runif(1) < HR){
      accept[i] = 1
      coefficients[i,] = proposal_coeff
      logpost[i] = logpostnew
      # if proposal is rejected, keep the values from the previous iteration
    }else{
      accept[i] = 0
      coefficients[i,] = coefficients[i-1,]
      logpost[i] = logpostold
      
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
    
  } # end of the MCMC loop
    
  # rm(proposal_innovation) #delete matrix of proposal innovations to save memory # not needed
  
  ###  Function output
  output = list(params = data.frame(A = coefficients[save_it,1],
                           dKA = coefficients[save_it,2],
                           M = coefficients[save_it,3],
                           B = coefficients[save_it,4],
                           sdy = sdy[save_it],
                           logpost = logpost[save_it])
                )
  if(!(is.null(obsmat))) {
    output$yestimate = yestimate[save_it,]
    output$sdyest = sdyest[save_it,]
    
  }
    
  
  if(any(!(is.na(obsmat$sd)))) {
    output$obs_yestimate = lapply(1:length(obs_yestimate), function(f) obs_yestimate[[f]][save_it,])
    output$obs_sdy_index = which(sapply(ylist_sd,function(f) any(!(is.na(f)))))
  }
    output$lat = x
    output$proposal_cov = proposal_cov
    
    output$call = mget(names(formals()),sys.frame(sys.nframe()))
  
  return(output)
}