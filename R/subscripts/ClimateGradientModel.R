############## Climate gradient model
############## 26/06/2022

gradient <- function(x, coeff, sdy) { # parametrise with difference between cold and hot end instead
  if(is.list(coeff) & !(is.data.frame(coeff) | is.matrix(coeff))) coeff = unlist(coeff)
  if(is.data.frame(coeff) | is.matrix(coeff)) {
    A = coeff[,1]
    DKA = coeff[,2]
    M = coeff[,3]
    Q = coeff[,4]
    
    lat = t(data.frame(lat=x))
    lat = lat[rep(1, each=length(A)),]

    if(sdy == 0) {out = A + DKA/((1+(exp(Q*(lat-M)))))
    } else {
      out = A + DKA/((1+(exp(Q*(lat-M)))))+ rnorm(length(x),0,sdy)
    }
    
  } else {
    A = coeff[1]
    DKA = coeff[2]
    M = coeff[3]
    Q = coeff[4]

  if(sdy == 0) {return(A + DKA/((1+(exp(Q*(x-M))))))
   } else {
   out = A + DKA/((1+(exp(Q*(x-M)))))+ rnorm(length(x),0,sdy)
   }
  }
  return(out)
}

loglik_norm <- function(x, yest, ymean, sdyest, coeff, sdy) {
  # extract regression coefficients
  coeff = unlist(coeff)
  A = coeff[1]
  DKA = coeff[2]
  M = coeff[3]
  Q = coeff[4]

  ll1 <- sum(dnorm(yest, ymean, sdyest,log=TRUE),na.rm=T)
  
  pred = A + DKA/((1+(exp(Q*(x-M)))))
  ll2 <- sum(dnorm(yest, mean = pred, sd = sdy, log = TRUE))
  return(c(ll2+ll1))
}


loglik_skew <- function(x, yest, mu, sigma, lambda, coeff, sdy) {
  # extract regression coefficients
  coeff = unlist(coeff)
  A = coeff[1]
  DKA = coeff[2]
  M = coeff[3]
  Q = coeff[4]

  ll1 <- sum(log(2/sigma)+dnorm((yest-mu)/sigma,log=T)+pnorm(lambda*(yest-mu)/sigma,log=T))
  
  pred = A + DKA/((1+(exp(Q*(x-M)))))
  ll2 <- sum(dnorm(yest, mean = pred, sd = sdy, log = TRUE))
  return(c(ll2+ll1))
}

logprior <- function(coeff) {
  coeff = unlist(coeff)
  return(sum(c(
    dsnorm(coeff[1],location = -2.7, scale = 16, alpha = 16, log = TRUE), # prior on A
    dtnorm(coeff[2], 0, Inf,25,12, log = TRUE), # prior on DKA
    dnorm(coeff[3], 45, 15, log = TRUE), # prior on M
    dlnorm(coeff[4], -2.2, 0.8, log = TRUE)))) # prior on Q
}

# function to generate truncated normal
dtnorm <- function(x,lower,upper,mean,sd, log = FALSE) {
  ret <- numeric(length(x))
  ret[x < lower | x > upper] <- if (log)
    -Inf
  else 0
  ret[upper < lower] <- NaN
  ind <- x >= lower & x <= upper
  if (any(ind)) {
    denom <- pnorm(upper, mean, sd) - pnorm(lower, mean,
                                            sd)
    xtmp <- dnorm(x, mean, sd, log)
    if (log)
      xtmp <- xtmp - log(denom)
    else xtmp <- xtmp/denom
    ret[x >= lower & x <= upper] <- xtmp[ind]
  }
  ret
} # from msm

dsnorm <- function(x,location,scale,alpha, log = TRUE) {
  if(log == TRUE) out = log(2/scale)+dnorm((x - location)/scale,log=T)+pnorm(alpha*(x - location)/scale,log=T)
  if(log == FALSE) out = (2/scale)*dnorm((x - location)/scale,log=F)*pnorm(alpha*(x - location)/scale,log=F)
  return(out)
}


logposterior_norm <- function(x, yest, ymean, sdyest, coeff, sdy){
  return (loglik_norm(x, yest, ymean, sdyest, coeff, sdy))
}

logposterior_skew <- function(x, yest, mu, sigma, lambda, coeff, sdy){
  return (loglik_skew(x, yest, mu, sigma, lambda, coeff, sdy))
}


MH_propose_multi <- function(nprop,coeff,proposal_cov) {
  
  
  mvnfast::rmvn(n = nprop, mu = 0.95*c(coeff[1:3],log(coeff[4])),
                sigma = 2.4/sqrt(4)*proposal_cov)+
    rnorm(n = 4*nprop, mean = 0.05*c(coeff[1:3],log(coeff[4])),
          sd = 0.001)
  
}


# Gibbs sampling of mu with skew normal likelihood and normal prior
### this is for single observations (sample mean), so no mean(x) or mean(y)
skew_mu <- function(x, y, sigma, rho, mu_prior, sigma_prior) {
  n1 = 1
  rnorm(length(x),(n1*sigma_prior^2*(x-sigma*rho*y)+sigma^2*(1-rho^2)*mu_prior)/
          (n1*sigma_prior^2+sigma^2*(1-rho^2)),
        sqrt((sigma_prior^2*sigma^2*(1-rho^2))/(n1*sigma_prior^2+sigma^2*(1-rho^2))) )

}

### generate weighted covariance
weighted_cov <- function(x, weights) { # takes 2.5 times as long as cov()
  dims = ncol(x)
  out = matrix(NA_real_,nrow = dims, ncol = dims)
  
  for(d in 1:dims) {
    out[d,d] <- sum(weights*((x[,d]-sum(weights*x[,d])/sum(weights))^2))/(sum(weights))
    if(d < dims) for(d2 in (d+1):dims) {
      out[d,d2] <- sum(weights*((x[,d]-sum(weights*x[,d])/sum(weights))*(x[,d2]-sum(weights*x[,d2])/sum(weights))))/(sum(weights))
      out[d2,d] <- out[d,d2]
    }
  }
  diag(out) <- diag(out) + 0.000001*diag(out) # add small increment to the two variances (but not to covariance) to ensure chol decomposition does not fail
  
  return(out)
}

#
# Main MCMCM function
run_MCMC <- function(nIter = 1000, obsmat = NULL, distrmat = NULL, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                     proposal_var_inits = c(2,2,2,0.2), adapt_sd = floor(0.1 * nIter),
                     adapt_sd_decay = max(floor(0.01*nIter),1), quiet = FALSE){
  ### Initialisation
  coefficients = array(dim = c(nIter,4)) # set up array to store coefficients
  coefficients[1,] = coeff_inits # initialise coefficients
  sdy = rep(NA_real_,nIter) # set up vector to store sdy
  sdy[1] = sdy_init # intialise sdy
  
  logpostold = NA
  accept = rep(NA,nIter)
  x = NULL
  if(!is.null(obsmat)) {
    ylist <- lapply(unique(obsmat$sample), function(x) obsmat$temperature[which(obsmat$sample == x)])
    x <- as.numeric(sapply(unique(obsmat$sample), function(x) unique(obsmat$latitude[which(obsmat$sample == x)])))
    n_p <- length(ylist)
    n_p_ind <- 1:n_p
  } else n_p <- 0
  
  if(!is.null(distrmat)) {
    
    x <- c(x,as.numeric(distrmat$latitude))
    
    ynorm_mu <- as.numeric(distrmat$location[which(distrmat$distribution == "normal")])
    ynorm_sd <- as.numeric(distrmat$scale[which(distrmat$distribution == "normal")])
    
    yskew_mu <- as.numeric(distrmat$location[which(distrmat$distribution == "skew-normal")])
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
  
  yestimate = array(dim = c(nIter,nbin)) # set up array to store coefficients
  yestimate[1,] = yest_inits # initialise coefficients
  
  A_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  shape_sdy <- A_sdy+nbin/2 # shape parameter for the inverse gamma
  
  sdyest = array(dim = c(nIter,n_p)) # set up vector to store sdy
  sdyest[1,] = sdyest_inits # intialise sdy
  
  ymean = NULL
  
  # for proposals / adaptation
  proposal_cov <- matrix(0,4,4)
  diag(proposal_cov) <- proposal_var_inits
  proposal_var_inits <- proposal_cov
  
  all_weights <- exp((-(adapt_sd-1)):0/adapt_sd_decay)
  
  proposal_factor <- 1 # to adjust acceptance rate
  
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
  
  logpost = rep(NA,nIter)
  # start progress bar
  if (!quiet) cli::cli_progress_bar('Sampling', total = nIter)
  ### The MCMC loop
  for (i in 2:nIter){
    # update progress bar
    if (!quiet) cli::cli_progress_update(set = i, status = paste0('iteration ', i))
    
    pred = gradient(x,coefficients[i-1,],0)
    
    ### 1. Gibbs steps for data that have multiple points (estimate global mean and sd)
    
    ## 1.1.a Gibbs step to estimate yestimate
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
                                      sqrt(1/sdy[i-1]^2 + 1/ynorm_sd^2))
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
    if(i <= adapt_sd) proposal_coeff = MH_propose_multi(1,coefficients[i-1,],proposal_cov =  proposal_cov) # new proposed values
    if(i > adapt_sd) proposal_coeff = c(coefficients[i-1,1:3],log(coefficients[i-1,4])) + proposal_innovation[i-adapt_sd,]
    proposal_coeff[4] <- exp(proposal_coeff[4])
    
    if(any(proposal_coeff[4] <= 0) | i == 2) {
      HR = 0
    } else {# Q needs to be >0
      # Hastings ratio of the proposal
      logpostold = 0
      
      if(n_p != 0) logpostold <- logpostold + logposterior_norm(x = x[n_p_ind], yest = yestimate[i,n_p_ind], ymean = yobs_mean,
                                                                sdyest = sdyest[i,], coeff = coefficients[i-1,],
                                                                sdy = sdy[i])
      if(n_norm != 0) logpostold <- logpostold + logposterior_norm(x = x[n_norm_ind], yest = yestimate[i,n_norm_ind], ymean = ynorm_mu,
                                                                   sdyest = ynorm_sd, coeff = coefficients[i-1,],
                                                                   sdy = sdy[i])
      if(n_skew != 0) logpostold <- logpostold +  logposterior_skew(x = x[n_skew_ind], yest = yestimate[i,n_skew_ind], mu = yskew_mu, yskew_sigma, yskew_lambda,
                                                                    coeff = coefficients[i-1,], sdy[i])
      
      logpostold = logpostold + logprior(coefficients[i-1,])
      
      logpostnew = 0
      
      if(n_p != 0) logpostnew <- logpostnew + logposterior_norm(x = x[n_p_ind], yest = yestimate[i,n_p_ind], ymean = yobs_mean,
                                                                sdyest = sdyest[i,], coeff = proposal_coeff,
                                                                sdy = sdy[i])
      if(n_norm != 0) logpostnew <- logpostnew + logposterior_norm(x = x[n_norm_ind], yest = yestimate[i,n_norm_ind], ymean = ynorm_mu,
                                                                   sdyest = ynorm_sd, coeff = proposal_coeff,
                                                                   sdy = sdy[i])
      if(n_skew != 0) logpostnew <- logpostnew +  logposterior_skew(x = x[n_skew_ind], yest = yestimate[i,n_skew_ind], mu = yskew_mu, yskew_sigma, yskew_lambda,
                                                                    coeff = proposal_coeff, sdy[i])
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
    if(i <= adapt_sd & i>=3) {
      weights = all_weights[(adapt_sd-i+1):adapt_sd]
      proposal_cov <- weighted_cov(cbind(coefficients[1:i,1:3],log(coefficients[1:i,4])),weights = weights)
      if(any(diag(proposal_cov)==0)) proposal_cov <- proposal_var_inits/i
      if(i %in% seq(2*floor(adapt_sd/10), adapt_sd, floor(adapt_sd/10))) {
        #
        if(mean(accept[(i-floor(adapt_sd/10)):i]) < 0.23) proposal_factor <- proposal_factor - 0.1*proposal_factor
        if(mean(accept[(i-floor(adapt_sd/10)):i]) > 0.45) proposal_factor <- proposal_factor + 0.1*proposal_factor
        proposal_cov <- proposal_cov * proposal_factor
        while(any(eigen(proposal_cov)$values <= 0.000001)) {
          diag(proposal_cov) <- 1.25 * diag(proposal_cov)
          # print(paste("stuck in while loop at it",i))}
        }
      
      }
    }
    # create matrix of proposal innovations as this is much faster than doing it anew at every it
    if(i == adapt_sd) proposal_innovation <-   mvnfast::rmvn(n = nIter-adapt_sd, mu = rep(0,4),
                                                             sigma = 2.4/sqrt(4)*proposal_cov)+
      rnorm(n = 4*(nIter-adapt_sd), mean = rep(0,4),
            sd = 0.001)
    
  } # end of the MCMC loop
    
  rm(proposal_innovation) #delete matrix of proposal innovations to save memory
  
  ###  Function output
  output = list(params = data.frame(A = coefficients[,1],
                           DKA = coefficients[,2],
                           M = coefficients[,3],
                           Q = coefficients[,4],
                           sdy = sdy,
                           logpost = logpost),
                yestimate = yestimate,
                sdyest = sdyest,
                lat = x)
  return(output)
}