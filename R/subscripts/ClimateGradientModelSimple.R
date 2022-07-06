######################################################
######################################################
### A Metropolis-Hastings algorithm for Latitudinal
### Temperature Gradients - Simple version
###
### first developed: February 2022
### revised: June - July 2022
######################################################
######################################################
### FUNCTIONS

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

loglik_s <- function(x, y,  coeff, sdy) {
  A = coeff[1]
  DKA = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  pred = A + DKA/((1+(exp(Q*(x-M)))))
  return(sum(dnorm(y, mean = pred, sd = sdy, log = TRUE)))
}

logprior_s <- function(coeff) {
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
logposterior_s <- function(x, y, coeff, sdy){
  return (loglik_s(x, y, coeff, sdy) + logprior_s(coeff)) # + 2*(coeff[4])
}


MH_propose_multi <- function(nprop,coeff,proposal_cov) {
  
  
  mvnfast::rmvn(n = nprop, mu = 0.95*c(coeff[1:3],log(coeff[4])),
                sigma = 2.4/sqrt(4)*proposal_cov)+
    rnorm(n = 4*nprop, mean = 0.05*c(coeff[1:3],log(coeff[4])),
          sd = 0.001)
  
}

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



# Main MCMCM function
run_MCMC_simple <- function(x, y, nIter, coeff_inits = NULL, sdy_init = NULL,
                            adapt_sd = floor(0.1 * nIter),
                            adapt_sd_decay = max(floor(0.01*nIter),1),
                            proposal_var_inits = c(2,2,2,0.2)){
  ### Initialisation
  coefficients = array(dim = c(nIter,4)) # set up array to store coefficients
  coefficients[1,] = coeff_inits # initialise coefficients
  sdy = rep(NA_real_,nIter) # set up vector to store sdy
  sdy[1] = sdy_init # intialise sdy
  A_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  n <- length(y)
  shape_sdy <- A_sdy+n/2 # shape parameter for the inverse gamma
  
  # for proposals / adaptation
  proposal_cov <- matrix(0,4,4)
  diag(proposal_cov) <- proposal_var_inits
  proposal_var_inits <- proposal_cov
  
  all_weights <- exp((-(adapt_sd-1)):0/adapt_sd_decay)
  
  proposal_factor <- 1 # to adjust acceptance rate
  
  accept = rep(NA,nIter)

  
  ### The MCMC loop
  for (i in 2:nIter){
    
    ## 1. Gibbs step to estimate sdy
    sdy[i] = sqrt(1/rgamma(
      1,shape_sdy,B_sdy+0.5*sum((y-gradient(x,coefficients[i-1,],0))^2)))
    
    ## 2. Metropolis-Hastings step to estimate the regression coefficients
    if(i <= adapt_sd) proposal_coeff = MH_propose_multi(1,coefficients[i-1,],proposal_cov =  proposal_cov) # new proposed values
    if(i > adapt_sd) proposal_coeff = c(coefficients[i-1,1:3],log(coefficients[i-1,4])) + proposal_innovation[i-adapt_sd,]
    proposal_coeff[4] <- exp(proposal_coeff[4])
    
    #if(any(proposal[4] <= 0)) HR = 0 else # Q needs to be >0
    # Hastings ratio of the proposal
    HR = exp(logposterior_s(x = x, y = y, coeff = proposal_coeff, sdy = sdy[i]) -
               logposterior_s(x = x, y = y, coeff = coefficients[i-1,], sdy = sdy[i]) +
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
          print(paste("stuck in while loop at it",i))}
      }
      
    }
    # create matrix of proposal innovations as this is much faster than doing it anew at every it
    if(i == adapt_sd) proposal_innovation <-   mvnfast::rmvn(n = nIter-adapt_sd, mu = rep(0,4),
                                                             sigma = 2.4/sqrt(4)*proposal_cov)+
      rnorm(n = 4*(nIter-adapt_sd), mean = rep(0,4),
            sd = 0.001)
    
    
  } # end of the MCMC loop
  
  ###  Function output
  output = data.frame(A = coefficients[,1],
                      DKA = coefficients[,2],
                      M = coefficients[,3],
                      Q = coefficients[,4],
                      sdy = sdy)
  return(output)
}
