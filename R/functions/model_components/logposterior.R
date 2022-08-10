# calculate log posterior based on prior and likelihood
logposterior <- function(x, y, coeff, sdy, logprior){
  if(!is.null(y)) {
    return (loglik(x, y, coeff, sdy) + logprior(coeff)) # + 2*(coeff[4])
  } else return(logprior(coeff))
}
