# calculate weighted covariance for proposal covariance adaptation
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