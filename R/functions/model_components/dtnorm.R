# density of truncated normal distribution, from msm
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
} 

