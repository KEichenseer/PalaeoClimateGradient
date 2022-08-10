# Calculate log likelihood
loglik <- function(x, y,  coeff, sdy) {
  A = coeff[1]
  dKA = coeff[2]
  M = coeff[3]
  B = coeff[4]
  pred = A + dKA/((1+(exp(B*(x-M)))))
  return(sum(dnorm(y, mean = pred, sd = sdy, log = TRUE)))
}
