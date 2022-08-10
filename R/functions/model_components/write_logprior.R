# creating the logprior() function based on prior input from the user
write_logprior <- function(prior_fun,log=TRUE) {
  out <- function(coeff) {
    coeff = unlist(coeff)
    
    return(sum(c(
      prior_fun$f1(x=coeff[1],log=log),
      prior_fun$f2(x=coeff[1]+coeff[2],lower=coeff[1],log=log),
      prior_fun$f3(x=coeff[3],log=log),
      prior_fun$f4(x=coeff[4],log=log))))
  }
  return(out)
}
