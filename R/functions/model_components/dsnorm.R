# density of skew-normal distribution
dsnorm <- function(x,location,scale,alpha, log = TRUE) {
  if(log == TRUE) out = log(2/scale)+dnorm((x - location)/scale,log=T)+pnorm(alpha*(x - location)/scale,log=T)
  if(log == FALSE) out = (2/scale)*dnorm((x - location)/scale,log=F)*pnorm(alpha*(x - location)/scale,log=F)
  return(out)
}