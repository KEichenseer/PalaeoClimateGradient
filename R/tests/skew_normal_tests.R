x <- seq(-5,60,.1)

scale <- 12
location <- -2.85
alpha <- 10

lskd <- function(x,location,scale,alpha) log(2/scale)+dnorm((x - location)/scale,log=T)+pnorm(alpha*(x - location)/scale,log=T)
lskt <- function(x,location,scale,alpha) log(2/scale)+dnorm((x - location)/scale,log=T)+pnorm(alpha*(x - location)/scale,log=T)

skd <- function(x,location,scale,alpha) (2/scale)*dnorm((x - location)/scale,log=F)*pnorm(alpha*(x - location)/scale,log=F)

plot(x,  log(2/scale)+dnorm((x - location)/scale,log=T)+pnorm(alpha*(x - location)/scale,log=T),
type = "l")
points(x,  lskd(x,-3,15,12),
     type = "l", col = "red")
cbind(x,lskd(x,-3.05,15,12))

plot(x,  skd(x,24,13,4), col = "blue", type = "l", yaxs = "i")
abline(v=0)     

points(x,dtnorm(x, 0, Inf,25,15, log = FALSE), col = "red", type = "l",)
abline(v=-4)
points(x,0.76*skd(x,22.8,10,4), col = "green", type = "l")

x[which.max(skd(x,22.8,10,4))]
x[which.max(skd(x,37.09,24,-2.5))]
