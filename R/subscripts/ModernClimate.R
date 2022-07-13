### load climate data
# Mean annual sea surface temperatures (Bio-Oracle)

sstm <- sdmpredictors::load_layers(c("BO21_tempmean_ss")) 

# sample points randomly from the raster (within specified latitudinal limits) 
source("R/subscripts/AuxiliaryFunctions.R")
coords <- fastRandomPoints_lat(sstm$BO21_tempmean_ss,10000,-90,90)

# extract temperatures from the sampled points
sstr <- raster::extract(sstm,coords)

plot(abs(coords[,2]),sstr, pch = 19
     , col = rgb(0,0,0,0.1))
     
### Test modern data
library(foreach)
source("R/subscripts/ClimateParallelSimple.R")
source("R/subscripts/ClimateGradientModelSimple.R")


priorvec <- 
c("dsnorm(x,location = -2.7, scale = 16, alpha = 16, log = TRUE)", # prior on A
"dtnorm(x, 0, Inf,25,12, log = TRUE)", # prior on DKA
"dnorm(x, 45, 15, log = TRUE)", # prior on M
"dlnorm(x, -2.2, 0.8, log = TRUE)") # prior on Q


priorvec <- 
  c("dnorm(x,10,15,log=TRUE)", # prior on A
    "dnorm(x,10,15,log = TRUE)", # prior on DKA
    "dnorm(x, 45, 25, log = TRUE)", # prior on M
    "dlnorm(x, -2, 1, log = TRUE)") # prior on Q

priorvec <- 
  c("dsnorm(x,location = 4.3, scale = 12, alpha = 3, log = TRUE)", # prior on A
    "dtnorm(x, 0, Inf,16,10, log = TRUE)", # prior on DKA
    "dnorm(x, 45, 15, log = TRUE)", # prior on M
    "dlnorm(x, -2.2, 0.75, log = TRUE)")#dsnorm(x,location = 0.025, scale = 0.3, alpha = 20, log = TRUE)") # prior on Q

par(mfrow=c(2,2), mar = c(4,4,.5,.5), mgp = c(2.5,0.75,0))
latx <- seq(-10,48,0.001)
dens <- prior_dens(latx,priorvec,1)
plot_dens(latx,dens,xlab="A")

latx <- seq(-3,55,0.1)
dens <- prior_dens(latx,priorvec,2)
plot_dens(latx,dens,xlab="K")

latx <- seq(-1,91,0.1)
dens <- prior_dens(latx,priorvec,3)
plot_dens(latx,dens,xlab="Q")

latx <- seq(-0.01,0.56,0.002)
dens <- prior_dens(latx,priorvec,4)
#dens2 <- exp(dsnorm(latx,location = 0.025, scale = 0.3, alpha = 20, log = TRUE))
plot_dens(latx,dens,xlab="M", add = F)
#plot_dens(latx,dens,xlab="M", add = T, col = rgb(0.75,0,0,0.2))

plot(latx,gradient(latx,c(20,30,45,0.2),0), type = "l")

cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)

mod0b <- climate_simple_parallel(3,40000,NULL,NULL,priorvec)
mod10c <- climate_simple_parallel(3,40000,abs(coords[,2][1:10]),c(sstr[1:10]), priorvec,
                                adapt_sd_decay = 100, adapt_sd = 4000)
doParallel::stopImplicitCluster()



plot_chains(mod4)
par(mfrow=c(1,1))
plot_gradient(mod0b[[2]], add = F, ylim = c(-5,55), line_col = rgb(0,0,0,0.5), confint_col = rgb(0,0,0,0.1))
plot_gradient(mod2[[2]], add = T, line_col = rgb(0,0.3,0.8,1), confint_col = rgb(0,0.3,0.8,0.2))
plot_gradient(mod5[[2]], add = T, line_col = rgb(0.8,0.3,0,1), confint_col = rgb(0.8,0.3,0,0.2))
plot_gradient(mod10[[2]], add = T, line_col = rgb(0.3,0.8,0,1), confint_col = rgb(0.3,0.8,0,0.2))
plot_gradient(mod10b[[2]], add = T, line_col = rgb(0.3,0,0.8,1), confint_col = rgb(0.3,0,0.8,0.2))
plot_gradient(mod10c[[2]], add = T, line_col = rgb(0,0.2,0.8,1), confint_col = rgb(0,0.2,0.8,0.2))

points(abs(coords[,2])[1:10],sstr[1:10], pch = 19, col = rgb(0.3,0.8,0,0.5))

points(abs(coords[,2]),sstr, pch = 19, col = rgb(0.8,0,0,0.03))

plot_dens(seq(0,90,0.1),prior_dens(seq(0,90,0.1),priorvec,1))

prior_dens(seq(0,90,0.1),priorvec,3)
exp(eval(parse(text = priorvec[3])))

make_prior(priorvec)
run_MCMC_simple(NULL,NULL,1000,coeff_inits = c(0,0,45,1),sdy_init = 1)
#coords <- fastRandomPoints(sstm$BO21_tempmean_ss,1000)

#sstr <- raster::extract(sstm,coords)


# looking good

# next: use naive priors