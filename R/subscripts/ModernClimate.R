### load climate data
# Mean annual sea surface temperatures (Bio-Oracle)

sstm <- read "data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc"

# sample points randomly from the raster (within specified latitudinal limits) 
source("R/subscripts/AuxiliaryFunctions.R")
coords <- fastRandomPoints_lat(sstm$BO21_tempmean_ss,200,-90,90)

# extract temperatures from the sampled points
sstr <- raster::extract(sstm,coords)

plot(abs(coords[,2]),sstr, pch = 19
     , col = rgb(0,0,0,0.1))



modm <- climate_simple_parallel(nChains = nClust, nIter = 10000, nThin = 10,
                                x = seq(0,90,10), y = seq(30,0,length.out = 10)+rnorm(10,0,2), 
                                coeff_inits = NULL, sdy_init = NULL, 
                                prior_fun = prior_fun,
                                proposal_var_inits = c(2,2,2,0.2), adapt_sd = 1000,
                                adapt_sd_decay = 100, start_adapt = 101, quiet = FALSE)

### Test modern data
library(foreach)
source("R/subscripts/ClimateParallelSimple.R")

cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)

modm <- climate_simple_parallel(3,40000,NULL,NULL,priorvec)
modm <- climate_simple_parallel(3,20000,abs(coords[,2]),c(sstr), priorvec,
                                adapt_sd_decay = 100, adapt_sd = 4000)
doParallel::stopImplicitCluster()



plot_chains(mod0b)
par(mfrow=c(1,1))
plot_gradient(modm[[1]], add = T, ylim = c(-5,55), line_col = rgb(0,0,0,0.5), confint_col = rgb(0,0,0,0.1))
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