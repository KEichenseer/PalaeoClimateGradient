### load climate data
# Mean annual sea surface temperatures (Bio-Oracle)

sstm <- sdmpredictors::load_layers(c("BO21_tempmean_ss")) 

fastRandomPoints <- function(r, n) {
  if(raster::nlayers(r) > 1) r <- r[[1]]
  v <- raster::getValues(r)
  v.notNA <- which(!is.na(v))
  x <- sample(v.notNA, n)
  pts <- raster::xyFromCell(r, x)
  return(pts)
}
fastRandomPoints_lat <- function(r, n, min, max) {
  #if(raster::nlayers(r) > 1) r <- r[[1]]
  #row = raster::rowFromY(r, c(min,max))
  #row = row[1]:row[2]
  #cells = raster::cellFromRow(r, row)
  #r <- raster::rasterFromCells(r, cells, values = FALSE)
  extent <- c(-180,180,min,max)
  r <- raster::crop(r,extent)
  v <- raster::getValues(r)
  v.notNA <- which(!is.na(v))
  x <- sample(v.notNA, n)
  pts <- raster::xyFromCell(r, x)
  return(pts)
}
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

coords <- fastRandomPoints_lat(sstm$BO21_tempmean_ss,1000,-90,90)

sstr <- raster::extract(sstm,coords)

plot(abs(coords[,2]),sstr, pch = 19
     , col = rgb(0,0,0,0.1))
     
### Test modern data
library(foreach)
source("R/subscripts/AuxiliaryFunctions.R")
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



latx <- seq(-40,60,0.1)
dens <- prior_dens(latx,priorvec,1)
plot_dens(latx,dens,xlab="A")

latx <- seq(-40,60,0.1)
dens <- prior_dens(latx,priorvec,2)
plot_dens(latx,dens,xlab="K")

latx <- seq(-45,135,0.1)
dens <- prior_dens(latx,priorvec,3)
plot_dens(latx,dens,xlab="Q")

latx <- seq(0,2,0.002)
dens <- prior_dens(latx,priorvec,4)
plot_dens(latx,dens,xlab="M")


cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)

mod3 <- climate_simple_parallel(3,40000,abs(coords[,2]),sstr, priorvec)
doParallel::stopImplicitCluster()



plot_chains(mod3)
par(mfrow=c(1,1))
plot_gradient(mod3[[1]])
points(abs(coords[,2]),sstr, pch = 19, col = rgb(0.8,0,0,0.33))

plot_dens(seq(0,90,0.1),prior_dens(seq(0,90,0.1),priorvec,1))

prior_dens(seq(0,90,0.1),priorvec,3)
exp(eval(parse(text = priorvec[3])))

make_prior(priorvec)
run_MCMC_simple(NULL,NULL,1000,coeff_inits = c(0,0,45,1),sdy_init = 1)
#coords <- fastRandomPoints(sstm$BO21_tempmean_ss,1000)

#sstr <- raster::extract(sstm,coords)


# looking good

# next: use naive priors