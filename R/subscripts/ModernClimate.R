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
  row = raster::rowFromY(r, c(min,max))
  row = row[1]:row[2]
  cells = raster::cellFromRow(r, row)
  r <- raster::rasterFromCells(r, cells, values = FALSE)
  v <- raster::getValues(r)
  raster::cel
  v.notNA <- which(!is.na(v))
  x <- sample(v.notNA, n)
  pts <- raster::xyFromCell(r, x)
  return(pts)
}
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

coords <- fastRandomPoints_lat(sstm$BO21_tempmean_ss,10,-25,25)

sstr <- raster::extract(sstm,coords)

plot(coords[,2],sstr, pch = 19
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


cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)

mod3 <- climate_simple_parallel(3,100000,abs(coords[,2]),sstr, priorvec)
doParallel::stopImplicitCluster()

plot_chains(mod3)
par(mfrow=c(1,1))
plot_gradient(mod3[[1]])
#coords <- fastRandomPoints(sstm$BO21_tempmean_ss,1000)

#sstr <- raster::extract(sstm,coords)

points(abs(coords[,2]),sstr, pch = 19, col = rgb(0.8,0,0,0.33))

# looking good

# next: use naive priors