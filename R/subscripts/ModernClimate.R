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

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

coords <- fastRandomPoints(sstm$BO21_tempmean_ss,10000)

sstr <- raster::extract(sstm,coords)

plot(coords[,2],sstr, pch = 19, col = rgb(0,0,0,0.1))
     
