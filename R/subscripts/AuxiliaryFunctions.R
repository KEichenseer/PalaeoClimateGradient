### auxiliary functions

error_polygon <- function(x,en,ep,color=rgb(0,0,0,0.2)) {
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
}

plot_gradient <- function(model_out, burnin = NULL, lat = seq(0,90,0.2), confint_n = NULL, add = F,
                          ylim = NULL, line_col = "black", confint_col = rgb(0,0,0,0.2)) {
  # select only params
  if("params" %in% names(model_out)) model_out <- model_out$params
  
  nIter <- nrow(model_out)
  if(is.null(burnin)) burnin = round(nrow(model_out)*12/100)+1
  if(is.null(confint_n)) confint_n = min(5000,length(burnin:nIter))
    
  sample_it <- sample(burnin:nIter,confint_n)
  grad_q <- apply(gradient(lat,model_out[sample_it,1:4],0),2,function(a) 
    quantile(a,probs = c(0.025,0.975)))

  med_grad <- gradient(lat,apply(model_out[burnin:nrow(model_out),1:4],2,median),0)
  
  if(is.null(ylim)) ylim <- 
    c(max(min(med_grad*0.75),min(grad_q["2.5%",])),min(max(med_grad*1.25),max(grad_q["97.5%",])))
  
  if(add == F) {
  plot(lat, med_grad, ylim = ylim,
       type = "l", lwd = 2, xlab = "|latitude|", ylab = "temperature", col = line_col)
  error_polygon(lat,grad_q["2.5%",],grad_q["97.5%",], col = confint_col)
  }
  if(add == T) {
    error_polygon(lat,grad_q["2.5%",],grad_q["97.5%",], col = confint_col)
    points(lat, med_grad,
           type = "l", lwd = 2, col = line_col)
  }
}

plot_sample_gradient <- function(model_out, burnin = NULL, lat = seq(0,90,0.2), confint_n = NULL, add = F,
                          ylim = NULL, line_col = rgb(0,0.33,1,0.33), n_samples = 8, return_data = F, plot = T) {
  # select only params
  if("params" %in% names(model_out)) model_out <- model_out$params
  
  nIter <- nrow(model_out)
  if(is.null(burnin)) burnin = round(nrow(model_out)*12/100)+1
  #if(is.null(confint_n)) confint_n = min(5000,length(burnin:nIter))
  
  sample_it <- sample(burnin:nIter,n_samples)
  
  
  grads <- gradient(lat,model_out[sample_it,1:4],0)
  
  if(plot==T) for(i in 1:n_samples) points(lat, grads[i,],
           type = "l", lwd = 2, col = line_col)
  if(return_data==T) return(grads)
}

plot_distr <- function(distrmat, lat = seq(0,90,0.2), trange = c(-10,100), distrwidth = 0.01,
                       col = rgb(0.15,0,0.7,0.25)) {
  
  tseq <- seq(trange[1],trange[2],0.1)
  
  for(i in 1:nrow(distrmat)) {
    if(distrmat$distribution[i] == "normal") {
      tdens <- dnorm(tseq,distrmat$location[i],distrmat$scale[i])
    }
    if(distrmat$distribution[i] == "skew-normal") {
      tdens <- dsnorm(tseq,distrmat$location[i],distrmat$scale[i],distrmat$shape[i], log = FALSE)
    }
    effrange <- which(tdens > 0.01*max(tdens))
  twidth <- distrwidth*(max(lat)-min(lat))
  y = tseq[effrange]
  x1 = distrmat$latitude[i] - twidth*tdens[effrange]/max(tdens[effrange])
  x2 = distrmat$latitude[i] + twidth*tdens[effrange]/max(tdens[effrange])
  
  polygon(c(x1,rev(x2)),c(y,rev(y)), col = col, border = NA)
  }
}

plot_data <- function(obsmat = NULL, distrmat = NULL, lat = seq(0,90,0.2), add = F, ylim = NULL,
                      col = rgb(0.75,0.45,0,0.4)) {
  if(add == T) {
    points(obsmat$latitude,obsmat$temperature, bg = col, pch = 21, col = NA)
  }
  if(add == F) {
    plot(obsmat$latitude,obsmat$temperature, bg = col, pch = 21, col = NA,
         xlim = range(lat), ylim = ylim,
         xlab = "|latitude|", ylab = "temperature")
  }
  
}

plot_posterior <- function(mod,burnin = NULL, lat = seq(0,90,0.2), confint_n = NULL, add = F,
                           ylim = NULL, col_obs = rgb(.85,0.3,0,0.75), col_dist = rgb(0.6,0,.9,0.75),
                           cex = 1, pch = 17) {
  nIter <- nrow(mod$params)
  if(is.null(burnin)) burnin = round(nrow(mod$params)*12/100)+1
  
  if(!(is.null(mod$sdyest))) {
  nobsloc <- dim(mod$sdyest)[2]
  samples <- 1:dim(mod$sdyest)[2]
  if(length(col_obs)!=nobsloc)  col_obs <- rep(col_obs[1],nobsloc)
  
  invisible(sapply(samples,function(x) points(mod$lat[x], mean(mod$yestimate[burnin:nIter,x]), pch = pch, col = col_obs[x], cex = cex)))
  
  invisible(sapply(samples,function(x) points(rep(mod$lat[x],2), quantile(mod$yestimate[burnin:nIter,x], probs = c(0.05,0.95)), 
                                    type = "l", col = col_obs[x])))
 
  } else nobsloc <- 0
  
  if(dim(mod$yestimate)[2] > nobsloc) {
    ndistrloc <- dim(mod$yestimate)[2]
  
    
  
  distr_ind <- (nobsloc+1):ndistrloc
  invisible(sapply(distr_ind,function(x) points(mod$lat[x], median(mod$yestimate[burnin:nIter,x]), pch = pch, col = col_dist, cex = cex)))
  invisible(sapply(distr_ind,function(x) points(rep(mod$lat[x],2), quantile(mod$yestimate[burnin:nIter,x], probs = c(0.05,0.95)), 
                                                type = "l", col = col_dist)))
  }
}

coral_distrmat <- function(data,stage) { # for PARED
  
  data_sub <- subset(data,early_stage == stage & !(is.na(pal_lat_scotese)))#table(iso$stage_2020)
  data_sub <- data_sub[with(data_sub, order(abs(pal_lat_scotese), longit)),]
  # prepare for use in the model
  if(nrow(data_sub) >= 1) {data_distrmat <- data.frame(latitude = abs(data_sub$pal_lat_scotese),
                                                         location = 22.8,
                                                         scale = 10,
                                                         shape = 4,
                                                         distribution = "skew-normal")
  } else data_distrmat <- data.frame(NULL)
  
}
  
iso_obsmat <- function(data, stage) { # for StabisoDB
  # select data from one stage to test, exclude NA data
  data_sub <- subset(data,stage_2020 == stage & !(is.na(paleolat)) & !(is.na(temperature)))
  data_sub <- data_sub[with(data_sub, order(abs(paleolat), longitude)),]
  # prepare for use in the model
  data_mod <- data.frame(sample = (paste(abs(data_sub$paleolat),data_sub$longitude)),
                        latitude = abs(data_sub$paleolat), temperature = data_sub$temperature)
  return(data_mod)
}


plot_chains <- function(mod, params = 1:4, nthin = NULL, logQ = TRUE) {
  if("params" %in% names(mod[[1]])) mod <- lapply(mod, function(x) x$params)
  op <- par()[c("mfrow","mar","mgp")]
  nplot <- length(params)
  cols <- c(rgb(0,0.5,0.75,0.7),
            rgb(0.75,0,0.5,0.7),
            rgb(0.77,0.67,0,0.7),
            rgb(0,0.75,0,0.7))
  par(mfrow = c(nplot,1), mar  = c(3.5,3.5,0.5,0.5), mgp = c(2.25,0.75,0), las = 1)
  if(!("data.frame" %in% class(mod))) {
    nchains <- length(mod)
    nIter <- nrow(mod[[1]])
    if(is.null(nthin)) nthin <- round(nIter/2000)
        iteration <- seq(1,nIter,nthin)
    for(j in 1:nplot) {
      if(j != 4 | logQ == FALSE){
      plot(iteration,mod[[1]][iteration,params[j]],type = "l",
                           col = cols[1], ylab = names(mod[[1]][j]))
      if(nchains >= 2)  for(i in 2:nchains) {
        points(iteration,mod[[i]][iteration,params[j]],type = "l",
             col = cols[i])
      }
      }
      if(j == 4 & logQ == TRUE){
        plot(iteration,log10(mod[[1]][iteration,params[j]]),type = "l",
             col = cols[1], ylab = "log10(Q)")
        if(nchains >= 2)  for(i in 2:nchains) {
          points(iteration,log10(mod[[i]][iteration,params[j]]),type = "l",
                 col = cols[i])
        }
      }

    }
  }
  par(op)
}

prior_dens <- function(x,priorvec,priorindex) {
  x = x
  exp(eval(parse(text = priorvec[priorindex])))
}

plot_dens <- function(x,dens,xlim = NULL,ylim = NULL,col = rgb(0,0.5,0.7,0.25),
                      xlab = "", ylab = "normalised density", add = FALSE, xaxs = "i") {
  dens <- dens/max(dens)
  if(is.null(ylim)) ylim = c(0,max(dens)*1.03)
  if(is.null(xlim)) xlim = range(x)
  
  if(add == FALSE) plot(0,0,type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, xaxs = xaxs, yaxs = "i")
  polygon(c(x[1],x,x[length(x)]),c(0,dens,0), col = col, border = NA)
  
}

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
dsnorm <- function(x,location,scale,alpha, log = TRUE) {
  if(log == TRUE) out = log(2/scale)+dnorm((x - location)/scale,log=T)+pnorm(alpha*(x - location)/scale,log=T)
  if(log == FALSE) out = (2/scale)*dnorm((x - location)/scale,log=F)*pnorm(alpha*(x - location)/scale,log=F)
  return(out)
}

plot_prior <- function(prior,xval,log=FALSE) {
  par(mfrow=c(2,2), mar = c(4.1,4.1,1,1), mgp = c(2.5,0.75,0))
  dens <- prior_fun$f1(xval[[1]],log=log)
  plot_dens(xval[[1]],dens,xlab = "A")
  
  dens <- prior_fun$f2(xval[[2]],xval[[1]][which.max(dens)],log=log)
  plot_dens(xval[[2]],dens,xlab = "K")
  
  dens <- prior_fun$f3(xval[[3]],log=log)
  plot_dens(xval[[3]],dens,xlab = "M")
  
  dens <- prior_fun$f4(xval[[4]],log=log)
  plot_dens(xval[[4]],dens,xlab = "Q")
  
  par(mfrow=c(1,1))
  
}

plot_accept <- function(variable,binsize,return_val = FALSE) {
  n <- floor(length(variable)/binsize)
  out <- rep(NA,n)
  for(i in 1:n) out[i] <- length(which(diff(variable[(binsize*(i-1)+1):(binsize*i)]) != 0))/binsize
  plot(seq(1,n*binsize,by=binsize),out,type = "o", pch = 19, cex = 0.5)
  if(return_val) return(out)
}
