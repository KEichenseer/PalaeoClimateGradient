### auxiliary functions

error_polygon <- function(x,en,ep,color=rgb(0,0,0,0.2)) {
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
}


gradient <- function(x, coeff, sdy) { # parametrise with difference between cold and hot end instead
  if(is.list(coeff) & !(is.data.frame(coeff) | is.matrix(coeff))) coeff = unlist(coeff)
  if(is.data.frame(coeff) | is.matrix(coeff)) {
    A = coeff[,1]
    DKA = coeff[,2]
    M = coeff[,3]
    Q = coeff[,4]
    
    lat = t(data.frame(lat=x))
    lat = lat[rep(1, each=length(A)),]
    
    if(sdy == 0) {out = A + DKA/((1+(exp(Q*(lat-M)))))
    } else {
      out = A + DKA/((1+(exp(Q*(lat-M)))))+ rnorm(length(x),0,sdy)
    }
    
  } else {
    A = coeff[1]
    DKA = coeff[2]
    M = coeff[3]
    Q = coeff[4]
    
    if(sdy == 0) {return(A + DKA/((1+(exp(Q*(x-M))))))
    } else {
      out = A + DKA/((1+(exp(Q*(x-M)))))+ rnorm(length(x),0,sdy)
    }
  }
  return(out)
}


plot_gradient <- function(model_out, burnin = NULL, lat = seq(0,90,0.2), confint_n = NULL, add = F,
                          ylim = NULL, line_col = "black", confint_col = rgb(0,0,0,0.2), lwd = 2) {
  # select only params
  if("params" %in% names(model_out)) model_out <- model_out$params
  
  n_iter <- nrow(model_out)
  if(is.null(burnin)) burnin = round(nrow(model_out)*12/100)+1
  if(is.null(confint_n)) confint_n = min(5000,length(burnin:n_iter))
    
  sample_it <- sample(burnin:n_iter,confint_n)
  grad_q <- apply(gradient(lat,model_out[sample_it,1:4],0),2,function(a) 
    quantile(a,probs = c(0.025,0.975)))

  med_grad <- gradient(lat,apply(model_out[burnin:nrow(model_out),1:4],2,median),0)
  
  if(is.null(ylim)) ylim <- 
    c(max(min(med_grad*0.75),min(grad_q["2.5%",])),min(max(med_grad*1.25),max(grad_q["97.5%",])))
  
  if(add == F) {
  plot(lat, med_grad, ylim = ylim,
       type = "l", lwd = lwd, xlab = "|latitude|", ylab = "temperature", col = line_col)
  error_polygon(lat,grad_q["2.5%",],grad_q["97.5%",], col = confint_col)
  }
  if(add == T) {
    error_polygon(lat,grad_q["2.5%",],grad_q["97.5%",], col = confint_col)
    points(lat, med_grad,
           type = "l", lwd = lwd, col = line_col)
  }
}

plot_sample_gradient <- function(model_out, burnin = NULL, lat = seq(0,90,0.2), confint_n = NULL, add = F,
                          ylim = NULL, line_col = rgb(0,0.33,1,0.33), n_samples = 8, return_data = F, plot = T) {
  # select only params
  if("params" %in% names(model_out)) model_out <- model_out$params
  
  n_iter <- nrow(model_out)
  if(is.null(burnin)) burnin = round(nrow(model_out)*12/100)+1
  #if(is.null(confint_n)) confint_n = min(5000,length(burnin:n_iter))
  
  sample_it <- sample(burnin:n_iter,n_samples)
  
  
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
      tdens <- dnorm(tseq,distrmat$mu[i],distrmat$scale[i])
    }
    if(distrmat$distribution[i] == "skew-normal") {
      tdens <- dsnorm(tseq,distrmat$mu[i],distrmat$scale[i],distrmat$shape[i], log = FALSE)
    }
    effrange <- which(tdens > 0.01*max(tdens))
  twidth <- distrwidth*(max(lat)-min(lat))
  y = tseq[effrange]
  x1 = distrmat$p_lat[i] - twidth*tdens[effrange]/max(tdens[effrange])
  x2 = distrmat$p_lat[i] + twidth*tdens[effrange]/max(tdens[effrange])
  
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
                           ylim = NULL, col_obs = rgb(.8,0.5,0,0.75), col_dist = rgb(0,0.5,.8,0.75),
                           cex = 1, pch = 17) {
  n_iter <- nrow(mod$params)
  if(is.null(burnin)) burnin = round(nrow(mod$params)*12/100)+1
  
  if(!(is.null(mod$sdyest))) {
  nobsloc <- dim(mod$sdyest)[2]
  samples <- 1:dim(mod$sdyest)[2]
  if(length(col_obs)!=nobsloc)  col_obs <- rep(col_obs[1],nobsloc)
  
  invisible(sapply(samples,function(x) points(mod$lat[x], mean(mod$yestimate[burnin:n_iter,x]), pch = pch, col = col_obs[x], cex = cex)))
  
  invisible(sapply(samples,function(x) points(rep(mod$lat[x],2), quantile(mod$yestimate[burnin:n_iter,x], probs = c(0.05,0.95)), 
                                    type = "l", col = col_obs[x])))
 
  } else nobsloc <- 0
  
  if(dim(mod$yestimate)[2] > nobsloc) {
    ndistrloc <- dim(mod$yestimate)[2]
  
    
  
  distr_ind <- (nobsloc+1):ndistrloc
  invisible(sapply(distr_ind,function(x) points(mod$lat[x], median(mod$yestimate[burnin:n_iter,x]), pch = pch, col = col_dist, cex = cex)))
  invisible(sapply(distr_ind,function(x) points(rep(mod$lat[x],2), quantile(mod$yestimate[burnin:n_iter,x], probs = c(0.05,0.95)), 
                                                type = "l", col = col_dist)))
  }
}



plot_chains <- function(mod, params = 1:4, n_thin = NULL, logB = TRUE) {
  if("params" %in% names(mod[[1]])) {mod <- lapply(mod, function(x) x$params)
  modnames <- names(mod[[1]])
  if(!is.numeric(params)) plot_it <- which(modnames %in% params) else plot_it <- params
  }
  op <- par()[c("mfrow","mar","mgp")]
  nplot <- length(params)
  cols <- c(rgb(0,0.5,0.75,0.7),
            rgb(0.75,0,0.5,0.7),
            rgb(0.77,0.5,0,0.7),
            rgb(0,0.75,0,0.7))
  par(mfrow = c(nplot,1), mar  = c(3.5,3.5,0.5,0.5), mgp = c(2.25,0.75,0), las = 1)
  if(!("data.frame" %in% class(mod))) {
    nchains <- length(mod)
    n_iter <- nrow(mod[[1]])
    if(is.null(n_thin)) n_thin <- max(c(1,round(n_iter/2000)))
        iteration <- seq(1,n_iter,n_thin)
    for(j in plot_it) {
      if(j != 4 | logB == FALSE){
      plot(iteration,mod[[1]][iteration,j],type = "l",
                           col = cols[1], ylab = modnames[j])
      if(nchains >= 2)  for(i in 2:nchains) {
        points(iteration,mod[[i]][iteration,j],type = "l",
             col = cols[i])
      }
      }
      if(j == 4 & logB == TRUE){
        plot(iteration,log10(mod[[1]][iteration,j]),type = "l",
             col = cols[1], ylab = "log10(B)")
        if(nchains >= 2)  for(i in 2:nchains) {
          points(iteration,log10(mod[[i]][iteration,j]),type = "l",
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

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}
dsnorm <- function(x,location,scale,alpha, log = TRUE) {
  if(log == TRUE) out = log(2/scale)+dnorm((x - location)/scale,log=T)+pnorm(alpha*(x - location)/scale,log=T)
  if(log == FALSE) out = (2/scale)*dnorm((x - location)/scale,log=F)*pnorm(alpha*(x - location)/scale,log=F)
  return(out)
}


# function to generate truncated normal
dtnorm <- function(x,lower,upper,mean,sd, log = FALSE) {
  ret <- numeric(length(x))
  ret[x < lower | x > upper] <- if (log)
    -Inf
  else 0
  ret[upper < lower] <- NaN
  ind <- x >= lower & x <= upper
  if (any(ind)) {
    denom <- pnorm(upper, mean, sd) - pnorm(lower, mean,
                                            sd)
    xtmp <- dnorm(x, mean, sd, log)
    if (log)
      xtmp <- xtmp - log(denom)
    else xtmp <- xtmp/denom
    ret[x >= lower & x <= upper] <- xtmp[ind]
  }
  ret
} # from msm


plot_prior <- function(prior,xval,log=FALSE) {
  par(mfrow=c(2,2), mar = c(4.1,4.1,1,1), mgp = c(2.5,0.75,0))
  dens <- prior$f1(xval[[1]],log=log)
  plot_dens(xval[[1]],dens,xlab = "A")
  
  dens <- prior$f2(xval[[2]],xval[[1]][which.max(dens)],log=log)
  plot_dens(xval[[2]],dens,xlab = "K")
  
  dens <- prior$f3(xval[[3]],log=log)
  plot_dens(xval[[3]],dens,xlab = "M")
  
  dens <- prior$f4(xval[[4]],log=log)
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

combine_posterior <- function(mod, burnin = NULL) {
  n_thin <- mod[[1]]$call$n_thin
  n_iter <- floor(mod[[1]]$call$n_iter/n_thin)
  if(is.null(burnin)) burnin <- 0 else burnin <- floor(burnin/n_thin)
  out <- do.call(rbind,lapply(1:length(mod),function(f) mod[[f]]$params[(burnin+1):n_iter,]))
}
