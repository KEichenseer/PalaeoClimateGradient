### auxiliary functions

error_polygon <- function(x,en,ep,color=rgb(0,0,0,0.2)) {
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
}

plot_gradient <- function(model_out, burnin = NULL, lat = seq(0,90,0.2), confint_n = NULL, add = F,
                          ylim = NULL) {
  
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
       type = "l", lwd = 2, xlab = "|latitude|", ylab = "temperature")
  error_polygon(lat,grad_q["2.5%",],grad_q["97.5%",])
  }
  if(add == T) {
    error_polygon(lat,grad_q["2.5%",],grad_q["97.5%",])
    points(lat, med_grad, ylim = ylims,
           type = "l", lwd = 2)
  }
}

plot_distr <- function(distrmat, lat = seq(0,90,0.2), trange = c(-10,100), distrwidth = 0.01,
                       col = rgb(0,0.5,0.75,0.33)) {
  
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
                      col = rgb(0.75,0.4,0,0.33)) {
  if(add == T) {
    points(obsmat$latitude,obsmat$temperature, bg = col, pch = 21, col = NA)
  }
  if(add == F) {
    plot(obsmat$latitude,obsmat$temperature, bg = col, pch = 21, col = NA,
         xlim = range(lat), ylim = ylim,
         xlab = "|latitude|", ylab = "temperature")
  }
  
}
