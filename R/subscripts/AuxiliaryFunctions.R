### auxiliary functions

error_polygon <- function(x,en,ep,color=rgb(0,0,0,0.2)) {
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
}

plot_gradient <- function(model_out, burnin = NULL, lat = seq(0,90,0.2), confint_n = NULL) {
  
  nIter <- nrow(model_out)
  if(is.null(burnin)) burnin = round(nrow(model_out)*12/100)+1
  if(is.null(confint_n)) confint_n = min(5000,length(burnin:nIter))
    
  sample_it <- sample(burnin:nIter,confint_n)
  grad_q <- apply(gradient(lat,model_out[sample_it,1:4],0),2,function(a) 
    quantile(a,probs = c(0.025,0.975)))

  med_grad <- gradient(lat,apply(model_out[burnin:nrow(model_out),1:4],2,median),0)
  
  ylims <- c(max(min(med_grad*0.75),min(grad_q["2.5%",])),min(max(med_grad*1.25),max(grad_q["97.5%",])))
  
  
  plot(lat, med_grad, ylim = ylims,
       type = "l", lwd = 2, xlab = "|latitude|", ylab = "temperature")
  error_polygon(lat,grad_q["2.5%",],grad_q["97.5%",])
}
