plot_dens <- function(x,dens,xlim = NULL,ylim = NULL,col = rgb(0,0.5,0.7,0.25),
                      xlab = "", ylab = "normalised density", add = FALSE, xaxs = "i",...) {
  dens <- dens/max(dens)
  if(is.null(ylim)) ylim = c(0,max(dens)*1.03)
  if(is.null(xlim)) xlim = range(x)
  
  if(add == FALSE) plot(0,0,type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, xaxs = xaxs, yaxs = "i",...)
  polygon(c(x[1],x,x[length(x)]),c(0,dens,0), col = col, border = NA)
  
}