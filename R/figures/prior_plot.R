### Define the priors
source("R/options.R")
source("R/functions/plotting/plot_dens.R")
source("R/functions/model_components/dsnorm.R")
source("R/functions/model_components/dtnorm.R")

# plot all priors
xval = list(seq(-6,36,0.01),
             seq(-5,65,0.01),
             seq(9,80,0.1),
             seq(0,0.4,0.001)
)

png("figures/priors.png",res = 600,width = 4,height = 3.8,units = "in")
par(mfrow=c(2,2), mar = c(4.1,4.1,1,1), mgp = c(2,0.7,0),las = 1)
log=F
col = "#67a9cf"
dens <- priors$f1(xval[[1]],log=log)
plot_dens(xval[[1]],dens,xlab = "A", yaxt = "n", xaxt = "n", ylab = "", col = col)
axis(2,c(0,0.5,1),c(0,0.5,1))
axis(1,c(-15,0,15,30,45))

dens <- priors$f2(xval[[2]],xval[[1]][which.max(dens)],log=log)
plot_dens(xval[[2]],dens,xlab = "K", yaxt = "n", xaxt = "n", ylab = "", col = col)
axis(2,c(0,0.5,1),c(0,0.5,1))
axis(1,c(-25,0,25,50,75))


dens <- priors$f3(xval[[3]],log=log)
plot_dens(xval[[3]],dens,xlab = "M", yaxt = "n", xaxt = "n", ylab = "", col = col)
axis(2,c(0,0.5,1),c(0,0.5,1))
axis(1,c(-15,0,15,45,75,105))

dens <- priors$f4(xval[[4]],log=log)
plot_dens(xval[[4]],dens,xlab = "B", yaxt = "n", xaxt = "n", ylab = "", col = col)
axis(2,c(0,0.5,1),c(0,0.5,1))
axis(1,c(-1,0,0.2,0.4,1),c(NA,0,0.2,0.4,NA))
par(mfrow=c(1,1))
dev.off()