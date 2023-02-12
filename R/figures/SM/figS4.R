### Supplementary figure: different gradient shapes
source("R/functions/model_components/climate_model_modern.R")
source("R/functions/model_components/gradient.R")
source("R/functions/model_processing/combine_posterior.R")
source("R/functions/auxiliary_functions.R")

# source options for analysis
source("R/options.R")
#

ylim = c(0,30)
xlim = c(0,90)
lat <- abs(c(-59.42, -59.10, -54.13, -47.84, -47.84, -47.47, -46.14, -45.80, -45.80, -24.25,  -8.57,  -2.05,  -1.86,   0.49,   4.29,  14.67,
  27.44,  31.00,  32.50,  32.80,  32.90,  35.54,  36.22,  39.02,  40.39,  40.63,  40.89,  43.79,  44.11,  50.04,  50.39,  58.45,
  76.45,  77.13))

## 1. Linear
set.seed(1)
#lat <- runif(30,0,90)
temp1 <- 30-lat/3
m1 <- run_MCMC_simple(n_iter = 20000,
                n_thin = 10,
                x = lat,
                y = temp1,
                prior_input = priors,
                adapt_sd = NULL)
r1 <- combine_posterior(list(m1))
plot_chains(list(m1))


## 2. Flat
set.seed(1)
#lat <- runif(30,0,90)
temp2 <- rep(25,length(lat))
m2 <- run_MCMC_simple(n_iter = 20000,
                      n_thin = 10,
                      x = lat,
                      y = temp2,
                      prior_input = priors,
                      adapt_sd = NULL)
r2 <- combine_posterior(list(m2))
plot_chains(list(m1))

## 3. Quadratic
set.seed(1)
#lat <- runif(30,0,90)
temp3 <- 30-lat^2/270
m3 <- run_MCMC_simple(n_iter = 20000,
                      n_thin = 10,
                      x = lat,
                      y = temp3,
                      prior_input = priors,
                      adapt_sd = NULL)
r3 <- combine_posterior(list(m3))
plot_chains(list(m1))

## 4. Quadratic, relaxed A prior
# relax prior on A
priors$f1 <- function(x,log) dnorm(x,0,12,log)
priors$f3 <- function(x,log) dnorm(x,45,25,log)
priors$f4 <- function(x,log) dgamma(x,shape = 3,rate = 21,log=log)
set.seed(1)
#lat <- runif(30,0,90)
temp4 <- 30-lat^2/270
m4 <- run_MCMC_simple(n_iter = 20000,
                      n_thin = 10,
                      x = lat,
                      y = temp4,
                      prior_input = priors,
                      adapt_sd = NULL)
r4 <- combine_posterior(list(m4))
plot_chains(list(m1))


# ## 4. Inverse Sigmoidal
# set.seed(1)
# #lat <- runif(30,0,90)
# temp4 <- 14-3*log(lat/(90-lat))
# plot(lat,temp4)
# m4 <- run_MCMC_simple(n_iter = 20000,
#                       n_thin = 10,
#                       x = lat,
#                       y = temp4,
#                       prior_input = priors,
#                       adapt_sd = NULL)
# r4 <- combine_posterior(list(m4))
# plot_chains(list(m4))


png("figures/SM/FigS4_different_gradients.png",width = 5,height = 4.5,unit = "in", res = 300)
par(mfrow=c(2,2), mar = c(4,3.75,.75,1), mgp = c(2,0.6,0), las = 1)

plot(lat,temp1, xlab = expression("|latitude| ("*degree*")"), 
     ylab = expression("temperature ("*degree*"C)"),
     pch = 21, col = NA, bg = rgb(0,0,0,0.5), ylim = ylim, xlim = xlim, cex = 1.4)
points(c(0,90),c(30,0),type = "l", lwd = 2)
plot_gradient(r1,add = T, line_col = "red", confint_col = rgb(0.95,0,0,0.33))
mtext("a",side=2,line=3,padj=-4.7,adj=0,cex=1.2)

plot(lat,temp2, xlab = expression("|latitude| ("*degree*")"), 
     ylab = expression("temperature ("*degree*"C)"),
     pch = 21, col = NA, bg = rgb(0,0,0,0.5), ylim = c(23,27), xlim = xlim, cex = 1.4)
points(c(0,90),c(25,25),type = "l", lwd = 2)
plot_gradient(r2,add = T, line_col = "red", confint_col = rgb(0.95,0,0,0.33))
mtext("b",side=2,line=3,padj=-4.7,adj=0,cex=1.2)

plot(lat,temp3, xlab = expression("|latitude| ("*degree*")"), 
     ylab = expression("temperature ("*degree*"C)"),
     pch = 21, col = NA, bg = rgb(0,0,0,0.5), ylim = ylim, xlim = xlim, cex = 1.4)
true_grad <- 30- seq(0,90,0.2)^2/270
points(seq(0,90,0.2),true_grad,type = "l", lwd = 2)
plot_gradient(r3,add = T, line_col = "red", confint_col = rgb(0.95,0,0,0.33))
mtext("c",side=2,line=3,padj=-4.7,adj=0,cex=1.2)

plot(lat,temp4, xlab = expression("|latitude| ("*degree*")"), 
     ylab = expression("temperature ("*degree*"C)"),
     pch = 21, col = NA, bg = rgb(0,0,0,0.5), ylim = ylim, xlim = xlim, cex = 1.4)
#true_grad <- 14-3*log(seq(0,90,0.2)/(90-seq(0,90,0.2)))
true_grad <- 30- seq(0,90,0.2)^2/270
points(seq(0,90,0.2),true_grad,type = "l", lwd = 2)
plot_gradient(r4,add = T, line_col = "red", confint_col = rgb(0.95,0,0,0.33))
mtext("d",side=2,line=3,padj=-4.7,adj=0,cex=1.2)
dev.off()


## cos test

lat <- seq(-90,90,1)
a = 10
b = -0.03
c = 0.0003
d = 5
t1 <- a + b*lat + d*cos(c*lat^2)
plot(lat,t1, type = "l")
