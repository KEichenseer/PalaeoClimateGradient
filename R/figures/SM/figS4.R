### Supplementary figure: different gradient shapes
source("R/functions/model_components/climate_model_modern.R")
source("R/functions/model_components/gradient.R")
source("R/functions/model_processing/combine_posterior.R")
source("R/functions/auxiliary_functions.R")

# source options for analysis
source("R/options.R")
#

ylim = c(-2,32)

## 1. Linear
set.seed(1)
lat <- runif(30,0,90)
temp <- 30-lat/3
m1 <- run_MCMC_simple(n_iter = 20000,
                n_thin = 10,
                x = lat,
                y = temp,
                prior_input = priors,
                adapt_sd = NULL)
r1 <- combine_posterior(list(m1))
plot_chains(list(m1))
par(mfrow=c(1,1))
plot(lat,temp, xlab = expression("|latitude| ("*degree*")"), 
     ylab = expression("temperature ("*degree*"C)"),
     pch = 21, col = NA, bg = rgb(0,0,0,0.5), ylim = ylim)
points(c(0,90),c(30,0),type = "l", lwd = 2)
plot_gradient(r1,add = T, line_col = "red", confint_col = rgb(1,0,0,0.2))

## 2. Flat
set.seed(1)
lat <- runif(30,0,90)
temp <- rep(25,30)
m1 <- run_MCMC_simple(n_iter = 20000,
                      n_thin = 10,
                      x = lat,
                      y = temp,
                      prior_input = priors,
                      adapt_sd = NULL)
r1 <- combine_posterior(list(m1))
plot_chains(list(m1))
par(mfrow=c(1,1))
plot(lat,temp, xlab = expression("|latitude| ("*degree*")"), 
     ylab = expression("temperature ("*degree*"C)"),
     pch = 21, col = NA, bg = rgb(0,0,0,0.5), ylim = c(23,27))
points(c(0,90),c(25,25),type = "l", lwd = 2)
plot_gradient(r1,add = T, line_col = "red", confint_col = rgb(1,0,0,0.2))

## 3. Quadratic
set.seed(1)
lat <- runif(30,0,90)
temp <- 30-lat^2/270
m1 <- run_MCMC_simple(n_iter = 20000,
                      n_thin = 10,
                      x = lat,
                      y = temp,
                      prior_input = priors,
                      adapt_sd = NULL)
r1 <- combine_posterior(list(m1))
plot_chains(list(m1))
par(mfrow=c(1,1))
plot(lat,temp, xlab = expression("|latitude| ("*degree*")"), 
     ylab = expression("temperature ("*degree*"C)"),
     pch = 21, col = NA, bg = rgb(0,0,0,0.5), ylim = ylim)
true_grad <- 30- seq(0,90,0.2)^2/270
points(seq(0,90,0.2),true_grad,type = "l", lwd = 2)
plot_gradient(r1,add = T, line_col = "red", confint_col = rgb(1,0,0,0.2))

## 4. Inverse Sigmoidal
set.seed(1)
lat <- runif(30,0,90)
temp <- 15-3*log(lat/(90-lat))
plot(lat,temp)
m1 <- run_MCMC_simple(n_iter = 20000,
                      n_thin = 10,
                      x = lat,
                      y = temp,
                      prior_input = priors,
                      adapt_sd = NULL)
r1 <- combine_posterior(list(m1))
plot_chains(list(m1))
par(mfrow=c(1,1))
plot(lat,temp, xlab = expression("|latitude| ("*degree*")"), 
     ylab = expression("temperature ("*degree*"C)"),
     pch = 21, col = NA, bg = rgb(0,0,0,0.5), ylim = ylim)
true_grad <- 30- seq(0,90,0.2)^2/270
points(seq(0,90,0.2),true_grad,type = "l", lwd = 2)
plot_gradient(r1,add = T, line_col = "red", confint_col = rgb(1,0,0,0.2))
