### Assess EECO model output
# load figures for plotting, distributions ect.
source("./R/functions/auxiliary_functions.R")
# read model output
mod <- readRDS("./results/eeco/eeco_climate_model_output.rds")
# read the input data for comparison
obsmat <- readRDS("./data/processed/observation_matrix_for_EECO_model.rds")
distrmat <- readRDS("./data/processed/distribution_matrix_for_EECO_model.rds")
# calculate location means of observation input data for easier plotting
ylist <- lapply(unique(obsmat$sample), function(x) obsmat$temperature[which(obsmat$sample == x)])
yobs_mean= c(sapply(ylist,mean))
n_obs <- length(yobs_mean)
# add distribution means
mean_t <- c(yobs_mean,distrmat$mu)
lats <- mod[[1]]$lat # latitudes of locations is saved in model output
  

### Assess output
# monitor the chains to check for convergence and mixing
plot_chains(mod) # default is A, dKA, M and log B
# plot chains of sdy and lof posterior
plot_chains(mod, params = c("sdy", "logpost"))

# combine chains, discard burnin
mode_all <- combine_posterior(mod,100000)
saveRDS(mode_all,"results/eeco/eeco_climate_model_output_params.rds")
# check effective sample size for A, dKA, M, and B
mcmcse::multiESS(mode_all[,1:4])

# Plot the gradient
plot_gradient(mode_all,ylim = c(-2,39))

# add median of the poster estimates of location means
plot_posterior(mod[[1]]) # select the first run, doesn't yet work for combined chains

## add the empirical locality means 
# colour by type 
cols = rep(rgb(.85,0,0,0.5),length(lats)) 
cols[(n_obs+1):length(lats)] <- rgb(0,0,.85,0.5) # set distribution colour to blue
# plot empirical means
points(lats,mean_t,pch = 1, col = cols, bg = NA, lwd = 2)

legend("topright",c("posterior mean (proxy obs.)", "posterior mean (distribution)",
                    "empirical mean (proxy obs.)", "assigned distribution mean"),
       pch = c(17,17,1,1), col = c(rgb(0.8,0.5,0,0.75),rgb(0,0.5,.8,0.75),
                                   rgb(0.85,0,0,0.75),rgb(0,0,.85,0.75)))

# add modern gradient:
# source("./R/functions/model_processing/temp_from_gradient.R")
# modern_sample <- readRDS("./results/modern/modern_climate_model_output.rds")
# points(0:90,apply(temp_from_gradient(0:90,modern_sample),1,median), type = "l", lwd = 2)

