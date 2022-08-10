nChains = 3
nIter = 1000
nThin = 1
obsmat = NULL
distrmat = NULL
coeff_inits = NULL
sdy_init = NULL
yest_inits = NULL
sdyest_inits = NULL
prior_fun = prior_fun
proposal_var_inits = c(2,2,2,0.2)
adapt_sd = floor(0.1 * nIter)
adapt_sd_decay = max(floor(0.01*nIter),1)
start_adapt = 101
quiet = FALSE




A_sdy = 1000# parameter for the prior on the inverse gamma distribution of sdy
B_sdy = 1000 # parameter for the prior on the inverse gamma distribution of sdy

nbin <- length(unique(obsmat$sample)) + nrow(distrmat)
nbin = 46131

shape_sdy <- A_sdy+nbin/2 # shape parameter for the inverse gamma


lats <- c(as.numeric(sapply(unique(obsmat$sample), function(x) unique(obsmat$p_lat[which(obsmat$sample == x)]))),
          distrmat$p_lat)
temp_mu_emp <- c(as.numeric(sapply(unique(obsmat$sample), function(x) mean(obsmat$p_lat[which(obsmat$sample == x)]))),
                 distrmat$p_lat)

temp_mu <- apply(mode_s10[[1]]$yestimate[1001:5000,],2,median)

temp_pred <- gradient(lats,apply(test1$params[1001:7000,1:4],2,median),0)

#plot(lats,temp_mu)
#points(lats,temp_pred, col = "red")

deviation <- temp_mu_emp-temp_pred
deviation = nbin
hist(sqrt(1/rgamma(10000,
              shape_sdy,
              B_sdy+0.5*sum(deviation^2))),100)

alpha <- 20+nbin/2
beta <- 20+0.5*sum(5*nbin)
beta^2 /((alpha-1)^2*(alpha-2))
beta/(alpha-1)

mean(sqrt(1/rgamma(10000,
                   shape = shape_sdy,
                   rate = B_sdy+0.5*sum(5*nbin))),100)
hist(sqrt(MCMCpack::rinvgamma(10000, shape_sdy, 1/(B_sdy+0.5*sum(5*nbin)))),100)


modm <- readRDS("results/modern/modern_sample_gradient.RDS")
source("R/functions/model_processing/temp_from_gradient.R")
ind <- runif(200,1,nrow(modm))
modt <- temp_from_gradient(lat[ind],modm)
mean((temp[ind]-modt$median)^2)

plot(lat[ind],modt$median)
points(lat[ind],temp[ind], col = "red",pch=2)
mods_all <- combine_posterior(mods,5000)

