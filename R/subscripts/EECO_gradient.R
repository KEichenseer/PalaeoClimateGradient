
### Hollis data

dat <- readRDS("data/processed/Hollis_processed_EECO_2022_07_19.rds")

# prepare for use in the model
obsmat <- data.frame(sample = (paste(abs(dat$p_lat),dat$longitude, dat$proxy)),
                     p_lat = abs(dat$p_lat), temperature = dat$temperature,
                     sd = dat$temperature_sd,
                     proxy = dat$proxy)

plot(obsmat$p_lat,obsmat$temperature,xlim = c(-90,90), col = rgb(0,0,0,0.33))

### Mangrove and Coral data
# read and assign proxy type
bioprox <- readRDS("data/processed/bio_proxies_2022_08_08.RDS")
bioprox$proxy <- rep("Reef",nrow(bioprox))
bioprox$proxy[which(bioprox$taxa=="Avicennia")] <- "Avicennia"
bioprox$proxy[which(bioprox$type=="mangrove" & bioprox$taxa!="Avicennia")] <- "Avicennia-Rhizophoraceae"

# define proxy distributions
proxy_distributions <- data.frame(name = c("Avicennia", "Avicennia-Rhizophoraceae", "Reef"),
               distribution = c("normal", "normal", "normal"),
                mean = c(mean(c(15.6,22.5)), mean(c(20.7,29.5)), 27.6),
                sd = c((22.5-15.6)/4, c(29.5-20.7)/4, (29.5-21)/4),
                shape = rep(NA,3))

# create distribution matrix for use with model
proxy_index <- sapply(bioprox$proxy, function(f) which(proxy_distributions$name==f))

distrmat = data.frame(p_lat = abs(bioprox$p_lat), 
                      mu = proxy_distributions$mean[proxy_index],
                      scale = proxy_distributions$sd[proxy_index],
                      shape = proxy_distributions$shape[proxy_index],
                      distribution = proxy_distributions$distribution[proxy_index])

points(distrmat$p_lat,distrmat$mu,pch = 19, col = rgb(0,0.8,0.2,0.75))
### Define the priors
prior_fun <- readRDS("results/modern/prior_from_modern_gradient.RDS")

xval <- list(seq(-5,35,0.1),
             seq(-5,60,0.1),
             seq(0,90,0.1),
             seq(0,0.6,0.01))
source("R/subscripts/AuxiliaryFunctions.R")
plot_prior(prior_fun,xval)
### Prepare for model run
# Source model script
source("R/subscripts/models/EECO_climate_model.R")

# set up clusters for parallel computing
library(doParallel)
nClust <- 4
cl <- parallel::makeCluster(nClust)
doParallel::registerDoParallel(cl)

nChains = nClust
### Run model
# source script
mode20 <- foreach(pc = 1:nChains) %dopar% {
  # call model functions
  source("R/subscripts/models/EECO_climate_model.R")
  # set random seed
  set.seed(pc)
  test <- hierarchical_model(n_iter = 35000, n_thin = 5,
                  obsmat = obsmat, distrmat = distrmat, 
                  logprior_input = prior_fun,adapt_sd = 3000,
                  adapt_sd_decay = 100, A_sdy = 30, B_sdy = 30)
}

### stop cluster
parallel::stopCluster(cl)


#######################################
#
# Assess output

plot_chains(test)
plot(test$params$sdy)

mode_all <- combine_posterior(mode,4000)
mode_all_s20 <- combine_posterior(mode20,4000)
mode_all_s10 <- combine_posterior(mode_s10,4000)

mode_all_nobio <- combine_posterior(mode_nobio,4000)
mode_all_nosd <- combine_posterior(mode_nosd,4000)


mcmcse::multiESS(mode_all_s10[,1:4])
plot_gradient(test,ylim = c(13,39))
plot_gradient(mode_all_s3,ylim = c(10,37), confint_col = rgb(0,.7,0.5,0.2),line_col = rgb(0,.7,0.5,1), add = T)
plot_gradient(mode_all_s10,ylim = c(10,37), confint_col = rgb(0.7,0,0.5,0.2),line_col = rgb(0.7,0,0.5,1), add = T)

plot_posterior(test)
plot_posterior(mode_s10[[3]], col_obs = rgb(.5,1,0,0.75), col_dist = rgb(0,0,1,0.75))

points(distrmat$p_lat, distrmat$mu, col = rgb(.5,0,.5,.75))
