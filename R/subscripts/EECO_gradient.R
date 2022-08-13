
### Hollis data
# 
 dat <- readRDS("data/processed/Hollis_processed_EECO_2022_07_19.rds")
# # select data from one stage to test, exclude NA data
# data_sub <- subset(dat,EECO == 1 & !(is.na(temperature)) & !(is.na(paleolat_Meredith)) &
#                      (is.na(depth_habitat) | depth_habitat %in% c("Mixed-layer",     "Mixed layer")))
# 
# data_sub <- subset(data_sub,preservation != "recrystallised" | proxy !="d18O")
# dat <- data_sub
# dat$p_lat <- dat$paleolat_Meredith
# prepare for use in the model
obsmat <- data.frame(sample = (paste(abs(dat$p_lat),dat$longitude, dat$proxy)),
                     p_lat = abs(dat$p_lat), temperature = dat$temperature,
                     sd = dat$temperature_sd,
                     proxy = dat$proxy)

# northern hemisphere:
obsmat <- subset(obsmat,p_lat > 0)

# southern hemisphere:
obsmat <- subset(obsmat,p_lat < 0)
obsmat$p_lat = abs(obsmat$p_lat)

#plot(obsmat$p_lat,obsmat$temperature,xlim = c(-90,90), col = rgb(0,0,0,0.33))

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

#points(distrmat$p_lat,distrmat$mu,pch = 19, col = rgb(0,0.8,0.2,0.75))
### Define the priors

xval <- list(seq(-5,35,0.1),
             seq(-5,60,0.1),
             seq(0,90,0.1),
             seq(0,0.3,0.002))
source("R/subscripts/AuxiliaryFunctions.R")
# plot_prior(priors,xval)

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
mode_south2 <- foreach(pc = 1:nChains) %dopar% {
  # call model functions
  source("R/subscripts/models/EECO_climate_model.R")
  # set random seed
  set.seed(pc)
  hierarchical_model(n_iter = 25000, n_thin = 5,
                  obsmat = obsmat, distrmat = NULL, 
                  logprior_input = priors,adapt_sd = 3000,
                  adapt_sd_decay = 100)
}

### stop cluster
parallel::stopCluster(cl)


#######################################
#
# Assess output

plot_chains(mode)
plot_chains(mode_north)
plot_chains(mode_south)

mode_all <- combine_posterior(mode,5000)
mode_all_north<- combine_posterior(mode_north,4000)
mode_all_south2<- combine_posterior(mode_south2,4000)

mode_all_nobio <- combine_posterior(mode_nobio,4000)
mode_all_nosd <- combine_posterior(mode_nosd,4000)


mcmcse::multiESS(mode_all[,1:4])
plot_gradient(mode_all,ylim = c(13,39))
plot_gradient(mode_all_north,ylim = c(10,37), confint_col = rgb(0,.7,0.5,0.2),line_col = rgb(0,.7,0.5,1), add = F)
plot_gradient(mode_all_south,ylim = c(23,40), confint_col = rgb(0.7,0,0.5,0.2),line_col = rgb(0.7,0,0.5,1), add = F)

plot_posterior(mode[[1]])
plot_posterior(mode_north[[3]], col_obs = rgb(0,.7,0,0.75), col_dist = rgb(0,0,1,0.75))
plot_posterior(mode_south[[3]], col_obs = rgb(1,.2,0,0.75), col_dist = rgb(0,0,1,0.75))

mean_obs_t <-sapply(unique(obsmat$sample), function(x) mean(obsmat$temperature[which(obsmat$sample==x)]))
quant_obs_t <-sapply(unique(obsmat$sample), function(x) quantile(obsmat$temperature[which(obsmat$sample==x)],probs=c(0.025,0.975)))

obs_lat <-sapply(unique(obsmat$sample), function(x) mean(obsmat$p_lat[which(obsmat$sample==x)]))

plot_distr(distrmat)

points(distrmat$p_lat, distrmat$mu, col = rgb(.5,0,.5,.75))
points(obs_lat, mean_obs_t, col = rgb(0,0.75,0,.75), pch = 19)
for(i in 1:length(obs_lat)) points(rep(obs_lat[i],2),quant_obs_t[,i],type="l",col=rgb(0,0.75,0,.75))
points(obs_lat, mean_obs_t, col = rgb(.75,0,0,.75), pch = 10)

legend("topright",c("wide prior on sdy", "narrow prior on sdy"))