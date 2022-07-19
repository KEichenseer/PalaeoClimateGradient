dat <- readRDS("data/processed/Hollis_processed_2022_07_19.rds")
hist(dat$temperature,100)

### generate obsmat


  # select data from one stage to test, exclude NA data
  data_sub <- subset(dat,EECO == 1 & !(is.na(temperature)) & !(is.na(paleolat_Meredith)) &
                       (is.na(depth_habitat) | depth_habitat %in% c("Mixed-layer",     "Mixed layer")))
  data_sub <- data_sub[with(data_sub, order(abs(paleolat_Meredith), longitude)),]
  
  
  plot(data_sub$latitude,data_sub$temperature,xlim=c(-90,90))
  points(data_sub$latitude[which(data_sub$preservation %in% c("recrystallised","Recrystallized"))],
         data_sub$temperature[which(data_sub$preservation %in% c("recrystallised","Recrystallized"))],
         col = NA, pch = 21, bg = rgb(1,0,0,0.5))
  
  # prepare for use in the model
  obsmat <- data.frame(sample = (paste(abs(data_sub$paleolat_Meredith),data_sub$longitude)),
                         latitude = abs(data_sub$paleolat_Meredith), temperature = data_sub$temperature,
                         sd = data_sub$temperature_sd)


  distrmat = NULL
  
  nIter = 50000
  obsmat = obsmat
  proposal_var_inits = c(2,2,2,0.2)
  adapt_sd = floor(0.1 * nIter)
  adapt_sd_decay = max(floor(0.01*nIter),1)
  quiet = FALSE
  yest_inits = NULL
  sdyest_inits = NULL
  # set random seed
  set.seed(1)
  
  # random setting of initial values for the regression parameters
  coeff_inits = rep(NA,4)
  coeff_inits[1] = rnorm(1,10,3) # 20,45),c(1,2,4.5)
  coeff_inits[2] = coeff_inits[1] + truncnorm::rtruncnorm(1,0,Inf,12,6)
  coeff_inits[3] = rnorm(1,45,7.5)
  coeff_inits[4] = exp(rnorm(1,-2.3,0.25))
  sdy_init = exp(rnorm(1,0.7,0.25))
  
  # deterministic setting of initial values for the temperature and temperature sd estimates
  if((!(is.null(distrmat)) | !(is.null(obsmat))) & is.null(yest_inits)) {
    yest_inits <- c(unlist(sapply(unique(obsmat$sample), function(x) mean(obsmat$temperature[which(obsmat$sample == x)]))),
                    distrmat$location + distrmat$scale * sqrt(2/pi) * 
                      distrmat$shape/sqrt(1+distrmat$shape^2) )
  }
  if(!(is.null(obsmat)) & is.null(sdyest_inits)) sdyest_inits <- rep(2,length(unique(obsmat$sample)))
  
  source("R/subscripts/ClimateGradientModelwithSDonObs.R")
  source("R/subscripts/AuxiliaryFunctions.R")
  
  modh1 <- run_MCMC_sd_obs(nIter = nIter, obsmat = obsmat, distrmat = NULL, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                               proposal_var_inits = c(2,2,2,0.2), adapt_sd = floor(0.2 * nIter),
                               adapt_sd_decay = max(floor(0.005*nIter),1), quiet = FALSE)
  # only mixed layer d18O
  modh3 <- run_MCMC_sd_obs(nIter = nIter, obsmat = obsmat, distrmat = NULL, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                           proposal_var_inits = c(2,2,2,0.2), adapt_sd = floor(0.2 * nIter),
                           adapt_sd_decay = max(floor(0.005*nIter),1), quiet = FALSE)
  plot(modh1$params$sdy)
plot_gradient(modh1, ylim = c(2,42))  
points(obsmat$latitude,obsmat$temperature, col = NA, bg = rgb(0,0.4,0.5,0.33), pch = 21, cex = 0.8)
plot_posterior(modh1)

plot_gradient(modh2, ylim = c(2,42), add = T, line_col = rgb(0,0.8,0,1), confint_col = rgb(0,0.8,0,0.2))  
points(obsmat$latitude,obsmat$temperature, col = NA, bg = rgb(0,0.4,0.5,0.33), pch = 21, cex = 0.8)
plot_posterior(modh2, col_obs  = rgb(0,1,0,0.5))


plot_gradient(modh3, ylim = c(2,42), add = T, line_col = rgb(0.8,0,0,1), confint_col = rgb(0.8,0,0,0.2))  
plot_posterior(modh3, col_obs  = rgb(1,0,0,0.75))

plot_sample_gradient(modh1, burnin = 25000)

graddat <- plot_sample_gradient(modh2, burnin = 25000, n_samples = 5000, return_data = T, plot = FALSE)

graddiff <- apply(graddat,1,function(x) abs(diff(range(x))))
hist(graddiff,seq(0,40,1))
abline(v=quantile(graddiff,probs = c(0.025,0.5,0.975)),col = "red", lwd = c(2,3,2), lty=c(3,2,3))
