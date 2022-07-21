dat <- readRDS("data/processed/Hollis_processed_2022_07_19.rds")
dat$proxy_value <- as.numeric(dat$proxy_value)
hist(dat$temperature,100)


dat_sub <- subset(dat,proxy=="d18O")
dat_sub$mineral = "calcite"
dat_sub$fossil_group = "foraminifera"
dat_sub$d18O_sw_ice <- -1.08

grosman_t <- get_temperature(d18O_sample = dat_sub$proxy_value,
                             d18O_sw_ice = dat_sub$d18O_sw_ice, lat = abs(dat_sub$palaeolat_ori),
                             mineral = dat_sub$mineral, fossil_group = dat_sub$organism)

plot(grosman_t,dat_sub$temperature)
     hist(dat$proxy_value)
     abline(lm(dat_sub$temperature~grosman_t), col = "red")
     
     dat$grosman_t <- NA
     dat$grosman_t[which(dat$proxy=="d18O")] <- grosman_t
     
  mean(grosman_t-dat_sub$temperature)   
     
  ### generate obsmat


  # select data from one stage to test, exclude NA data
  data_sub <- subset(dat,EECO == 1 & !(is.na(temperature)) & !(is.na(paleolat_Meredith)) &
                       (is.na(depth_habitat) | depth_habitat %in% c("Mixed-layer",     "Mixed layer")))
  data_sub <- data_sub[with(data_sub, order(abs(paleolat_Meredith), longitude)),]
  
  data_sub <- subset(data_sub,proxy=="d18O")
  plot(data_sub$latitude,data_sub$temperature,xlim=c(-90,90), pch = 17, col = rgb(1,0,0,0.5), ylim = c(0,45))
  points(data_sub$latitude[which(data_sub$preservation %in% c("recrystallised","Recrystallized"))],
         data_sub$temperature[which(data_sub$preservation %in% c("recrystallised","Recrystallized"))],
         col = NA, pch = 21, bg = rgb(1,0,0,0.5))
  points(data_sub$latitude,
         data_sub$grosman_t,
         col = NA, pch = 21, bg = rgb(0,1,0,0.5))
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
  
  mod_nod18O <- run_MCMC_sd_obs(nIter = nIter, obsmat = obsmat, distrmat = NULL, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                               proposal_var_inits = c(2,2,2,0.2), adapt_sd = floor(0.2 * nIter),
                               adapt_sd_decay = max(floor(0.005*nIter),1), quiet = FALSE)
  # only mixed layer d18O
  mod_dO18 <- run_MCMC_sd_obs(nIter = nIter, obsmat = obsmat, distrmat = NULL, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                           proposal_var_inits = c(2,2,2,0.2), adapt_sd = floor(0.2 * nIter),
                           adapt_sd_decay = max(floor(0.005*nIter),1), quiet = FALSE)
  plot(modh1$params$sdy)
plot_gradient(mod_nod18O, ylim = c(2,42), line_col = rgb(0,0.8,0,1), confint_col = rgb(0,0.8,0,0.2) )  
points(obsmat2$latitude,obsmat2$temperature, col = NA, bg = rgb(0,0.7,0.2,0.33), pch = 21, cex = 1)
plot_posterior(mod_nod18O,col_obs = rgb(0,0.5,0.8,0.7))

plot_gradient(mod_dO18, ylim = c(2,42), line_col = rgb(0.8,0,0.1,1), confint_col = rgb(0.8,0,0.1,0.175), add = T)  
points(obsmat$latitude,obsmat$temperature, col = NA, bg = rgb(0.7,0,0.1,0.33), pch = 22, cex = 1)
plot_posterior(mod_dO18,col_obs = rgb(0.85,0.75,0,0.7),cex = 1.1)

legend("topright",c("d18O", "MgCa, D47 and TEX86"), lty = c(1,1), lwd = c(2,2), col = c(rgb(0.8,0,0,1),rgb(0,0.8,0,1)))




plot_gradient(modh2, ylim = c(2,42), add = T, line_col = rgb(0,0.8,0,1), confint_col = rgb(0,0.8,0,0.2))  

cols <- c(rgb(1,0,0,0.4),rgb(0.8,0.5,0,0.5),rgb(0,1,0,0.4), rgb(0,0,1,0.4))
pchs = c(21,22,23,24)
colindex <- sapply(data_sub$proxy,function(x) which(x==c("d18O","d47","MgCa","TEX86")))

plot_gradient(modh3, ylim = c(2,42))  
mtext("EECO - Hollis et al. 2019 data, Grosman comparison", cex = 1.2)
points(obsmat$latitude,obsmat$temperature, col = NA, bg = cols[colindex], pch = pchs[colindex], cex = 0.8)
plot_posterior(modh3, col_obs  = rgb(0,0,0,0.35))

legend("topright",c("d18O", "MgCa", "D47", "TEX86", "post. estimate"), pch = c(21,23,22,24,24), 
       col = NA, pt.bg  = c(cols[c(1,3,2,4)],rgb(0,0,0,0.35)),pt.cex = 1.3)


plot_gradient(modh_dO18, ylim = c(2,42), add = T, line_col = rgb(1,0,0,1), confint_col = rgb(1,0,0,0.2))  

plot_gradient(modh3, ylim = c(2,42), add = T, line_col = rgb(0.8,0,0,1), confint_col = rgb(0.8,0,0,0.2))  
plot_posterior(modh3, col_obs  = rgb(1,0,0,0.75))

plot_sample_gradient(modh1, burnin = 25000)

graddat <- plot_sample_gradient(modh2, burnin = 25000, n_samples = 5000, return_data = T, plot = FALSE)

graddiff <- apply(graddat,1,function(x) abs(diff(range(x))))
hist(graddiff,seq(0,40,1))
abline(v=quantile(graddiff,probs = c(0.025,0.5,0.975)),col = "red", lwd = c(2,3,2), lty=c(3,2,3))
