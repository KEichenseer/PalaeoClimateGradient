# compare gradient sdy
source("./R/functions/model_processing/combine_posterior.R")
source("./R/functions/model_processing/temp_from_gradient.R")


modern_gradient_eeco_lat <- readRDS("./results/modern/modern_sample_eocene_p_lat_gradient.RDS")
modern_gradient_eeco_lat <- combine_posterior(mod = modern_gradient_eeco_lat, burnin = 5000)
eeco_gradient <- readRDS("./results/eeco/eeco_climate_model_output_combined.rds")
modern_gradient<- readRDS("results/modern/modern_climate_model_output.rds")

modern_eeco_sample_sdy <- as.list(format(round(quantile(modern_gradient_eeco_lat$sdy,probs = c(.025,.5,.975)),1),nsmall = 1))

eeco_sdy <- as.list(format(round(quantile(eeco_gradient$sdy,probs = c(.025,.5,.975)),1),nsmall = 1))

modern_sdy <- as.list(format(round(quantile(modern_gradient$sdy,probs = c(.025,.5,.975)),1),nsmall = 1))


saveRDS(modern_eeco_sample_sdy, "results/modern/modern_eeco_sample_sdy.rds")
saveRDS(eeco_sdy, "results/eeco/eeco_sdy.rds")
saveRDS(modern_sdy, "results/modern/modern_sdy.rds")

# just proxy sdy
eeco_prox <- readRDS("results/eeco/eeco_climate_model_output_just_proxy.rds")
eeco_prox <- combine_posterior(eeco_prox, burnin = 100000)
eeco_prox_sdy <- as.list(format(round(quantile(eeco_prox$sdy,probs = c(.025,.5,.975)),1),nsmall = 1))

# residual sum of squares modern gradient modeled vs zonal means
mod_emp_med <- readRDS("results/modern/empirical_median_1-deg-lat.rds")

lat <- seq(0.5,89.5,1)
mod_med <- (temp_from_gradient(lat,modern_gradient))$median

# coefficient of determination
R2_modern <- 1-sum((mod_emp_med$SST - mod_med)^2)/sum((mod_emp_med$SST-mean(mod_emp_med$SST))^2)
R2_modern <- format(round(R2_modern,3),nsmall = 1)
saveRDS(R2_modern,"results/modern/modern_model_empirical_R2.rds")

# Gelman's R2
lats <- seq(0.5,89.5,1)

mod_t <- gradient(lats,modern_gradient)

var_fit <- apply(mod_t,1,function(x) var(x))
var_res <- apply(mod_t,1,function(x) var(x-mod_emp_med$SST))

R_squared_gelman_modern <- var_fit/(var_fit+var_res)
quantile(R_squared_gelman_modern, probs = c(0.025,0.5,0.975))
saveRDS(format(round(median(R_squared_gelman_modern),3),nsmall = 1) ,"results/modern/modern_model_empirical_R2_gelman.rds")

# modern empirical gradient
saveRDS(format(round(mod_emp_med$SST[1]-mod_emp_med$SST[90],1),nsmall = 1),"results/modern/modern_empirical_gradient.rds")
