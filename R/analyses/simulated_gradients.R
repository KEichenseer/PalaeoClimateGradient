### Supplementary figure: different gradient shapes

library(dplyr)

# source functions
source("R/functions/model_components/climate_model_modern.R")
source("R/functions/model_components/gradient.R")
source("R/functions/model_processing/combine_posterior.R")
source("R/functions/auxiliary_functions.R")
# source options for analysis
source("R/options.R")
#

# number of samples (repetitions)
n_rep <- 1000 # careful, with 1000 reps it takes hours to run the whole script
# sample sizes in the repetitions (4 different)
sample_sizes <- c(5,10,20,50)

### Model 4 different gradients
##
## 1. Linear
m_lin <- list(NULL)
for(ss in 1:length(sample_sizes)) {
  m_lin[[ss]] <- list()
  print(paste("sample_size =",sample_sizes[ss]))
  set.seed(1)
  for(rep in 1:n_rep) {
    lat <- 180/pi * asin(runif(sample_sizes[ss],0,1))
    temp <- 30-lat/3
    m_lin[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                          n_thin = 50,
                          x = lat,
                          y = temp,
                          prior_input = priors,
                          adapt_sd = NULL)
    print(paste("repetition =",rep))
    
  }
  
  
}

m_lin_proc <- list()
for(i in 1:length(sample_sizes)) m_lin_proc[[i]] <- combine_posterior(m_lin[[i]],burnin = 5000)

lats <- seq(0,90,0.5)
temp_lin <- list()
lats <- seq(0,90,0.5)
for(i in 1:4) {
  temp <- gradient(lats,m_lin_proc[[i]])
  temp_lin[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

par(mfrow=c(4,4))
for(i in 1:4) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,30), xlab = "abs. latitude", ylab = "temperature")
  points(lats, 30-lats/3, type = "l", lwd = 2)
  #plot_gradient(m_lin_proc[[i]], add = T)
  error_polygon(lats,temp_lin[[i]][2,],temp_lin[[i]][3,], col = rgb(1,0,0,0.10))
  error_polygon(lats,temp_lin[[i]][4,],temp_lin[[i]][5,], col = rgb(1,0,0,0.15))
  error_polygon(lats,temp_lin[[i]][6,],temp_lin[[i]][7,], col = rgb(1,0,0,0.20))
  
  points(lats,temp_lin[[i]][1,], type = "l", lwd = 2, col = rgb(1,0,0,0.75))
}


## 2. Flat
m_flat <- list(NULL)
for(ss in 1:length(sample_sizes)) {
  m_flat[[ss]] <- list()
  print(paste("sample_size =",sample_sizes[ss]))
  set.seed(1)
  for(rep in 1:n_rep) {
    lat <- 180/pi * asin(runif(sample_sizes[ss],0,1))
    temp <- rep(25,length(lat))
    m_flat[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                                          n_thin = 50,
                                          x = lat,
                                          y = temp,
                                          prior_input = priors,
                                          adapt_sd = NULL)
    print(paste("repetition =",rep))
    
  }
  
  
}

m_flat_proc <- list()
for(i in 1:length(sample_sizes)) m_flat_proc[[i]] <- combine_posterior(m_flat[[i]],burnin = 5000)

lats <- seq(0,90,0.5)
temp_flat <- list()
lats <- seq(0,90,0.5)
for(i in 1:4) {
  temp <- gradient(lats,m_flat_proc[[i]])
  temp_flat[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

for(i in 1:4) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,30), xlab = "abs. latitude", ylab = "temperature")
  points(lats, rep(25,length(lats)), type = "l", lwd = 2)
  #plot_gradient(m_flat_proc[[i]], add = T)
  error_polygon(lats,temp_flat[[i]][2,],temp_flat[[i]][3,], col = rgb(1,0,0,0.12))
  error_polygon(lats,temp_flat[[i]][4,],temp_flat[[i]][5,], col = rgb(1,0,0,0.17))
  error_polygon(lats,temp_flat[[i]][6,],temp_flat[[i]][7,], col = rgb(1,0,0,0.20))
  
  points(lats,temp_flat[[i]][1,], type = "l", lwd = 2, col = rgb(1,0,0,0.75))
}


## 3. Quadratic
m_quad <- list(NULL)
for(ss in 1:length(sample_sizes)) {
  m_quad[[ss]] <- list()
  print(paste("sample_size =",sample_sizes[ss]))
  set.seed(1)
  for(rep in 1:n_rep) {
    lat <- 180/pi * asin(runif(sample_sizes[ss],0,1))
    temp <- 30-lat^2/270
    m_quad[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                                           n_thin = 50,
                                           x = lat,
                                           y = temp,
                                           prior_input = priors,
                                           adapt_sd = NULL)
    print(paste("repetition =",rep))
    
  }
  
  
}

m_quad_proc <- list()
for(i in 1:length(sample_sizes)) m_quad_proc[[i]] <- combine_posterior(m_quad[[i]],burnin = 5000)

lats <- seq(0,90,0.5)
temp_quad <- list()
lats <- seq(0,90,0.5)
for(i in 1:4) {
  temp <- gradient(lats,m_quad_proc[[i]])
  temp_quad[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

for(i in 1:4) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,33), xlab = "abs. latitude", ylab = "temperature")
  points(lats, 30-lats^2/270, type = "l", lwd = 2)
  #plot_gradient(m_quad_proc[[i]], add = T)
  error_polygon(lats,temp_quad[[i]][2,],temp_quad[[i]][3,], col = rgb(1,0,0,0.12))
  error_polygon(lats,temp_quad[[i]][4,],temp_quad[[i]][5,], col = rgb(1,0,0,0.17))
  error_polygon(lats,temp_quad[[i]][6,],temp_quad[[i]][7,], col = rgb(1,0,0,0.20))
  
  points(lats,temp_quad[[i]][1,], type = "l", lwd = 2, col = rgb(1,0,0,0.75))
}


## 4. Steep
m_steep <- list(NULL)
for(ss in 1:length(sample_sizes)) {
  m_steep[[ss]] <- list()
  print(paste("sample_size =",sample_sizes[ss]))
  set.seed(1)
  for(rep in 1:n_rep) {
    lat <-180/pi * asin(runif(sample_sizes[ss],0,1))
    temp <- 8.5*atan(-0.5*lat+15)+12.5
    m_steep[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                                           n_thin = 50,
                                           x = lat,
                                           y = temp,
                                           prior_input = priors,
                                           adapt_sd = NULL)
    print(paste("repetition =",rep))
    
  }
  
  
}

m_steep_proc <- list()
for(i in 1:length(sample_sizes)) m_steep_proc[[i]] <- combine_posterior(m_steep[[i]],burnin = 5000)

lats <- seq(0,90,0.5)
temp_steep <- list()
lats <- seq(0,90,0.5)
for(i in 1:4) {
  temp <- gradient(lats,m_steep_proc[[i]])
  temp_steep[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

for(i in 1:4) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,35), xlab = "abs. latitude", ylab = "temperature")
  #plot_gradient(m_steep_proc[[i]], add = T)
  error_polygon(lats,temp_steep[[i]][2,],temp_steep[[i]][3,], col = rgb(1,0,0,0.12))
  error_polygon(lats,temp_steep[[i]][4,],temp_steep[[i]][5,], col = rgb(1,0,0,0.17))
  error_polygon(lats,temp_steep[[i]][6,],temp_steep[[i]][7,], col = rgb(1,0,0,0.20))
  
  #for(j in 1:100) points(lats,gradient(lats,m_steep_proc[[i]][j*500,1:4]), type = "l", col = rgb(0,0.3,0.9,0.2)) 
  points(lats,temp_steep[[i]][1,], type = "l", lwd = 2, col = rgb(1,0,0,0.75))
  points(lats, 8.5*atan(-0.5*lats+15)+12.5, type = "l", lwd = 2, col = "black")
  
  }

m_lin <- four_gradients_model[[1]]
m_flat <- four_gradients_model[[2]]
m_quad <- four_gradients_model[[3]]
m_steep <- four_gradients_model[[4]]

### save list of all temperature outputs
four_gradients_model <- list(m_lin_proc,m_flat_proc,m_quad_proc,m_steep_proc)
four_gradients_temp <- list(temp_lin,temp_flat,temp_quad,temp_steep)

saveRDS(four_gradients_model, "D://OneDrive - Durham University/projects/four_gradients_model_proc.rds")
saveRDS(four_gradients_temp, "D://OneDrive - Durham University/projects/four_gradients_temp.rds")

gradient_formulas <- list( function(x) 30-x/3, function(x) rep(25,length(x)), function(x) 30-x^2/270, function(x) 8.5*atan(-0.5*x+15)+12.5)

### get temperatures from combining all gradients
t_list <- list()
for(i in 1:4) {
  t_list[[i]] <- list()
  for(j in 1:4) {

    grad <- gradient(lats,four_gradients_model[[i]][[j]],0)
    t_list[[i]][[j]] <- data.frame(lat = lats,
                                   median = apply(grad,2,median),
                                   mean = apply(grad,2,mean),
                                   p_005 = apply(grad,2, function(x) quantile(x,0.005)),
                                   p_025 = apply(grad,2, function(x) quantile(x,0.025)),
                                   p_100 = apply(grad,2, function(x) quantile(x,0.1)),
                                   p_900 = apply(grad,2, function(x) quantile(x,0.9)),
                                   p_975 = apply(grad,2, function(x) quantile(x,0.975)),
                                   p_995 = apply(grad,2, function(x) quantile(x,0.995)),
                                   known = gradient_formulas[[i]](lats)
                                   
    )
  }
  names(t_list[[i]]) <- c(5,10,20,50)
}
names(t_list) <- c("linear","flat","quadratic","steep")

# save in long format
t_long <- bind_rows(lapply(t_list, function(x) bind_rows(x,.id = "sample_size")),.id="gradient_type")
t_long$sample_size <- as.numeric(t_long$sample_size)

saveRDS(t_long,"./results/simulation/simulated_gradients_temperature_estimates.rds")
