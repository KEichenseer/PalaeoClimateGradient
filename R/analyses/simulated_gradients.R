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

gradient_formulas <- rev(list( function(x) 30-x^2/1080, function(x) 30-x^2/405, 
                               function(x)  11.25*atan(-0.1*x+3.5)+15.25, function(x) 8.5*atan(-0.5*x+15)+12.5))


### Model 4 different gradients
##
## 1. extreme icehouse
grad_type <- 1
m_1 <- list(NULL)
set.seed(1)
for(ss in 1:length(sample_sizes)) {
  m_1[[ss]] <- list()
  print(paste("sample_size =",sample_sizes[ss]))
  for(rep in 1:n_rep) {
    lat <- 180/pi * asin(runif(sample_sizes[ss],0,1))
    temp <- gradient_formulas[[grad_type]](lat)
    m_1[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                          n_thin = 50,
                          x = lat,
                          y = temp,
                          prior_input = priors,
                          adapt_sd = NULL)
    print(paste("repetition =",rep))
    
  }
  
  
}

m_1_proc <- list()
for(i in 1:length(sample_sizes)) m_1_proc[[i]] <- combine_posterior(m_1[[i]],burnin = 2500)

lats <- seq(0,90,0.5)
temp_1 <- list()
lats <- seq(0,90,0.5)
for(i in 1:4) {
  temp <- gradient(lats,m_1_proc[[i]])
  temp_1[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

par(mfrow=c(4,4))
for(i in 1:4) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,40), xlab = "abs. latitude", ylab = "temperature")
  points(lats, gradient_formulas[[grad_type]](lats), type = "l", lwd = 2)
  #plot_gradient(m_1_proc[[i]], add = T)
  error_polygon(lats,temp_1[[i]][2,],temp_1[[i]][3,], col = rgb(1,0,0,0.10))
  error_polygon(lats,temp_1[[i]][4,],temp_1[[i]][5,], col = rgb(1,0,0,0.15))
  error_polygon(lats,temp_1[[i]][6,],temp_1[[i]][7,], col = rgb(1,0,0,0.20))
  
  points(lats,temp_1[[i]][1,], type = "l", lwd = 2, col = rgb(1,0,0,0.75))
}


##
## 2. icehouse
grad_type <- 2
m_2 <- list(NULL)
set.seed(1)
for(ss in 1:length(sample_sizes)) {
  m_2[[ss]] <- list()
  print(paste("sample_size =",sample_sizes[ss]))
  for(rep in 1:n_rep) {
    lat <- 180/pi * asin(runif(sample_sizes[ss],0,1))
    temp <- gradient_formulas[[grad_type]](lat)
    m_2[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                                        n_thin = 50,
                                        x = lat,
                                        y = temp,
                                        prior_input = priors,
                                        adapt_sd = NULL)
    print(paste("repetition =",rep))
    
  }
  
  
}

m_2_proc <- list()
for(i in 1:length(sample_sizes)) m_2_proc[[i]] <- combine_posterior(m_2[[i]],burnin = 2500)

lats <- seq(0,90,0.5)
temp_2 <- list()
lats <- seq(0,90,0.5)
for(i in 1:4) {
  temp <- gradient(lats,m_2_proc[[i]])
  temp_2[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

for(i in 1:4) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,40), xlab = "abs. latitude", ylab = "temperature")
  points(lats, gradient_formulas[[grad_type]](lats), type = "l", lwd = 2)
  #plot_gradient(m_2_proc[[i]], add = T)
  error_polygon(lats,temp_2[[i]][2,],temp_2[[i]][3,], col = rgb(1,0,0,0.10))
  error_polygon(lats,temp_2[[i]][4,],temp_2[[i]][5,], col = rgb(1,0,0,0.15))
  error_polygon(lats,temp_2[[i]][6,],temp_2[[i]][7,], col = rgb(1,0,0,0.20))
  
  points(lats,temp_2[[i]][1,], type = "l", lwd = 2, col = rgb(1,0,0,0.75))
}



##
## 3. greenhouse
##
grad_type <- 3
m_3 <- list(NULL)
set.seed(1)
for(ss in 1:length(sample_sizes)) {
  m_3[[ss]] <- list()
  print(paste("sample_size =",sample_sizes[ss]))
  for(rep in 1:n_rep) {
    lat <- 180/pi * asin(runif(sample_sizes[ss],0,1))
    temp <- gradient_formulas[[grad_type]](lat)
    m_3[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                                        n_thin = 50,
                                        x = lat,
                                        y = temp,
                                        prior_input = priors,
                                        adapt_sd = NULL)
    print(paste("repetition =",rep))
    
  }
  
  
}

m_3_proc <- list()
for(i in 1:length(sample_sizes)) m_3_proc[[i]] <- combine_posterior(m_3[[i]],burnin = 2500)

lats <- seq(0,90,0.5)
temp_3 <- list()
lats <- seq(0,90,0.5)
for(i in 1:4) {
  temp <- gradient(lats,m_3_proc[[i]])
  temp_3[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

for(i in 1:4) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,40), xlab = "abs. latitude", ylab = "temperature")
  points(lats, gradient_formulas[[grad_type]](lats), type = "l", lwd = 2)
  #plot_gradient(m_3_proc[[i]], add = T)
  error_polygon(lats,temp_3[[i]][2,],temp_3[[i]][3,], col = rgb(1,0,0,0.10))
  error_polygon(lats,temp_3[[i]][4,],temp_3[[i]][5,], col = rgb(1,0,0,0.15))
  error_polygon(lats,temp_3[[i]][6,],temp_3[[i]][7,], col = rgb(1,0,0,0.20))
  
  points(lats,temp_3[[i]][1,], type = "l", lwd = 2, col = rgb(1,0,0,0.75))
}


##
## 4. extreme greenhouse
##
grad_type <- 4
m_4 <- list(NULL)
set.seed(1)
for(ss in 1:length(sample_sizes)) {
  m_4[[ss]] <- list()
  print(paste("sample_size =",sample_sizes[ss]))
  for(rep in 1:n_rep) {
    lat <- 180/pi * asin(runif(sample_sizes[ss],0,1))
    temp <- gradient_formulas[[grad_type]](lat)
    m_4[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                                        n_thin = 50,
                                        x = lat,
                                        y = temp,
                                        prior_input = priors,
                                        adapt_sd = NULL)
    print(paste("repetition =",rep))
    
  }
  
  
}

m_4_proc <- list()
for(i in 1:length(sample_sizes)) m_4_proc[[i]] <- combine_posterior(m_4[[i]],burnin = 2500)

lats <- seq(0,90,0.5)
temp_4 <- list()
lats <- seq(0,90,0.5)
for(i in 1:4) {
  temp <- gradient(lats,m_4_proc[[i]])
  temp_4[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

for(i in 1:4) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,40), xlab = "abs. latitude", ylab = "temperature")
  points(lats, gradient_formulas[[grad_type]](lats), type = "l", lwd = 2)
  #plot_gradient(m_4_proc[[i]], add = T)
  error_polygon(lats,temp_4[[i]][2,],temp_4[[i]][3,], col = rgb(1,0,0,0.10))
  error_polygon(lats,temp_4[[i]][4,],temp_4[[i]][5,], col = rgb(1,0,0,0.15))
  error_polygon(lats,temp_4[[i]][6,],temp_4[[i]][7,], col = rgb(1,0,0,0.20))
  
  points(lats,temp_4[[i]][1,], type = "l", lwd = 2, col = rgb(1,0,0,0.75))
}


m_lin <- four_gradients_model[[1]]
m_flat <- four_gradients_model[[2]]
m_quad <- four_gradients_model[[3]]
m_steep <- four_gradients_model[[4]]

### save list of all temperature outputs
four_gradients_model <- list(m_1_proc,m_2_proc,m_3_proc,m_4_proc)
four_gradients_temp <- list(temp_1,temp_2,temp_3,temp_4)

saveRDS(four_gradients_model, "D://OneDrive - Durham University/projects/four_gradients_model_proc.rds")
saveRDS(four_gradients_temp, "D://OneDrive - Durham University/projects/four_gradients_temp.rds")

#gradient_formulas <- list( function(x) 30-x/3, function(x) rep(25,length(x)), function(x) 30-x^2/270, function(x) 8.5*atan(-0.5*x+15)+12.5)


plot(0,0,type = "n", xlim = c(0,90), ylim = c(-2,33), xlab = "abs. latitude", ylab = "temperature")
cols = rev(c("red", "orange", "turquoise", "blue"))
for(i in 1:4) points(lats,gradient_formulas[[i]](lats),type = "l", col = cols[i], lwd = 2)
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
names(t_list) <- c("extreme_icehouse","icehouse","greenhouse","extreme_greenhouse")

# save in long format
t_long <- bind_rows(lapply(t_list, function(x) bind_rows(x,.id = "sample_size")),.id="gradient_type")
t_long$sample_size <- as.numeric(t_long$sample_size)

saveRDS(t_long,"./results/simulation/simulated_gradients_temperature_estimates.rds")
