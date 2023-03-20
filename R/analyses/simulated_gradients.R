### Supplementary figure: different gradient shapes
source("R/functions/model_components/climate_model_modern.R")
source("R/functions/model_components/gradient.R")
source("R/functions/model_processing/combine_posterior.R")
source("R/functions/auxiliary_functions.R")

# source options for analysis
source("R/options.R")
#

ylim = c(0,30)
xlim = c(0,90)
lat <- abs(c(-59.42, -59.10, -54.13, -47.84, -47.84, -47.47, -46.14, -45.80, -45.80, -24.25,  -8.57,  -2.05,  -1.86,   0.49,   4.29,  14.67,
             27.44,  31.00,  32.50,  32.80,  32.90,  35.54,  36.22,  39.02,  40.39,  40.63,  40.89,  43.79,  44.11,  50.04,  50.39,  58.45,
             76.45,  77.13))

n_rep <- 100
sample_sizes <- c(5,10,20,50)

## 1. Linear
m_lin <- list(NULL)
for(ss in 1:length(sample_sizes)) {
  m_lin[[ss]] <- list()
  print(paste("sample_size =",sample_sizes[ss]))
  set.seed(1)
  for(rep in 1:n_rep) {
    lat <- runif(sample_sizes[ss],0,90)
    temp <- 30-lat/3
    m_lin[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                          n_thin = 10,
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
  temp_lin[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.025,0.975)))
}

par(mfrow=c(4,4))
for(i in 1:4) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,30), xlab = "abs. latitude", ylab = "temperature")
  points(lats, 30-lats/3, type = "l", lwd = 2)
  #plot_gradient(m_lin_proc[[i]], add = T)
  error_polygon(lats,temp_lin[[i]][2,],temp_lin[[i]][3,], col = rgb(1,0,0,0.2))
  points(lats,temp_lin[[i]][1,], type = "l", lwd = 2, col = rgb(1,0,0,0.75))
}


## 2. Flat
m_flat <- list(NULL)
for(ss in 1:length(sample_sizes)) {
  m_flat[[ss]] <- list()
  print(paste("sample_size =",sample_sizes[ss]))
  set.seed(1)
  for(rep in 1:n_rep) {
    lat <- runif(sample_sizes[ss],0,90)
    temp <- 30-lat/3
    m_flat[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                                          n_thin = 10,
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
  temp_flat[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.025,0.975)))
}

for(i in 1:4) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,30), xlab = "abs. latitude", ylab = "temperature")
  points(lats, 30-lats/3, type = "l", lwd = 2)
  #plot_gradient(m_flat_proc[[i]], add = T)
  error_polygon(lats,temp_flat[[i]][2,],temp_flat[[i]][3,], col = rgb(1,0,0,0.2))
  points(lats,temp_flat[[i]][1,], type = "l", lwd = 2, col = rgb(1,0,0,0.75))
}


## 3. Quadratic
m_quad <- list(NULL)
for(ss in 1:length(sample_sizes)) {
  m_quad[[ss]] <- list()
  print(paste("sample_size =",sample_sizes[ss]))
  set.seed(1)
  for(rep in 1:n_rep) {
    lat <- runif(sample_sizes[ss],0,90)
    temp <- 30-lat/3
    m_quad[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                                           n_thin = 10,
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
  temp_quad[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.025,0.975)))
}

for(i in 1:4) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,30), xlab = "abs. latitude", ylab = "temperature")
  points(lats, 30-lats/3, type = "l", lwd = 2)
  #plot_gradient(m_quad_proc[[i]], add = T)
  error_polygon(lats,temp_quad[[i]][2,],temp_quad[[i]][3,], col = rgb(1,0,0,0.2))
  points(lats,temp_quad[[i]][1,], type = "l", lwd = 2, col = rgb(1,0,0,0.75))
}


## 4. Steep
m_steep <- list(NULL)
for(ss in 1:length(sample_sizes)) {
  m_steep[[ss]] <- list()
  print(paste("sample_size =",sample_sizes[ss]))
  set.seed(1)
  for(rep in 1:n_rep) {
    lat <- runif(sample_sizes[ss],0,90)
    temp <- 8.5*atan(-0.5*lat+15)+12.5
    m_steep[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                                           n_thin = 10,
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
  temp_steep[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.025,0.975)))
}

for(i in 1:4) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,40), xlab = "abs. latitude", ylab = "temperature")
  points(lats, 8.5*atan(-0.5*lats+15)+12.5, type = "l", lwd = 2, col = "black")
  #plot_gradient(m_steep_proc[[i]], add = T)
  error_polygon(lats,temp_steep[[i]][2,],temp_steep[[i]][3,], col = rgb(1,0,0,0.2))
  points(lats,temp_steep[[i]][1,], type = "l", lwd = 2, col = rgb(1,0,0,0.75))
}



for(i in 1:4) { plot(lats,tempq[[i]][1,], type = "l", lwd = 2)
  error_polygon(lats,tempq[[i]][2,],tempq[[i]][3,])
}

m1 <- run_MCMC_simple(n_iter = 20000,
                      n_thin = 10,
                      x = lat,
                      y = temp1,
                      prior_input = priors,
                      adapt_sd = NULL)
r1 <- combine_posterior(list(m1))
plot_chains(list(m1))

yt <- 8.5*atan(-0.5*lats+15)+12.5
plot(lats,yt,type = "l")
abline(h=0)
abline(h=25)
