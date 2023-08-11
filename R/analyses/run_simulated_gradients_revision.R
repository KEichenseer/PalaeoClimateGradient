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
eocene_lat <- c(-59.42, -59.10, -54.13, -47.84, -47.84, -47.47, -46.14, -45.80, -45.80, -24.25, -8.57, -2.05, -1.86, 0.49, 4.29, 14.67, 27.44, 31.00, 32.50, 32.80, 32.90, 35.54, 36.22, 39.02, 40.39, 40.63, 40.89, 43.79, 44.11, 50.04, 50.39, 58.45, 76.45, 77.13)
# 
generate_lat_samples <- function(n) 180/pi * asin(runif(n,0,1))

# spatially clustered

sd1 <- sd(eocene_lat[which(eocene_lat < -40)])
n1 <- length(eocene_lat[which(eocene_lat < -40)])
p1 <- n1/length(eocene_lat) 

sd2 <- sd(eocene_lat[which(eocene_lat < 20 & eocene_lat > -40)])
n2 <- length(eocene_lat[which(eocene_lat < 20 & eocene_lat > -40)])
p2 <- n2/length(eocene_lat) 

sd3 <- sd(eocene_lat[which(eocene_lat < 70 & eocene_lat > 20)])
n3 <- length(eocene_lat[which(eocene_lat < 70 & eocene_lat > 20)])
p3 <- n3/length(eocene_lat) 

sd4 <- sd(eocene_lat[which(eocene_lat > 70)])
n4 <- length(eocene_lat[which(eocene_lat > 70)])
p4 <- n4/length(eocene_lat) 

sd_vec <- c(sd1,sd2,sd3,sd4)
prob_vec <- c(p1,p2,p3,p4)
generate_lat_samples <- function(n, sd_vec, prob_vec) {
  
  n_centers <- 4
  lat_centers <-    180/pi * asin(runif(n_centers,0,1)) * sample(c(1,-1),n_centers, replace = T)

  center_index <- sample(1:4,size = n, replace = T, prob = prob_vec)
  
  lats <- rnorm(n,mean = lat_centers[center_index], sd = sd_vec[center_index])
  
  lats[which(lats < -90)] <- -180 - lats[which(lats < -90)] 
  lats[which(lats > 90)] <- 180 - lats[which(lats > 90)] 
  # return absolute latitude
  abs(lats)
}

set.seed(1)


# number of samples (repetitions)
n_rep <- 200 # careful, with 1000 reps it takes hours to run the whole script
# sample sizes in the repetitions (4 different)
sample_sizes <- c(5,10,20,50)

gradient_formulas_static <- rev(list( function(x) 30-x^2/1080, # extreme greenhouse
                               function(x) 30-x^2/360, 
                               function(x) gradient(x, c(-1.7748065, 30.1152231, 42.0699179,  0.1057308)), # now use modern gradient instead # 11.25*atan(-0.1*x+3.5)+15.25, 
                               function(x) 8.75*atan(-0.5*x+15)+12)) # extreme icehouse

# with noise
gradient_formulas <- rev(list( function(x) 30-x^2/1080 + rnorm(length(x), 0, 4.39), # extreme greenhouse
                               function(x) 30-x^2/360 + rnorm(length(x), 0, 4.39), 
                               function(x) gradient(x, c(-1.7748065, 30.1152231, 42.0699179,  0.1057308)) + rnorm(length(x), 0, 4.39), # now use modern gradient instead # 11.25*atan(-0.1*x+3.5)+15.25, 
                               function(x) 8.75*atan(-0.5*x+15)+12 + rnorm(length(x), 0, 4.39))) # extreme icehouse


# read and process temperature and latitude data
modern_samples <- readRDS("./results/modern/modern_sample.RDS")
eocene_lat <- abs(modern_samples[[1]][,"p_lat"])
modern_1deg <- readRDS("results/modern/empirical_median_1-deg-lat.rds")
modern_temp_eocene_lat <- modern_1deg[sapply(floor(eocene_lat)+0.5,function(x) which(x==modern_1deg[,"lat"])),2] # sapply(1:34, function(y) mean(sapply(1:100, function(x) modern_samples[[x]][y,2])))
# interpolate so modern temperature is available at same resolution for plotting
linear_interp <- linear_interp <- function(a, b) {
  # Sort the matrix by the x column
  b <- b[order(b[,1]),]
  
  # Determine the minimum and maximum x values in b
  min_x <- min(b[,1])
  max_x <- max(b[,1])
  
  # Initialize a vector to store the interpolated y values
  y_interp <- numeric(length(a))
  
  # Interpolate the y values for each element of a
  for (i in seq_along(a)) {
    # If the x value is outside the range of b, use the first or last y value
    if (a[i] < min_x) {
      y_interp[i] <- b[1,2]
    } else if (a[i] > max_x) {
      y_interp[i] <- b[nrow(b),2]
    } else {
      # Find the indices of the two closest x values in b
      idx <- findInterval(a[i], b[,1], left.open = TRUE, all.inside = T)
      
      # Perform linear interpolation between the two closest points
      x1 <- b[idx,1]
      x2 <- b[idx+1,1]
      y1 <- b[idx,2]
      y2 <- b[idx+1,2]
      y_interp[i] <- y1 + (y2 - y1) * (a[i] - x1) / (x2 - x1)
    }
  }
  
  return(y_interp)
}
modern_1deg_interpolated <- linear_interp(a = seq(0,90,0.5),b = modern_1deg)


eocene_sample_nIter <- 30000 # 250000
random_sample_nIter <- 2000 # 10000

burnin_random_sample <- 1000 # 5000
burnin_eocene_sample <- 10000 # 50000

### Model 4 different gradients
##
## 1. extreme icehouse
grad_type <- 1
m_1 <- list(NULL)
set.seed(1)
for(ss in 1:length(sample_sizes)) {
  m_1[[ss]] <- list()
  for(rep in 1:n_rep) {
    lat <- generate_lat_samples(sample_sizes[ss], sd_vec = sd_vec, prob_vec = prob_vec) 
    temp <- gradient_formulas[[grad_type]](lat)
    m_1[[ss]][[rep]] <- run_MCMC_simple(n_iter = random_sample_nIter,
                          n_thin = 25,
                          x = lat,
                          y = temp,
                          prior_input = priors,
                          adapt_sd = NULL)
    print(paste("gradient =",grad_type,"sample_size =",sample_sizes[ss],"repetition =",rep))
    
  }
  
  
}

# slot 5 is the eocene distribution
m_1[[5]] <- list() 

m_1[[5]][[1]] <-  run_MCMC_simple(n_iter = eocene_sample_nIter,
                                  n_thin = 20,
                                  x = eocene_lat,
                                  y = gradient_formulas[[grad_type]](eocene_lat),
                                  prior_input = priors,
                                  adapt_sd = NULL)

m_1_proc <- list()
for(i in 1:length(sample_sizes)) m_1_proc[[i]] <- combine_posterior(m_1[[i]],burnin = burnin_random_sample)
m_1_proc[[5]] <- combine_posterior(m_1[[5]],burnin = burnin_eocene_sample)
dim(m_1_proc[[5]])

lats <- seq(0,90,0.5)
temp_1 <- list()
lats <- seq(0,90,0.5)
for(i in 1:5) {
  temp <- gradient(lats,m_1_proc[[i]])
  temp_1[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

par(mfrow=c(5,5), mar = c(4,4,.5,.5))
for(i in 1:5) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,40), xlab = "abs. latitude", ylab = "temperature")
  points(lats, gradient_formulas[[grad_type]](lats), type = "l", lwd = 2)
  #plot_gradient(m_1_proc[[i]], add = T)
  error_polygon(lats,temp_1[[i]][2,],temp_1[[i]][3,], col = rgb(1,0,0,0.10))
  error_polygon(lats,temp_1[[i]][4,],temp_1[[i]][5,], col = rgb(1,0,0,0.15))
  error_polygon(lats,temp_1[[i]][6,],temp_1[[i]][7,], col = rgb(1,0,0,0.20))
  
  points(lats,temp_1[[i]][1,], type = "p", cex = .5, col = rgb(1,0,0,0.75))
}


##
## 2. icehouse
grad_type <- 2
m_2 <- list(NULL)
set.seed(1)
for(ss in 1:length(sample_sizes)) {
  m_2[[ss]] <- list()
  for(rep in 1:n_rep) {
    lat <- generate_lat_samples(sample_sizes[ss], sd_vec = sd_vec, prob_vec = prob_vec) 
    temp <- gradient_formulas[[grad_type]](lat)
    m_2[[ss]][[rep]] <- run_MCMC_simple(n_iter = random_sample_nIter,
                                        n_thin = 25,
                                        x = lat,
                                        y = temp,
                                        prior_input = priors,
                                        adapt_sd = NULL)
    print(paste("gradient =",grad_type,"sample_size =",sample_sizes[ss],"repetition =",rep))
    
  }
  
  
}


m_2[[5]] <- list() 

m_2[[5]][[1]] <-  run_MCMC_simple(n_iter = eocene_sample_nIter,
                                  n_thin = 20,
                                  x = eocene_lat,
                                  y = gradient_formulas[[grad_type]](eocene_lat),
                                  prior_input = priors,
                                  adapt_sd = NULL)

m_2_proc <- list()
for(i in 1:length(sample_sizes)) m_2_proc[[i]] <- combine_posterior(m_2[[i]],burnin = burnin_random_sample)
m_2_proc[[5]] <- combine_posterior(m_2[[5]],burnin = burnin_eocene_sample)
dim(m_2_proc[[5]])

lats <- seq(0,90,0.5)
temp_2 <- list()
lats <- seq(0,90,0.5)
for(i in 1:5) {
  temp <- gradient(lats,m_2_proc[[i]])
  temp_2[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

for(i in 1:5) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,40), xlab = "abs. latitude", ylab = "temperature")
  points(lats, gradient_formulas[[grad_type]](lats), type = "l", lwd = 2)
  #plot_gradient(m_2_proc[[i]], add = T)
  error_polygon(lats,temp_2[[i]][2,],temp_2[[i]][3,], col = rgb(1,0,0,0.10))
  error_polygon(lats,temp_2[[i]][4,],temp_2[[i]][5,], col = rgb(1,0,0,0.15))
  error_polygon(lats,temp_2[[i]][6,],temp_2[[i]][7,], col = rgb(1,0,0,0.20))
  
  points(lats,temp_2[[i]][1,], type = "p", cex = .5, col = rgb(1,0,0,0.75))
}



##
## 3. greenhouse
##
grad_type <- 3
m_3 <- list(NULL)
set.seed(1)
for(ss in 1:length(sample_sizes)) {
  m_3[[ss]] <- list()
  for(rep in 1:n_rep) {
    lat <- generate_lat_samples(sample_sizes[ss], sd_vec = sd_vec, prob_vec = prob_vec) 
    temp <- gradient_formulas[[grad_type]](lat)
    m_3[[ss]][[rep]] <- run_MCMC_simple(n_iter = random_sample_nIter,
                                        n_thin = 25,
                                        x = lat,
                                        y = temp,
                                        prior_input = priors,
                                        adapt_sd = NULL)
    print(paste("gradient =",grad_type,"sample_size =",sample_sizes[ss],"repetition =",rep))
    
  }
  
  
}


m_3[[5]] <- list() 

m_3[[5]][[1]] <-  run_MCMC_simple(n_iter = eocene_sample_nIter,
                                  n_thin = 20,
                                  x = eocene_lat,
                                  y = gradient_formulas[[grad_type]](eocene_lat),
                                  prior_input = priors,
                                  adapt_sd = NULL)

m_3_proc <- list()
for(i in 1:length(sample_sizes)) m_3_proc[[i]] <- combine_posterior(m_3[[i]],burnin = burnin_random_sample)
m_3_proc[[5]] <- combine_posterior(m_3[[5]],burnin = burnin_eocene_sample)
dim(m_3_proc[[5]])

lats <- seq(0,90,0.5)
temp_3 <- list()
lats <- seq(0,90,0.5)
for(i in 1:5) {
  temp <- gradient(lats,m_3_proc[[i]])
  temp_3[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

for(i in 1:5) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,40), xlab = "abs. latitude", ylab = "temperature")
  points(lats, gradient_formulas[[grad_type]](lats), type = "l", lwd = 2)
  #plot_gradient(m_3_proc[[i]], add = T)
  error_polygon(lats,temp_3[[i]][2,],temp_3[[i]][3,], col = rgb(1,0,0,0.10))
  error_polygon(lats,temp_3[[i]][4,],temp_3[[i]][5,], col = rgb(1,0,0,0.15))
  error_polygon(lats,temp_3[[i]][6,],temp_3[[i]][7,], col = rgb(1,0,0,0.20))
  
  points(lats,temp_3[[i]][1,], type = "p", cex = .5, col = rgb(1,0,0,0.75))
}


##
## 4. extreme greenhouse
##
grad_type <- 4
m_4 <- list(NULL)
set.seed(1)
for(ss in 1:length(sample_sizes)) {
  m_4[[ss]] <- list()
  for(rep in 1:n_rep) {
    lat <- generate_lat_samples(sample_sizes[ss], sd_vec = sd_vec, prob_vec = prob_vec) 
    temp <- gradient_formulas[[grad_type]](lat)
    m_4[[ss]][[rep]] <- run_MCMC_simple(n_iter = random_sample_nIter,
                                        n_thin = 25,
                                        x = lat,
                                        y = temp,
                                        prior_input = priors,
                                        adapt_sd = NULL)
    print(paste("gradient =",grad_type,"sample_size =",sample_sizes[ss],"repetition =",rep))
    
  }
  
  
}



m_4[[5]] <- list() 

m_4[[5]][[1]] <-  run_MCMC_simple(n_iter = eocene_sample_nIter,
                                  n_thin = 20,
                                  x = eocene_lat,
                                  y = gradient_formulas[[grad_type]](eocene_lat),
                                  prior_input = priors,
                                  adapt_sd = NULL)

m_4_proc <- list()
for(i in 1:length(sample_sizes)) m_4_proc[[i]] <- combine_posterior(m_4[[i]],burnin = burnin_random_sample)
m_4_proc[[5]] <- combine_posterior(m_4[[5]],burnin = burnin_eocene_sample)
dim(m_4_proc[[5]])

lats <- seq(0,90,0.5)
temp_4 <- list()
lats <- seq(0,90,0.5)
for(i in 1:5) {
  temp <- gradient(lats,m_4_proc[[i]])
  temp_4[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

for(i in 1:5) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,40), xlab = "abs. latitude", ylab = "temperature")
  points(lats, gradient_formulas[[grad_type]](lats), type = "l", lwd = 2)
  #plot_gradient(m_4_proc[[i]], add = T)
  error_polygon(lats,temp_4[[i]][2,],temp_4[[i]][3,], col = rgb(1,0,0,0.10))
  error_polygon(lats,temp_4[[i]][4,],temp_4[[i]][5,], col = rgb(1,0,0,0.15))
  error_polygon(lats,temp_4[[i]][6,],temp_4[[i]][7,], col = rgb(1,0,0,0.20))
  
  points(lats,temp_4[[i]][1,], type = "p", cex = .5, col = rgb(1,0,0,0.75))
}



##
## 5. modern
##
grad_type <- 5
m_5 <- list(NULL)
set.seed(1)
for(ss in 1:length(sample_sizes)) {
  m_5[[ss]] <- list()
  for(rep in 1:n_rep) {
    lat <- floor(generate_lat_samples(sample_sizes[ss], sd_vec = sd_vec, prob_vec = prob_vec))+0.5
    temp <-  modern_1deg[sapply(lat,function(x) which(x==modern_1deg[,"lat"])),2]
    m_5[[ss]][[rep]] <- run_MCMC_simple(n_iter = 10000,
                                        n_thin = 25,
                                        x = lat,
                                        y = temp,
                                        prior_input = priors,
                                        adapt_sd = NULL)
    print(paste("gradient =",grad_type,"sample_size =",sample_sizes[ss],"repetition =",rep))
    
  }
  
  
}



m_5[[5]] <- list() 

m_5[[5]][[1]] <-  run_MCMC_simple(n_iter = eocene_sample_nIter,
                                  n_thin = 20,
                                  x = eocene_lat,
                                  y = modern_1deg[sapply(floor(eocene_lat)+0.5,function(x) which(x==modern_1deg[,"lat"])),2],
                                  prior_input = priors,
                                  adapt_sd = NULL)

m_5_proc <- list()
for(i in 1:length(sample_sizes)) m_5_proc[[i]] <- combine_posterior(m_5[[i]],burnin = burnin_random_sample)
m_5_proc[[5]] <- combine_posterior(m_5[[5]],burnin = burnin_eocene_sample)
dim(m_5_proc[[5]])

lats <- seq(0,90,0.5)
temp_5 <- list()
lats <- seq(0,90,0.5)
for(i in 1:5) {
  temp <- gradient(lats,m_5_proc[[i]])
  temp_5[[i]] <- apply(temp,2,function(x) quantile(x,probs = c(0.5,0.005,0.995,0.025,0.975,0.1,0.9)))
}

for(i in 1:5) {
  plot(0,0,type = "n", xlim = c(0,90), ylim = c(0,40), xlab = "abs. latitude", ylab = "temperature")
  points(modern_1deg$lat,modern_1deg$SST, type = "l", lwd = 2)
  #plot_gradient(m_5_proc[[i]], add = T)
  error_polygon(lats,temp_5[[i]][2,],temp_5[[i]][3,], col = rgb(1,0,0,0.10))
  error_polygon(lats,temp_5[[i]][4,],temp_5[[i]][5,], col = rgb(1,0,0,0.15))
  error_polygon(lats,temp_5[[i]][6,],temp_5[[i]][7,], col = rgb(1,0,0,0.20))
  
  points(lats,temp_5[[i]][1,], type = "p", cex = .5, col = rgb(1,0,0,0.75))
}



m_lin <- four_gradients_model[[1]]
m_flat <- four_gradients_model[[2]]
m_quad <- four_gradients_model[[3]]
m_steep <- four_gradients_model[[4]]

### save list of all temperature outputs
four_gradients_model <- list(m_1_proc,m_2_proc,m_3_proc,m_4_proc,m_5_proc)
four_gradients_temp <- list(temp_1,temp_2,temp_3,temp_4,temp_5)

#saveRDS(four_gradients_model, "D://OneDrive - Durham University/projects/four_gradients_model_proc.rds")
#saveRDS(four_gradients_temp, "D://OneDrive - Durham University/projects/four_gradients_temp.rds")

#gradient_formulas <- list( function(x) 30-x/3, function(x) rep(25,length(x)), function(x) 30-x^2/270, function(x) 8.5*atan(-0.5*x+15)+12.5)

gradient_formulas_static2 <- (gradient_formulas_static)
gradient_formulas_static2[[5]] <- function(x) modern_1deg_interpolated

# gradient_formulas <- rev(list(
#   function(x) modern_1deg_interpolated,
#   function(x) 30-x^2/1080,
#   function(x) 30-x^2/360,
#   function(x)  11.25*atan(-0.1*x+3.5)+15.25,
#   function(x) 8.75*atan(-0.5*x+15)+12))



plot(0,0,type = "n", xlim = c(0,90), ylim = c(-2,33), xlab = "abs. latitude", ylab = "temperature")
cols = rev(c("green","red", "orange", "white", "blue"))
for(i in 1:5) points(lats,gradient_formulas_static2[[i]](lats),type = "l", col = cols[i], lwd = 2)
### get temperatures from combining all gradients
t_list <- list()
for(i in 1:5) { # this can take a few minutes
  t_list[[i]] <- list()
  for(j in 1:5) {

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
  names(t_list[[i]]) <- c(5,10,20,50,"eocene")
}
names(t_list) <- c("extreme_icehouse","icehouse","greenhouse","extreme_greenhouse","modern")

# save in long format
t_long <- bind_rows(lapply(t_list, function(x) bind_rows(x,.id = "sample_size")),.id="gradient_type")
#t_long$sample_size <- as.numeric(t_long$sample_size) # no longer works with `eocene`

saveRDS(t_long,"./results/simulation/simulated_gradients_temperature_estimates_with_noise_and_geo_clustering.rds")
