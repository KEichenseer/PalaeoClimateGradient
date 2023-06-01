### goodness of fit of the modelled, simulated gradients

library(dplyr)

source("R/functions/model_components/gradient.R")
source("R/functions/auxiliary_functions.R")

# read simulated gradient models
four_gradients_model <- readRDS("D://OneDrive - Durham University/projects/four_gradients_model_proc.rds")

lats <- seq(0,90,1)
grad <- gradient(lats, four_gradients_model[[1]][[1]])

gradient_formulas <- rev(list( function(x) 30-x^2/1080, 
                               function(x) 30-x^2/360, 
                               function(x)  11.25*atan(-0.1*x+3.5)+15.25, 
                               function(x) 8.75*atan(-0.5*x+15)+12))

plot(lats,gradient_formulas[[1]](lats), type = "l")
points(lats, grad[30000,], col = "orange")


# basic R squared

R_squared <- list(NULL)

for(g in 1:4) {
  R_squared[[g]] <- list(NULL)
  for(s in 1:4) {
    grad <- gradient(lats, four_gradients_model[[g]][[s]])
    true_grad <- gradient_formulas[[g]](lats)
    SSres <- apply(grad,1,function(x) sum((x-true_grad)^2))
    SStot <- sum((true_grad-mean(true_grad))^2)
    R_squared[[g]][[s]] <- 1-SSres/SStot
  }
}

### Try weighted R squared (by latitudinal area)

lats <- seq(0.5,89.5,1)
alpha1 <- seq(1,90,1)
alpha2 <- seq(0,89,1)
latweight <- sin(pi*alpha1/180) - sin(pi*alpha2/180) # sums to 1 - careful when using lat. subsets

R_squared_weighted <- list(NULL)

for(g in 1:4) {
  R_squared_weighted[[g]] <- list(NULL)
  for(s in 1:4) {
    grad <- gradient(lats, four_gradients_model[[g]][[s]])
    true_grad <- gradient_formulas[[g]](lats)
    SSres <- apply(grad,1,function(x) sum(latweight*(x-true_grad)^2))
    SStot <- sum(latweight*(true_grad-mean(true_grad))^2)
    R_squared_weighted[[g]][[s]] <- 1-SSres/SStot
    print(c(g,s))
  }
}

# Use Gelman's R squared for bayesian regression models
# http://www.stat.columbia.edu/~gelman/research/unpublished/bayes_R2.pdf
lats <- seq(0,90,1)

R_squared_gelman <- list(NULL)
for(g in 1:4) {
  R_squared_gelman[[g]] <- list(NULL)
  for(s in 1:4) {
    
    grad <- gradient(lats, four_gradients_model[[g]][[s]])
    if(s == 4)     grad <- gradient(lats, four_gradients_model[[g]][[s+1]])

    true_grad <- gradient_formulas[[g]](lats)
    
    var_fit <- apply(grad,1,function(x) var(x))
    var_res <- apply(grad,1,function(x) var(x-true_grad))
    
    #SSres <- apply(grad,1,function(x) sum(latweight*(x-true_grad)^2))
    #SStot <- sum(latweight*(true_grad-mean(true_grad))^2)
    R_squared_gelman[[g]][[s]] <- var_fit/(var_fit+var_res)
    print(c(g,s))
  }
}

# 
# 
# hist(R_squared,1000, xlim = c(0,1))
# quantile(R_squared[[3]][[1]],probs = c(0.025,0.5,0.975))
# quantile(R_squared_weighted[[3]][[1]],probs = c(0.025,0.5,0.975))
# 
# plot(true_grad, ylim = c(0,45))
# its <-  sample(1:200000,50)
# for(i in its) points(grad[i,], col = "blue", type = "l")

# Compare standard R squared and Gelman's R squared
par(mfrow=c(4,4), mar = c(4,4,.5,.5))
for(g in 1:4) {
  for(s in 1:4) {
    plot(1,quantile(R_squared[[g]][[s]],probs = c(0.5)), ylim = c(0,1.05), cex = 2, lwd = 2, ylab = "R squared", xlab = NA)   
    points(c(1,1),quantile(R_squared[[g]][[s]],probs = c(0.025,0.975)),type = "l")
    points(c(1,1),quantile(R_squared[[g]][[s]],probs = c(0.05,0.95)),type = "l", lwd = 2)
    points(c(1,1),quantile(R_squared[[g]][[s]],probs = c(0.1,0.9)),type = "l", lwd = 3)
    points(1.1,quantile(R_squared_gelman[[g]][[s]],probs = c(0.5)), col = "red", pch = 5, cex = 1.5, lwd = 2)   
    points(c(1.1,1.1),quantile(R_squared_gelman[[g]][[s]],probs = c(0.025,0.975)),type = "l", col = "red")
    points(c(1.1,1.1),quantile(R_squared_gelman[[g]][[s]],probs = c(0.05,0.95)),type = "l", lwd = 2, col = "red")
    points(c(1.1,1.1),quantile(R_squared_gelman[[g]][[s]],probs = c(0.1,0.9)),type = "l", lwd = 3, col = "red")
    
  }
}

# Gelman's R squared is bound to [0,1], which is desirable --> we use this one.


# organise list for figure

ssize <- c(5,10,20,"eocene")

R_squared_list <- list(NULL)
for(g in 1:4) {
  R_squared_list[[g]] <- list(NULL)
  
  for(s in 1:4) {
    R_squared_list[[g]][[s]] <- list(mean = mean(R_squared_gelman[[g]][[s]]),
                                     median = median(R_squared_gelman[[g]][[s]]),
                                     p_005 = quantile(R_squared_gelman[[g]][[s]],0.005),
                                     p_025 = quantile(R_squared_gelman[[g]][[s]],0.025),
                                     p_05 = quantile(R_squared_gelman[[g]][[s]],0.05),
                                     p_100 = quantile(R_squared_gelman[[g]][[s]],0.1),
                                     p_900 = quantile(R_squared_gelman[[g]][[s]],0.9),
                                     p_95 = quantile(R_squared_gelman[[g]][[s]],0.95),
                                     p_975 = quantile(R_squared_gelman[[g]][[s]],0.975),
                                     p_995 = quantile(R_squared_gelman[[g]][[s]],0.995))
  }
  names(R_squared_list[[g]]) <- ssize
  
}
names(R_squared_list) <- c("extreme_icehouse","icehouse","greenhouse","extreme_greenhouse")



# save in long format
r2_long <- bind_rows(lapply(R_squared_list, function(x) bind_rows(x,.id = "sample_size")),.id="gradient_type")
saveRDS(r2_long,"./results/simulation/simulated_gradients_r_squared_gelman_unweighted.rds")

# set par mfrow back to 1,1
par(mfrow = c(1,1))
