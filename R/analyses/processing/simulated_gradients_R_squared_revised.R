### goodness of fit of the modelled, simulated gradients

library(dplyr)

source("R/functions/model_components/gradient.R")
source("R/functions/auxiliary_functions.R")

# read simulated gradient models
four_gradients_model <- readRDS("D://OneDrive - Durham University/projects/four_gradients_revised_model_proc_M_wider_noise.rds")

lats <- seq(0,90,1)
grad <- gradient(lats, four_gradients_model[[1]][[1]])

gradient_formulas_static <- rev(list( function(x) 30-x^2/1080, # extreme greenhouse
                                      function(x) 30-x^2/360, 
                                      function(x) {
                                        x[which(x>=90)] <- 89.5
                                        modern_1deg[sapply(floor(x)+0.5,function(y) which(y==modern_1deg[,"lat"])),2]
                                      }, #gradient(x, c(-1.775, 30.115, 42.070,  0.106)), # now use modern gradient instead # 11.25*atan(-0.1*x+3.5)+15.25, 
                                      function(x) 8.75*atan(-0.5*x+15)+12)) # extreme icehouse

# Use Gelman's R squared for bayesian regression models
# http://www.stat.columbia.edu/~gelman/research/unpublished/bayes_R2.pdf
lats <- seq(0,90,1)

R_squared_gelman <- list(NULL)
for(g in 1:4) {
  R_squared_gelman[[g]] <- list(NULL)
  for(s in 1:4) {
    
    grad <- gradient(lats, four_gradients_model[[g]][[s]])
    #if(s == 4)     grad <- gradient(lats, four_gradients_model[[g]][[s+1]])

    true_grad <- gradient_formulas_static[[g]](lats)
    
    var_fit <- apply(grad,1,function(x) var(x))
    var_res <- apply(grad,1,function(x) var(x-true_grad))
    
    #SSres <- apply(grad,1,function(x) sum(latweight*(x-true_grad)^2))
    #SStot <- sum(latweight*(true_grad-mean(true_grad))^2)
    R_squared_gelman[[g]][[s]] <- var_fit/(var_fit+var_res)
    print(c(g,s))
  }
}

#

# organise list for figure

ssize <- c(5,10,20,34)

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
saveRDS(r2_long,"./results/simulation/simulated_gradients_r_squared_gelman_unweighted_revised_M_wider_noise.rds")

# set par mfrow back to 1,1
par(mfrow = c(1,1))

i2 <- ceiling(which(R_squared_gelman[[4]][[4]] < 0.1) / 300)
i1 <- which(R_squared_gelman[[4]][[4]] < 0.1) %% 200 + 1
par1 <- (vapply(1:length(i1), function(x) unlist(m_4[[4]][[i2[x]]]$params[i1[x],]), numeric(6)))
grad <- gradient(lats,t(par1))
plot(lats,gradient_formulas_static[[4]](lats), type = "l", lwd = 3)
points(lats, grad[1,], type = "l")
R_squared_gelman[[4]][[4]][which(R_squared_gelman[[4]][[4]] < 0.1)][1]
