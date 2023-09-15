# Function to generate median temperatures and confidence intervals based on 
# a gradient from the model output

temp_from_gradient <- function(lat,model_out, burnin = NULL) {
  
  # load gradient function
  source("R/functions/model_components/gradient.R")
  
  if (is.null(burnin)) first_it <- 1 else first_it = burnin
  last_it <- nrow(model_out)
  
  grad <- gradient(lat,model_out[first_it:last_it,c("A","dKA","M","B")],sdy=0)
  
  out <- data.frame(lat = lat, median = apply(grad,2,median), 
                    l_ci_95 = apply(grad,2, function(x) quantile(x,0.025)), 
                    u_ci_95 = apply(grad,2, function(x) quantile(x,0.975)),
                    l_ci_68 = apply(grad,2, function(x) quantile(x,0.16)), 
                    u_ci_68 = apply(grad,2, function(x) quantile(x,0.84)),
                    l_ci_90 = apply(grad,2, function(x) quantile(x,0.05)),
                    u_ci_90 = apply(grad,2, function(x) quantile(x,0.95)),
                    l_ci_99 = apply(grad,2, function(x) quantile(x,0.005)),
                    u_ci_99 = apply(grad,2, function(x) quantile(x,0.995)),
                    l_ci_80 = apply(grad,2, function(x) quantile(x,0.1)),
                    u_ci_80 = apply(grad,2, function(x) quantile(x,0.9))
                    
  )
  return(out)
  }
