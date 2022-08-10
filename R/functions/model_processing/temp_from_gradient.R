# Function to generate median temperatures and confidence intervals based on 
# a gradient from the model output

temp_from_gradient <- function(lat,model_out, burnin = NULL) {
  
  # load gradient function
  source("R/functions/model_components/gradient.R")
  
  if (is.null(burnin)) first_it <- 1 else first_it == burnin
  last_it <- nrow(model_out)
  
  grad <- gradient(lat,model_out[first_it:last_it,c("A","DKA","M","Q")],sdy=0)
  
  out <- data.frame(lat = lat, median = apply(grad,2,median), 
                    l_ci_95 = apply(grad,2, function(x) quantile(x,0.025)), 
                    u_ci_95 = apply(grad,2, function(x) quantile(x,0.975)))
  return(out)
  }
