

### Options for the Eocene climate model
#define options
n_chains = 4 # number of separate chains to run
#
n_iter = 100000 # number of iterations
# 
n_thin = 10 # thinning of the chains (save only every n_thinth iteration)
# 
adapt_sd = 5000 # number of iterations to adapt the proposal standard deviation
#
parallel = TRUE # should the chains run in parallel (using doParallel and foreach)
# 
# Define the priors
priors <- list(  # a list with one function per parameter of the logistic regression
  f1 = function(x,log) dsnorm(x,location = -3.0, scale = 12, alpha = 30, log = log), # prior on A (lower asymptote)
  f2 = function(x,lower,log) dtnorm(x, lower, upper = Inf, mean = 28, sd = 10, log = log), # prior on K (upper asymptote)
  f3 = function(x,log) dnorm(x, mean = 42, sd = 10, log = log), # prior on M (latitude of steepest point of the gradient)
  f4 = function(x,log) dgamma(x, shape = 4.3, rate = 30, log = log)# prior on B (growth rate)
)
#
# Define temperature distributions for Eocene ecological data
proxy_distributions <- data.frame(name = c("Avicennia", "Avicennia-Rhizophoraceae", "Reef"),
                                  distribution = c("normal", "normal", "normal"),
                                  mean = c(mean(c(15.6,22.5)), mean(c(20.7,29.5)), 27.6),
                                  sd = c((22.5-15.6)/4, c(29.5-20.7)/4, (29.5-21)/4),
                                  shape = rep(NA,3))
