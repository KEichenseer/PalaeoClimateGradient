

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
  f1 = function(x,log) dsnorm(x,location = -3.03, scale = 12, alpha = 30, log = log), # prior on A (lower asymptote)
  f2 = function(x,lower,log) dtnorm(x, lower, upper = Inf, mean = 28.3, sd = 10, log = log), # prior on K (upper asymptote)
  f3 = function(x,log) dnorm(x, mean = 42.1, sd = 10, log = log), # prior on M (latitude of steepest point of the gradient)
  f4 = function(x,log) dgamma(x, shape = 4.3, rate = 30, log = log)# prior on B (growth rate)
)
#
# 
