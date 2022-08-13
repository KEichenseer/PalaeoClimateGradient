### combine the posteriors of a list of model outputs, optionally discarding burnin each time
combine_posterior <- function(mod, burnin = NULL) {
  n_thin <- mod[[1]]$call$n_thin
  n_iter <- floor(mod[[1]]$call$n_iter/n_thin)
  if(is.null(burnin)) burnin <- 0 else burnin <- floor(burnin/n_thin)
  out <- do.call(rbind,lapply(1:length(mod),function(f) mod[[f]]$params[(burnin+1):n_iter,]))
}
