### combine the posteriors of a list of model outputs, optionally discarding burnin each time
combine_posterior <- function(mod, burnin = NULL) {
  nThin <- mod[[1]]$call$nThin
  nIter <- floor(mod[[1]]$call$nIter/nThin)
  if(is.null(burnin)) burnin <- 0 else burnin <- floor(burnin/nThin)
  out <- do.call(rbind,lapply(1:length(mod),function(f) mod[[f]]$params[(burnin+1):nIter,]))
}
