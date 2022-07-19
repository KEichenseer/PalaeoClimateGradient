### Straightforward test
# check for missing packages and install if needed
list.of.packages <- c("cli", "doParallel","foreach","parallel","mvnfast","truncnorm")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


# load isotope data
iso <- readRDS("data/processed/StabisoDB_processed_26_06_2022.rds")
# load coral reef data
coral <- readRDS("data/processed/PARED_coralreefs_processed_02_07_2022.rds")

source("R/subscripts/AuxiliaryFunctions.R")


stage <- "Ypresian"

obsmat <- iso_obsmat(iso,stage)
distrmat <- coral_distrmat(coral,stage)

#### Climate_parallel
####
source("R/subscripts/ClimateGradientModel.R")

source("R/subscripts/ClimateParallel.R")

library(foreach)
library(doParallel)

cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)

mod1 <- climate_parallel(nChains = 3, nIter = 20000, obsmat = obsmat, distrmat = NULL)
parallel::stopCluster(cl)


### check chains
plot_chains(mod1)
par(mfrow=c(1,1),mar = c(4,4,1,1))

### Plot result
mod <-  mod1[[1]]
plot_gradient(mod, ylim = c(-5,40))
plot_data(obsmat,add = T)
plot_distr(distrmat)
plot_posterior(mod)

