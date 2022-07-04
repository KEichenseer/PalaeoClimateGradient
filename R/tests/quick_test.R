### Straightforward test

# load isotope data
iso <- readRDS("data/processed/StabisoDB_processed_26_06_2022.rds")
# load coral reef data
coral <- readRDS("data/processed/PARED_coralreefs_processed_02_07_2022.rds")

source("R/subscripts/AuxiliaryFunctions.R")


stage <- "Cenomanian"

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

mod3 <- climate_parallel(nChains = 3, nIter = 20000, obsmat = obsmat, distrmat = distrmat)
stopImplicitCluster()

mod <-  mod3[[1]]
### Plot
plot_gradient(mod, ylim = c(-5,40))
plot_data(obsmat,add = T)
plot_distr(distrmat)
plot_posterior(mod)

