### Straightforward test
# check for missing packages and install if needed
list.of.packages <- c("cli", "doParallel","foreach","parallel","mvnfast","truncnorm")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


# load isotope data from Hollis 2019
iso <- readRDS("data/processed/Hollis_2019_dO18_2022_07_15.rds")

source("R/subscripts/AuxiliaryFunctions.R")


stage <- "Ypresian"

iso_sub <-  subset(iso,Timeslice %in% c("EECO","EECO?","late EECO") & Depth_habitat == "Mixed-layer")
iso_sub <- iso_sub[with(iso_sub, order(abs(Palaeolat_ori), Longitude)),]


# prepare for use in the model
data_mat <- data.frame(sample = (paste(abs(iso_sub$Palaeolat_ori),iso_sub$Longitude)),
                       latitude = abs(iso_sub$Palaeolat_ori), temperature = iso_sub$Temperature)
distrmat <- NULL # coral_distrmat(coral,stage)

#### Climate_parallel
####
source("R/subscripts/ClimateGradientModel.R")

source("R/subscripts/ClimateParallel.R")

library(foreach)
library(doParallel)

cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)

mod34 <- climate_parallel(nChains = 3, nIter = 20000, obsmat = data_mat, distrmat = NULL)
parallel::stopCluster(cl)


### check chains
plot_chains(mod3)
par(mfrow=c(1,1),mar = c(4,4,1,1))

### Plot result
mod <-  mod3[[1]]
plot_gradient(mod, ylim = c(-5,40))
plot_data(obsmat,add = T)
plot_distr(distrmat)
plot_posterior(mod)

plot_gradient(mod34[[1]],line_col = rgb(.8,0,0.5,.8), confint_col = rgb(.8,0,0.5,0.2), add = T)

plot_gradient(mod10c[[1]],line_col = rgb(0,0.75,0,.8), confint_col = rgb(0,0.75,0,0.125), add = T)

nIter = 20000
nChains = 3
obsmat = data_mat
distrmat = NULL
coeff_inits = NULL
sdy_init = NULL

yest_inits = NULL
sdyest_inits = NULL

proposal_var_inits = c(2,2,2,0.2)
adapt_sd = floor(0.1 * nIter)

adapt_sd_decay = max(floor(0.005*nIter),1)
quiet = FALSE

