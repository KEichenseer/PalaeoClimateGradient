### Function to generate climate gradient across all stages

### Does not quite work yet

## vector with stages
###### 
stages <- c("Pre-Cambrian boundary", "Fortunian", "Stage 2", "Stage 3", "Stage 4", "Wuliuan", "Drumian", 
            "Guzhangian", "Paibian", "Jiangshanian", "Stage 10", "Tremadocian", 
            "Floian", "Dapingian", "Darriwilian", "Sandbian", "Katian", "Hirnantian", 
            "Rhuddanian", "Aeronian", "Telychian", "Sheinwoodian", "Homerian", 
            "Gorstian", "Ludfordian", "Pridoli", "Lochkovian", "Pragian", 
            "Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", 
            "Visean", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", 
            "Gzhelian", "Asselian", "Sakmarian", "Artinskian", "Kungurian", 
            "Roadian", "Wordian", "Capitanian", "Wuchiapingian", "Changhsingian", 
            "Induan", "Olenekian", "Anisian", "Ladinian", "Carnian", "Norian", 
            "Rhaetian", "Hettangian", "Sinemurian", "Pliensbachian", "Toarcian", 
            "Aalenian", "Bajocian", "Bathonian", "Callovian", "Oxfordian", 
            "Kimmeridgian", "Tithonian", "Berriasian", "Valanginian", "Hauterivian", 
            "Barremian", "Aptian", "Albian", "Cenomanian", "Turonian", "Coniacian", 
            "Santonian", "Campanian", "Maastrichtian", "Danian", "Selandian", 
            "Thanetian", "Ypresian", "Lutetian", "Bartonian", "Priabonian", 
            "Rupelian", "Chattian", "Aquitanian", "Burdigalian", "Langhian", 
            "Serravallian", "Tortonian", "Messinian", "Zanclean", "Piacenzian", 
            "Gelasian", "Calabrian", "Chibanian", "Upper Pleistocene", "Holocene")
######

# load isotope data
iso <- readRDS("data/processed/StabisoDB_processed_26_06_2022.rds")

# load coral reef data
coral <- readRDS("data/processed/PARED_coralreefs_processed_02_07_2022.rds")


# prepare parallel running
source("R/subscripts/ClimateParallel.R")

library(foreach)
library(doParallel)

nChains = 3
nClusters = nChains
cl <- parallel::makeCluster(nClusters)
doParallel::registerDoParallel(cl)

# Iterations and burnin
nIter = 100000
burnin = 12000
iterOut = 4000
seqstep <- round((nIter-burnin+1)/iterOut)

nstage <- 10

# list to save output
mod <- vector(mode = "list", length = nstage)

cli::cli_progress_bar('Sampling', total = nstage)

# loop across stages
for(st in 10) {
  # start progress bar
  ### The MCMC loop

    # select stage
  mystage <- stages[st]
  
  # update progress bar
  cli::cli_progress_update(set = st, status = paste0('stage ', st, ' - ', mystage, ' - '))
  
  ##
  ## subset isotope data
  ##
  iso_sub <- subset(iso,stage_2020 == mystage & !(is.na(paleolat)) & !(is.na(temperature)))#table(iso$stage_2020)
  iso_sub <- iso_sub[with(iso_sub, order(abs(paleolat), longitude)),]
  # prepare for use in the model
  iso_mod <- data.frame(sample = (paste(abs(iso_sub$paleolat),iso_sub$longitude)),
                        latitude = abs(iso_sub$paleolat), temperature = iso_sub$temperature)
  samples1 <- unique(iso_mod$sample)
  for(i in 1:length(samples1)) iso_mod$sample[which(iso_mod$sample==samples1[i])] <- i
  
  nsamples <- length(samples1)
  samples <-  1:(max(1,nsamples))
  sample_lats <- sapply(samples, function(x) unique(iso_mod$lat[which(iso_mod$sample==x)]))
  
  ##
  ## subset coral reef data
  ##
  coral_sub <- subset(coral,early_stage == mystage & !(is.na(pal_lat_scotese)))
  coral_sub <- coral_sub[with(coral_sub, order(abs(pal_lat_scotese), longit)),]
  # prepare for use in the model
 if(nrow(coral_sub) >= 1) {coral_distrmat <- data.frame(latitude = abs(coral_sub$pal_lat_scotese),
                               location = 22.8,
                               scale = 10,
                               shape = 4,
                               distribution = "skew-normal")
 } else coral_distrmat <- data.frame(NULL)
  
  # run model in parallel chains
  smod <- climate_parallel(nChains = nChains, nIter = nIter, obsmat = iso_mod, distrmat = coral_distrmat)
  
  # store only part of the output to save memory
  mod[[st]] <- do.call("rbind",lapply(1:nChains,function(l) 
    smod[[l]]$params[seq(burnin+seqstep,nIter,seqstep),1:5]))
}

doParallel::stopImplicitCluster()

### check effective sample size
# should aim for ~ 10,000 for the publication
mcmcse::multiESS(mod[[10]][,1:4])

### plot output
source("R/subscripts/AuxiliaryFunctions.R")
plot_gradient(mod[[10]])

