### Read data
## Hollis data
dat <- readRDS("data/processed/Hollis_processed_EECO_2022_07_19.rds")
#
# prepare for use in the model
obsmat <- data.frame(sample = (paste(abs(dat$p_lat),dat$longitude, dat$proxy)),
                     p_lat = abs(dat$p_lat), temperature = dat$temperature,
                     sd = dat$temperature_sd,
                     proxy = dat$proxy)
# order observation matrix for easily keeping tracks of mean estimates
obsmat <- obsmat[with(obsmat, order(p_lat, sample)),]

plot(obsmat$p_lat,obsmat$sd)
summary(lm(obsmat$sd~obsmat$p_lat))

loc_sd <- sapply(unique(obsmat$sample), function(x) sd(obsmat$temperature[which(obsmat$sample==x)]))
loc_lat <- sapply(unique(obsmat$sample), function(x) unique(obsmat$p_lat[which(obsmat$sample==x)]))
plot(loc_lat,loc_sd, ylab = "sd of temperature estimates", )
