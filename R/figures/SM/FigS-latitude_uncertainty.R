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



png("figures/SM/FigS_temp_location_sd.png", width = 3.5, height = 3.4, units = "in", res = 400)
par(mar = c(4.1,4.1,1,1), las = 1)
plot(loc_lat,loc_sd, ylab = "sd of proxy temperatures", xlab =expression("|latitude| ("*degree*")"), pch = 21,
     col = NA, bg = rgb(0,0,0,0.5))
abline(lm(loc_sd~loc_lat))
sumlm <- summary(lm(loc_sd~loc_lat))
sumlm$r.squared
sumlm$coefficients[2,4]

legend("topleft",legend = bquote("R"^2~"="~.(round(sumlm$r.squared,2))*","~"p ="~.(round(sumlm$coefficients[2,4],3))), bty = "n",
       cex = 0.96)
dev.off()