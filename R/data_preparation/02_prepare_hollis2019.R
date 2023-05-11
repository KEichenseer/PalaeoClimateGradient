d18O <- readRDS("data/processed/Hollis_2019_dO18_2022_07_15.rds")
MgCa <- readRDS("data/processed/Hollis_2019_MgCa_2022_07_17.rds")
d47 <- readRDS("data/processed/Hollis_2019_D47_2022_07_18.rds")
TEX86 <- readRDS("data/processed/Hollis_2019_TEX86_2022_07_18.rds")

head(d18O)

sum(table(d18O$Timeslice))/nrow(d18O)
sum(table(MgCa$Timeslice))/nrow(MgCa)
sum(table(d47$Age))/nrow(d47)
sum(table(TEX86$Timeslice))/nrow(TEX86)

range(as.numeric(d18O$GTS2012[which(d18O$Timeslice=="EECO?")]),na.rm=T)

### assign d18O data to EECO if they are within the range of 49.11000 53.76049
### simplify to 49.1 - 53.8
d18O$EECO <- rep(0,nrow(d18O))
d18O$EECO[which(d18O$GTS2012 >= 49.1 & d18O$GTS2012 <=53.8 &!is.na(as.numeric(d18O$GTS2012)))] <- 1
table(d18O$Timeslice[which(is.na(d18O$GTS2012))])
## we should have captured all EECO data

### assign MgCa data to EECO if they are within the range of 49.1 - 53.8
MgCa$EECO <- rep(0,nrow(MgCa))
MgCa$EECO[which(MgCa$Age >= 49.1 & MgCa$Age <=53.8 & !(is.na(as.numeric(MgCa$Age))))] <- 1
table(MgCa$Timeslice[which(is.na(MgCa$Age))])
## we should have captured all EECO data

### assign d47 data to EECO if they are within the range of 49.1 - 53.8
d47$EECO <- rep(0,nrow(d47))
d47$EECO[which(d47$Age >= 49.1 & d47$Age <=53.8)] <- 1
table(MgCa$Timeslice[which(is.na(MgCa$Age))])
## all of them are EECO

### assign TEX86 data to EECO if they are within the range of 49.1 - 53.8
TEX86$EECO <- rep(0,nrow(TEX86))
TEX86$EECO[which(TEX86$Age >= 49.1 & TEX86$Age <=53.8 & !is.na(as.numeric(TEX86$Age)))] <- 1
TEX86$EECO[which(is.na(as.numeric(TEX86$Age)) & TEX86$Timeslice=="EECO")] <- 1
table(TEX86$Age[which(TEX86$Timeslice=="EECO" & TEX86$EECO==1)])
TEX86$Age[which(TEX86$EECO==1)]
## we should have captured all EECO data

### collate the data in one nice dataframe
dat <- data.frame(
location = c(d18O$Site, MgCa$Location, d47$Location, TEX86$Location),
region = c(d18O$Region, MgCa$Region, d47$Region, TEX86$Region),
latitude = c(d18O$Latitude, MgCa$Latitude, d47$Latitude, TEX86$Latitude),
longitude = c(d18O$Longitude, MgCa$Longitude, d47$Longitude, TEX86$Longitude),
palaeolat_ori = c(d18O$Palaeolat_ori, MgCa$Palaeolat_ori, d47$Palaeolat_ori, rep(NA,nrow(TEX86))),
EECO = c(d18O$EECO,MgCa$EECO,d47$EECO,TEX86$EECO),
Age = c(d18O$GTS2012,MgCa$Age,d47$Age,TEX86$Age),
proxy = c(rep("d18O",(nrow(d18O))),rep("MgCa",(nrow(MgCa))), 
               rep("d47",(nrow(d47))), rep("TEX86",(nrow(TEX86)))),
proxy_value = c(d18O$d18O_measured,MgCa$MgCa,d47$D47,TEX86$TEX),
d18Osw_adj = c(d18O$d18Osw_adj,rep(NA,nrow(MgCa)+nrow(d47)+nrow(TEX86))),
organism = c(rep("planktic foraminifera",(nrow(d18O))), rep("planktic foraminifera",(nrow(MgCa))),
             d47$Sample_type, rep("Nitrososphaerota-derived lipids", nrow(TEX86))),
species = c(d18O$Species, MgCa$Species, d47$Species,  rep("unknown", nrow(TEX86))),
preservation = c(d18O$Preservation, MgCa$Preservation, rep(NA, nrow(d47)),  rep(NA, nrow(TEX86))),
setting = c(d18O$Setting,MgCa$Setting,d47$Setting, TEX86$Setting),
depth_habitat = c(d18O$Depth_habitat,rep(NA,nrow(MgCa)), rep(NA,nrow(d47)), rep(NA,nrow(TEX86))),
temperature = c(d18O$Temperature, MgCa$SST_median, d47$Temp, TEX86$SST_median),
temperature_detail = c(rep("estimate",nrow(d18O)), rep("median and 2.5 and 97.5 percentiles", nrow(MgCa)),
                       rep("mean value and SD?",nrow(d47)), rep("mean and 5 and 95 percentiles",nrow(TEX86))),
temperature_lower = c(rep(NA,nrow(d18O)), MgCa$SST_025,
                      d47$Temp-2*d47$Temp_unc_neg, TEX86$SST_05),
temperature_upper = c(rep(NA,nrow(d18O)), MgCa$SST_975,
                      d47$Temp+2*d47$Temp_unc_pos, TEX86$SST_95),
temperature_sd =  c(rep(NA,nrow(d18O)),rep(NA,nrow(MgCa)),
                    (d47$Temp_unc_neg+d47$Temp_unc_pos)/2, rep(NA,nrow(TEX86)))

)

# correct a few wrong upper temperatures from the MgCa data
dat$temperature_upper[which(dat$temperature_upper == dat$temperature)] <-
  (2*dat$temperature-dat$temperature_lower)[which(dat$temperature_upper == dat$temperature)]

# check what the upper and lower temperatures look like
sub <- subset(dat,proxy == "TEX86")
sub <- sub[order(sub$temperature),]
plot(sub$temperature)
points(sub$temperature_lower, col = rgb(0.75,0,0,0.5), type = "o", lwd = 2)
points(sub$temperature_upper, col = rgb(0.75,0,0,0.5), type = "o", lwd = 2)

plot(sub$temperature_upper-sub$temperature, type = "l", lwd = 2)
points(sub$temperature_upper-sub$temperature, type = "l", lwd = 2, col = rgb(1,0,0,0.5))

### Ok so the Quantile ranges are symmetrical --> we use the quantiles directly as mean and derive SD from them

# MgCa offers 2.5 and 97.5 percentiles, so the range between 97.5 and 50 denotes 1.96 * SD
dat$temperature_sd[which(dat$proxy=="MgCa")] <- 
  (dat$temperature_upper-dat$temperature)[which(dat$proxy=="MgCa")]/1.96

# TEX86 offers 5 and 95 percentiles, so the range between 95 and 50 denotes 1.645 * SD
dat$temperature_sd[which(dat$proxy=="TEX86")] <- 
  (dat$temperature_upper-dat$temperature)[which(dat$proxy=="TEX86")]/1.645

# fix alternative spellings, make numeric
dat$proxy_value <- as.numeric(dat$proxy_value)
dat$preservation[which(dat$preservation=="Recrystallized")] <- "recrystallised"
dat$preservation[which(dat$preservation=="Recrystallization")] <- "recrystallised"


# select data from one stage to test, exclude NA data
data_sub <- subset(dat,EECO == 1 & !(is.na(temperature))  &
                     (is.na(depth_habitat) | depth_habitat %in% c("Mixed-layer",     "Mixed layer")))
# number of d18O samples
table(data_sub$proxy)


# exclude recrystallised data
data_sub <- subset(data_sub,preservation != "recrystallised" | proxy !="d18O")

saveRDS(data_sub,"data/processed/Hollis_processed_EECO_2022_07_19.rds")

