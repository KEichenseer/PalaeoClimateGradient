### d18O data

#datr <- gdata::read.xls("data/raw/Eocene/Hollis2019SI/SDF03_d18O.xlsx", sheet = 3, pattern = "Site")
datr <- lapply(3:24,function(x) readxl::read_excel("data/raw/Eocene/Hollis2019SI/SDF03_d18O.xlsx", sheet = x, skip = 11, trim_ws = T))

for(i in 1:length(datr)) {
  if(all(is.na(datr[[i]][1,]))) datr[[i]] <- datr[[i]][2:nrow(datr[[i]]),]
  datr[[i]] <- datr[[i]][1:(which(sapply(1:nrow(datr[[i]]), function(x) all(is.na(datr[[i]][x,]))))[1]-1),]
}

for(i in 1:length(datr)) {
  colnames(datr[[i]])[1:8] <- c("Site", "Hole", "Core", "Section", "Interval_cm", "Depth", "Timeslice", "GTS2012")
  if(i %in% c(1:6,8:11,13:14,16:22)) {colnames(datr[[i]])[9:19] <- c("d18O_measured", "d18Osw_adj", "Temperature", "Species", 
                                                                   "N_specimen", "Size_fraction_um", "Depth_habitat", 
                                                                   "Preservation", "Pres_assessment", "Source", 
                                                                   "Comment")
  datr[[i]]$d18Osw_not_adj <- NA
  datr[[i]]$SE_d18O <- NA
  datr[[i]]$d13C_measured <- NA
  
  
  }
  if(i ==7) { colnames(datr[[i]])[9:20] <- c("d18O_measured","d18Osw_not_adj", "d18Osw_adj", "Temperature", "Species", 
                                                                    "N_specimen", "Size_fraction_um", "Depth_habitat", 
                                                                    "Preservation", "Pres_assessment", "Source", 
                                                                    "Comment")
  datr[[i]]$SE_d18O <- NA
  datr[[i]]$d13C_measured <- NA
  }
  if(i==12){ colnames(datr[[i]])[9:21] <- c("d18O_measured","d18Osw_not_adj", "d18Osw_adj", "Temperature", "Species", 
                                           "N_specimen", "Size_fraction_um", "Depth_habitat", 
                                           "Preservation", "Pres_assessment", "Source", "SE_d18O",
                                           "Comment")
  
  datr[[i]]$d13C_measured <- NA
  }
  if(i==15) {colnames(datr[[i]])[9:20] <- c("d13C_measured", "d18O_measured", "d18Osw_adj", "Temperature", "Species", 
                                           "N_specimen", "Size_fraction_um", "Depth_habitat", 
                                           "Preservation", "Pres_assessment", "Source", 
                                           "Comment")
  datr[[i]]$d18Osw_not_adj <- NA
  datr[[i]]$SE_d18O <- NA
  }
}
  
col_order <- colnames(datr[[1]])
for(i in 1:length(datr)) {
  datr[[i]] <- datr[[i]][,col_order]
}
test <- sapply(datr,function(x) colnames(x))
apply(test,1,function(y) all(sapply(1:22, function(z) identical(y[1],y[z]))))
## good, all the columns have identical names and are in the right order

## add latitudes, longitudes, location, reference, setting 

locations <- readxl::read_excel("data/raw/Eocene/Hollis2019SI/Hollis2019_d18O_locations.xlsx")
colnames(locations) %in% colnames(dat)

for(i in 1:length(datr)) {
  datr[[i]]$Location <- rep(locations$Location[i],nrow(datr[[i]]))
  datr[[i]]$Region <- rep(locations$Region[i],nrow(datr[[i]]))
  datr[[i]]$Latitude <- rep(locations$Latitude[i],nrow(datr[[i]]))
  datr[[i]]$Longitude <- rep(locations$Longitude[i],nrow(datr[[i]]))
  datr[[i]]$Palaeolat_ori <- rep(locations$Palaeolat_ori[i],nrow(datr[[i]]))
  datr[[i]]$Setting <- rep(locations$Setting[i],nrow(datr[[i]]))
  datr[[i]]$Reference <- rep(locations$Reference[i],nrow(datr[[i]]))
}

dat <- do.call(rbind,datr)
dat$Temperature <- as.numeric(dat$Temperature)

saveRDS(dat,"data/processed/Hollis_2019_dO18_2022_07_15.rds")

### Mg/Ca data
datr <- lapply(2:12,function(x) readxl::read_excel("data/raw/Eocene/Hollis2019SI/SDF04_Mg-Ca.xlsx", sheet = x, skip = 10, trim_ws = T))

for(i in 1:length(datr)) {
  if(all(is.na(datr[[i]][1,]))) datr[[i]] <- datr[[i]][2:nrow(datr[[i]]),]
  datr[[i]] <- datr[[i]][1:(which(sapply(1:nrow(datr[[i]]), function(x) all(is.na(datr[[i]][x,]))))[1]-1),]
}

dat <- data.frame(matrix(NA,nrow = sum(sapply(datr,nrow)),ncol = 7))
colnames(dat) <- c("Timeslice", "Age", "MgCa", "SST_median", "SST_025", "SST_975", "Species")
ind <- c(0,cumsum(sapply(datr,nrow)))
for(i in 1:length(datr)) {
  if("Age" %in% colnames(datr[[i]])) {
    dat$Age[(ind[i]+1):(ind[i+1])] <- datr[[i]]$Age
  }
  if("Time-slice" %in% colnames(datr[[i]])) {
    
    dat$Timeslice[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`Time-slice`
  }
  dat$MgCa[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`Mg/Ca`
  dat$SST_median[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`50`
  dat$SST_025[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`2.5`
  dat$SST_975[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`97.5`
  dat$Species[(ind[i]+1):(ind[i+1])] <- datr[[i]]$Species
  
}  

dat$Timeslice[which(is.na(as.numeric(dat$Age)))] <- dat$Age[which(is.na(as.numeric(dat$Age)))]
dat$Age[which(is.na(as.numeric(dat$Age)))] <- NA 

## add latitudes, longitudes, location, reference, setting 

locations <- readxl::read_excel("data/raw/Eocene/Hollis2019SI/Hollis2019_MgCa_locations.xlsx")
colnames(locations) %in% colnames(dat)
dat$Location <- NA
dat$Region <- NA
dat$Latitude <- NA
dat$Longitude <- NA
dat$Palaeolat_ori <- NA
dat$Setting <- NA
dat$Reference <- NA

for(i in 1:length(datr)) {
  dat$Location[(ind[i]+1):(ind[i+1])] <- rep(locations$Location[i],nrow(datr[[i]]))
  dat$Region[(ind[i]+1):(ind[i+1])] <- rep(locations$Region[i],nrow(datr[[i]]))
  dat$Latitude[(ind[i]+1):(ind[i+1])] <- rep(locations$Latitude[i],nrow(datr[[i]]))
  dat$Longitude[(ind[i]+1):(ind[i+1])] <- rep(locations$Longitude[i],nrow(datr[[i]]))
  dat$Palaeolat_ori[(ind[i]+1):(ind[i+1])] <- rep(locations$Palaeolat_ori[i],nrow(datr[[i]]))
  dat$Setting[(ind[i]+1):(ind[i+1])] <- rep(locations$Setting[i],nrow(datr[[i]]))
  dat$Reference[(ind[i]+1):(ind[i+1])] <- rep(locations$Reference[i],nrow(datr[[i]]))
}
saveRDS(dat,"data/processed/Hollis_2019_MgCa_2022_07_17.rds")


### TEX86 data
datr <- lapply(2:20,function(x) readxl::read_excel("data/raw/Eocene/Hollis2019SI/SDF06_TEX86.xlsx", sheet = x, skip = 11, trim_ws = T))

for(i in 1:length(datr)) {
  if(all(is.na(datr[[i]][1,]))) datr[[i]] <- datr[[i]][2:nrow(datr[[i]]),]
  datr[[i]] <- datr[[i]][1:(which(sapply(1:nrow(datr[[i]]), function(x) all(is.na(datr[[i]][x,]))))[1]-1),]
}

table(unlist(sapply(datr,colnames)))
(sapply(datr,function(x) "Interval" %in% colnames(x)))

dat <- data.frame(matrix(NA,nrow = sum(sapply(datr,nrow)),ncol = 6))
colnames(dat) <- c("Timeslice", "Age", "TEX", "SST_median", "SST_05", "SST_95")
ind <- c(0,cumsum(sapply(datr,nrow)))
for(i in 1:length(datr)) {
  if("Age" %in% colnames(datr[[i]])) {
    dat$Age[(ind[i]+1):(ind[i+1])] <- datr[[i]]$Age
  }
  if("Age (GTS2012)" %in% colnames(datr[[i]])) {
    dat$Age[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`Age (GTS2012)`
  }
  if("Age (Ma)" %in% colnames(datr[[i]])) {
    dat$Age[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`Age (Ma)`
  }
  if("Time slice" %in% colnames(datr[[i]])) {
    
    dat$Timeslice[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`Time slice`
  }
  if("T-slice" %in% colnames(datr[[i]])) {
    
    dat$Timeslice[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`T-slice`
  }
  if("Interval" %in% colnames(datr[[i]])) {
    
    dat$Timeslice[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`Interval`
  }
  if("GTS2012(B)" %in% colnames(datr[[i]])) {
    dat$Age[(ind[i]+1):(ind[i+1])] <- apply(cbind(datr[[i]]$`GTS2012(B)`,datr[[i]]$`GTS2012(D)`),1,mean)
  }
  
  dat$TEX[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`TEX`
  dat$SST_median[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`0.5`
  dat$SST_05[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`0.05`
  dat$SST_95[(ind[i]+1):(ind[i+1])] <- datr[[i]]$`0.95`

}  
View(dat)
dat$Timeslice[which(is.na(as.numeric(dat$Age)))] <- dat$Age[which(is.na(as.numeric(dat$Age)))]
dat$Age[which(is.na(as.numeric(dat$Age)))] <- NA 

