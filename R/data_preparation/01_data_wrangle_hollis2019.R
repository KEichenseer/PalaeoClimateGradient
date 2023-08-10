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

dat <- data.frame(matrix(NA,nrow = sum(sapply(datr,nrow)),ncol = 8))
colnames(dat) <- c("Timeslice", "Age", "MgCa", "SST_median", "SST_025", "SST_975", "Species", "Preservation")
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
  dat$Preservation[(ind[i]+1):(ind[i+1])] <- datr[[i]]$Preservation
  
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
dat$Timeslice[which(is.na(as.numeric(dat$Age)) & !is.na(dat$Age) & is.na(dat$Timeslice))] <- 
  dat$Age[which(is.na(as.numeric(dat$Age)) & !is.na(dat$Age) & is.na(dat$Timeslice))]
dat$Age[which(is.na(as.numeric(dat$Age)))] <- NA 

datloc <- lapply(2:20,function(x) readxl::read_excel("data/raw/Eocene/Hollis2019SI/SDF06_TEX86.xlsx", sheet = x, skip = 0, trim_ws = T, n_max = 9))

locpar <- matrix(NA,nrow=length(datloc), ncol = 6)
for(i in 1:length(datloc)) {
locpar[i,] <- unlist(c(colnames(datloc[[i]])[2],datloc[[i]][1:5,2]))

}

# clean references
locpar[16,3] <- "Hollis et al. 2009; Hollis et al. 2012"
locpar[1,3] <- "Sluijs et al. 2006; 2008; 2009"    
locpar[9,3] <- "Keating-Bitonti et al. 2011" 
locpar[,3] <- gsub(",","",locpar[,3])
locpar[,2] <- gsub("Eastern USA","USA",locpar[,2]) # use USA to conform to the other data

locpar <- data.frame(locpar)
locpar[,4] <- as.numeric(locpar[,4])
locpar[,5] <- as.numeric(locpar[,5])

dat$Location <- NA
dat$Region <- NA
dat$Latitude <- NA
dat$Longitude <- NA
dat$Setting <- NA
dat$Reference <- NA

for(i in 1:nrow(locpar)) {
  dat$Location[(ind[i]+1):(ind[i+1])] <-  locpar[i,1]
  dat$Region[(ind[i]+1):(ind[i+1])] <-  locpar[i,2]
  dat$Latitude[(ind[i]+1):(ind[i+1])] <-  locpar[i,4]
  dat$Longitude[(ind[i]+1):(ind[i+1])] <-  locpar[i,5]
  dat$Setting[(ind[i]+1):(ind[i+1])] <-  locpar[i,6]
  dat$Reference[(ind[i]+1):(ind[i+1])] <-  locpar[i,3]
}

### EECO could be 53.3 - 49.2 Ma (according to the figures in Hollis et al 2019)
### However, data between range(as.numeric(dat$Age[which(dat$Timeslice=="EECO")]),na.rm=T)
### i.e. from 48.57246 to 53.48746 has been assigned to EECO in their database. Hence, we use that
# not needed anymore, age selection will be done in 02_prepare_hollis2019.R:
#dat$Timeslice[which(is.na(dat$Timeslice) & (as.numeric(dat$Age) >= 48.57246 | as.numeric(dat$Age) <= 53.48746))] <- "EECO"
# 

saveRDS(dat,"data/processed/Hollis_2019_TEX86_2022_07_18.rds")

#########################################################
### D47 data (already manually brought into a nice table)

d47 <- readxl::read_excel("data/raw/Eocene/Hollis2019SI/Hollis2019_D47_locations_and_values.xlsx")
saveRDS(d47,"data/processed/Hollis_2019_D47_2022_07_18.rds")

