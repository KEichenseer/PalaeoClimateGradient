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

