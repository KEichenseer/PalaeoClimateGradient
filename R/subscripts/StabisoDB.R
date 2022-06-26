###
## Transfer function for calcite from Kim and O'Neil (1997),
#  re-expressed in Leng and Marshall (2004):
#
calcite_temp <- function(d18O_calcite,d18O_sw) {
  13.8 - 4.58 * (d18O_calcite - d18O_sw) + 0.08 * (d18O_calcite - d18O_sw)^2
}

calcite_temp(-1.6,-0.78)
calcite_temp(-1.17,-0.67)

## belemnite correction: -1.5
## aragonite correction: +0.6

## latitude correction
lat_correction <- function(lat, d18O_sw_ice) {
  lat = abs(lat)
  d18O_sw_lat_mod = -6.650 * 10^(-4)*lat^2 + 3.363 * 10^(-2) * lat + 0.174
  d18O_sw_lat_pg = -4.944 * 10^(-4)*lat^2 + 2.492 * 10^(-2) * lat + 0.102
  out = d18O_sw_lat_mod * (d18O_sw_ice--1.08)/(0.45--1.08) + d18O_sw_lat_pg * (1-(d18O_sw_ice--1.08)/(0.45--1.08))
  return(out)
}

lat_correction(40.21,-1.08)
# test - seems to be roughly correct, according to the table from Grossman and Joachmiski 2022

###
## Transfer function for phosphate from Puceat et al. (2010),
#  which version?
#
apatite_temp <- function(d18O_apatite,d18O_sw) {
118.7 - 4.22 * (d18O_apatite - d18O_sw)
}
# do we need to correct for standard d18O in the apatite samples? Ask Grossman / Joachimski

get_temperature <- function(d18O_sample, d18O_sw_ice, lat, mineral, fossil_group) {
  out = rep(NA_real_,length(d18O_sample))
  
  index1 <- which(mineral == "calcite" & fossil_group != "Belemnite")
  if(length(index1)>=1) out[index1] <- calcite_temp(
    d18O_sample[index1], d18O_sw_ice[index1] + lat_correction(lat[index1],d18O_sw_ice[index1]))
  
  index2 <- which(mineral == "calcite" & fossil_group == "Belemnite")
  if(length(index2)>=1) out[index2] <- calcite_temp(
    d18O_sample[index2]- 1.5, d18O_sw_ice[index2] + lat_correction(lat[index2],d18O_sw_ice[index2]))
  
  index3 <- which(mineral == "aragonite" & fossil_group != "Belemnite")
  if(length(index3)>=1) out[index3] <- calcite_temp(
    d18O_sample[index3]- 0.6, d18O_sw_ice[index3] + lat_correction(lat[index3],d18O_sw_ice[index3]))
  
  index4 <- which(mineral == "apatite")
  if(length(index4)>=1) out[index4] <- apatite_temp(
    d18O_sample[index4], d18O_sw_ice[index4] + lat_correction(lat[index4],d18O_sw_ice[index4]))
 
return(out)
  }

# An example with ID 16333 (a sample from Grossman and Joachmiski 2022)
get_temperature(d18O_sample = c(-2.58), d18O_sw_ice=c(0), lat = c(7.99), mineral=c("aragonite"), fossil_group="something else")
# it looks like Grossman and Joachimski 2022 have not actually applied the -0.6 aragonite correction.
# Also, my temperature values are generally 0.1 - 0.2 deg C warmer than theirs, I can't figure out why.


###
## Read StabisoDB and clean
#
# stab1 <- sapply(0:11,function(x) openxlsx::read.xlsx(paste("data/raw/download (",x,").xlsx",sep = ""))) 
# fails so read individually

dat0 <- openxlsx::read.xlsx(paste("data/raw/download (",0,").xlsx",sep = ""))
dat1 <- openxlsx::read.xlsx(paste("data/raw/download (",1,").xlsx",sep = ""))
dat2 <- openxlsx::read.xlsx(paste("data/raw/download (",2,").xlsx",sep = ""))
dat3 <- openxlsx::read.xlsx(paste("data/raw/download (",3,").xlsx",sep = ""))
dat4 <- openxlsx::read.xlsx(paste("data/raw/download (",4,").xlsx",sep = ""))
dat5 <- openxlsx::read.xlsx(paste("data/raw/download (",5,").xlsx",sep = ""))
dat6 <- openxlsx::read.xlsx(paste("data/raw/download (",6,").xlsx",sep = ""))
dat7 <- openxlsx::read.xlsx(paste("data/raw/download (",7,").xlsx",sep = ""))
dat8 <- openxlsx::read.xlsx(paste("data/raw/download (",8,").xlsx",sep = ""))
dat9 <- openxlsx::read.xlsx(paste("data/raw/download (",9,").xlsx",sep = ""))
dat10 <- openxlsx::read.xlsx(paste("data/raw/download (",10,").xlsx",sep = ""))
dat11 <- openxlsx::read.xlsx(paste("data/raw/download (",11,").xlsx",sep = ""))

datlist  <- list(dat0,dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10,dat11)
for(i in 1:12) {
  
  colnames(datlist[[i]])[c(44,51:54)]<- c("d13c_permille","d47_permille", "d47_method", "d47_laboratory", "d47_comment")
  
}
dat <- do.call(rbind, datlist)
rm(dat0,dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10,dat11)

table(dat$period) # looks like it worked

### Assign seawater ice corrected d18O to the samples 
icetab <- read.csv("data/raw/seawater_ice_d18O_grossman_joachimski_2022.csv")
dat$sw_d18O_ice <- rep(NA,nrow(dat))
for(i in 1:nrow(icetab)) {
  dat$sw_d18O_ice[which(dat$age < icetab$start_age[i] & dat$age >=icetab$end_age[i])] <- icetab$d18O[i]
}
hist(dat$sw_d18O_ice)
hist(dat$age[which(dat$sw_d18O_ice==0.45)],1000)
# looks ok

### correct mineralogies:
dat$mineralogy[which(dat$mineralogy=="Aragonite")] <- "aragonite"
dat$mineralogy[which(dat$mineralogy=="Calcite")] <- "calcite"
dat$mineralogy[which(dat$mineralogy=="na")] <- NA

### add column for belemnite
table(dat$fossil_group)

### assign temperatures (use select d18O for now)
dat$temperature <- get_temperature(dat$select_d18o_permille, dat$sw_d18O_ice, dat$paleolat, dat$mineralogy, dat$fossil_group)
hist(dat$temperature)
### assign stages (use palaeoverse::time_bins)
dat$stage_2020[which(dat$stage_2020 == "Ionian")] <- "Chibanian"
dat$stage_2020[which(dat$stage_2020 == "Tarantian")] <- "Upper Pleistocene"
dat$stage_2020[which(dat$stage_2020 == "Cambrian")] <- "Pre-Cambrian boundary"
dat$stage_2020[which(dat$stage_2020 == "Age 3")] <- "Stage 2" # from boundary pack in prev interval (is empty anyway)
dat$stage_2020[which(dat$stage_2020 == "Ordovician")] <- "Stage 10" # from boundary pack in prev interval (is empty anyway)

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
unistage <- unique(dat$stage_2020)
unistage[which(!(unistage %in% stages))] # Ludlow Devonian Modern - not use


#get_temperature(d18O_sample = -0.82,d18O_sw_ice = -0.83,lat = -0, mineral = "calcite", fossil_group = "brachiopod")
#small_dat <- dat[,c("stage_2020","fossil_group","mineralogy","select_d18o_permille","sw_d18O_ice","temperature", "paleolat") ]

saveRDS(dat,"data/processed/StabisoDB_processed_26_06_2022.rds")
