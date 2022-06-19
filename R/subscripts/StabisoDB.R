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

get_temperature <- function(d18O_sample, d18O_sw_ice, lat, mineral, belemnite) {
  out = rep(NA_real_,length(d18O_sample))
  
  out[which(mineral == "calcite" & "belemnite" == F)] <- calcite_temp(
    d18O_sample, d18O_sw_ice + lat_correction(lat,d18O_sw_ice))
  out[which(mineral == "calcite" & "belemnite" == T)] <- calcite_temp(
    d18O_sample - 1.5, d18O_sw_ice + lat_correction(lat,d18O_sw_ice))
  out[which(mineral == "aragonite")] <- calcite_temp(
    d18O_sample - 0.6, d18O_sw_ice + lat_correction(lat,d18O_sw_ice))
  out[which(mineral == "apatite")] <- apatite_temp(
    d18O_sample, d18O_sw_ice + lat_correction(lat,d18O_sw_ice))
return(out)
  }

# An example with ID 16333 (a sample from Grossman and Joachmiski 2022)
get_temperature(d18O_sample = c(-2.58), d18O_sw_ice=c(0), lat = c(7.99), mineral=c("aragonite"), belemnite=c(F))
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

### assign temperatures (use select d18O for now)