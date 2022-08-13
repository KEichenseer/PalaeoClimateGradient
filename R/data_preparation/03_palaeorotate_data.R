# Load functions --------------------------------------------------------------
source("./R/functions/palaeorotate.R")
# Read data -------------------------------------------------------------------
data_bio <- read.csv(file = "./data/raw/Eocene/bio_proxies.csv")
data_chem <- readRDS(file = "./data/processed/Hollis_processed_EECO_2022_07_19.rds")
# Rotate data -----------------------------------------------------------------
# Add age for rotation
data_bio$age <- 51.2
data_chem$age <- 51.2
# Update column names for consistency
colnames(data_chem)[which(colnames(data_chem) == "longitude")] <- "lng"
colnames(data_chem)[which(colnames(data_chem) == "latitude")] <- "lat"
# Update lng of point falling off plate
data_bio$lng[which(data_bio$location == "Site 343 (214.30-213.99 mbsf)")] <- 9
# Rotate data
data_bio <- palaeorotate(x = data_bio)
data_chem <- palaeorotate(x = data_chem)
# Save data
saveRDS(data_bio, "./data/processed/bio_proxies_2022_08_08.RDS")
saveRDS(data_chem, "./data/processed/Hollis_processed_EECO_2022_07_19.rds")
