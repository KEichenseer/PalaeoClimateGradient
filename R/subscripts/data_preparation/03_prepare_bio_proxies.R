# Load functions --------------------------------------------------------------
source("./R/functions/palaeorotate.R")
# Read data -------------------------------------------------------------------
data <- read.csv(file = "./data/raw/Eocene/bio_proxies.csv")
# Add rotation age
data$age <- 51.45
# Modify coordinates slightly off the shelf
data[which(data$location == "Site 343 (214.30-213.99 mbsf)"), c("lng")] <- 8.9
# Palaeorotate data
data <- palaeorotate(x = data, model = "Merdith2021")
# Save data
saveRDS(data, "./data/processed/bio_proxies_2022_08_08.RDS")
