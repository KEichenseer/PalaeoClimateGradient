# Load libraries --------------------------------------------------------------
library(chronosphere)
# Read data -------------------------------------------------------------------
data_bio <- read.csv(file = "./data/raw/Eocene/bio_proxies.csv")
data_chem <- readRDS(file = "./data/processed/Hollis_processed_EECO_2022_07_19.rds")
# Rotate data -------------------------------------------------------------------
# Get plate model
pm <- fetch(dat = "paleomap", var = "model", 
           datadir = system.file("extdata", package = "chronosphere"))
# Rotate bio data
data_bio[, c("p_lng", "p_lat")] <- reconstruct(
  x = data_bio[, c("lng", "lat")],
  age = 51.2,
  model = pm)

# Rotate bio data
data_chem[, c("p_lng", "p_lat")] <- reconstruct(
  x = data_chem[, c("longitude", "latitude")],
  age = 51.2,
  model = pm)

# Save data
saveRDS(data_bio, "./data/processed/bio_proxies_2022_08_08.RDS")
saveRDS(data_chem, "./data/processed/Hollis_processed_EECO_2022_07_19.rds")
