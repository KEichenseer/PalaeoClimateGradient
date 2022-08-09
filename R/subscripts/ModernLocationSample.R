# Load libraries --------------------------------------------------------------
library(raster)
# Set options
reps = 100
# Read data -------------------------------------------------------------------
temp <- raster(
  "./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")
chem <- readRDS("./data/processed/Hollis_processed_2022_07_19.rds")
bio <- readRDS("./data/processed/bio_proxies_2022_08_08.RDS")
# Data processing -------------------------------------------------------------
# Unique locations
chem <- unique(chem[, c("longitude", "latitude")])
bio <- unique(bio[, c("lng", "lat")])
# Rename cols to conform
names(chem) <- c("lng", "lat")
# Bind data
locs <- rbind.data.frame(chem, bio)
# Round off data
locs <- round(locs, digits = 2)
locs <- unique(locs)
# Order by latitude
locs <- locs[order(locs$lat),]
# Rename rows
row.names(locs) <- 1:nrow(locs)
# Sample data
samples <- lapply(1:nrow(locs), function(i) {
  # Set extent values
  xmin <- -180
  xmax <- 180
  ymin <- floor(locs$lat[i])
  ymax <- ceiling(locs$lat[i])
  # Set extent for sampling
  e <- raster::extent(c(xmin, xmax, ymin, ymax))
  # Sample raster for given extent
  raster::sampleRandom(x = temp, 
                       ext = e,
                       size = reps,
                       na.rm = TRUE)
})
# Reformat --------------------------------------------------------------------
# Create empty list
l <- list()
# Run over number of reps
for(i in 1:reps){
  # Generate empty matrix
  m <- matrix(nrow = nrow(locs), ncol = 2)
  # Add lats
  m[, 1] <- locs$lat 
  # Add sampled temp
  m[, 2] <- sapply(X = 1:length(samples), function(x){
    samples[[x]][i]
  })
  # Add column names
  colnames(m) <- c("lat", "temp")
  # Add to list
  l[[i]] <- m
}
# Save data
saveRDS(l, "./results/modern/modern_sample.RDS")

