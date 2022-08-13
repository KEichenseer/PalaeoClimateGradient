# Load libraries --------------------------------------------------------------
library(raster)
# Set options
reps = 100
# Read data -------------------------------------------------------------------
temp <- raster(
  "./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")
chem <- readRDS("./data/processed/Hollis_processed_EECO_2022_07_19.rds")
bio <- readRDS("./data/processed/bio_proxies_2022_08_08.RDS")
# Data processing -------------------------------------------------------------
# Reduce resolution of raster
r <- raster(res = 1)
temp <- resample(x = temp, y = r)
# Unique locations
chem <- unique(chem[, c("p_lng", "p_lat", "proxy")])
bio <- unique(bio[, c("p_lng", "p_lat", "type")])
# Rename cols to conform
names(chem) <- c("p_lng", "p_lat", "type")
# Bind data
locs <- rbind.data.frame(chem, bio)
# Round off data
locs[, c("p_lng", "p_lat")] <- round(locs[, c("p_lng", "p_lat")], digits = 2)
locs <- unique(locs)
# Order by latitude
locs <- locs[order(locs$p_lat),]
# Rename rows
row.names(locs) <- 1:nrow(locs)
# Sample data
samples <- lapply(1:nrow(locs), function(i) {
  # Set extent values
  xmin <- -180
  xmax <- 180
  ymin <- floor(locs$p_lat[i])
  ymax <- ceiling(locs$p_lat[i])
  # Set extent for sampling
  e <- raster::extent(c(xmin, xmax, ymin, ymax))
  # Sample raster for given extent
  raster::sampleRandom(x = temp, 
                       ext = e,
                       size = reps,
                       na.rm = TRUE)
})
# Median and range calculation ------------------------------------------------
#Create empty dataframe
df <- data.frame(p_lat = locs$p_lat,
                 median = rep(NA, nrow(locs)),
                 min = rep(NA, nrow(locs)),
                 max = rep(NA, nrow(locs)))
# Run over number of samples
for(i in 1:length(samples)){
  df$median[i] <- median(samples[[i]])
  df$min[i] <- min(samples[[i]])
  df$max[i] <- max(samples[[i]])
}
saveRDS(df, "./results/modern/modern_sample_summary_stats.RDS")
# Reformat --------------------------------------------------------------------
# Create empty list
l <- list()
# Run over number of reps
for(i in 1:reps){
  # Generate empty matrix
  m <- matrix(nrow = nrow(locs), ncol = 2)
  # Add lats
  m[, 1] <- locs$p_lat 
  # Add sampled temp
  m[, 2] <- sapply(X = 1:length(samples), function(x){
    samples[[x]][i]
  })
  # Add column names
  colnames(m) <- c("p_lat", "temp")
  # Add to list
  l[[i]] <- m
}
# Save data
saveRDS(l, "./results/modern/modern_sample.RDS")

