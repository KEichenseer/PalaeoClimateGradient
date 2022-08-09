# Load libraries --------------------------------------------------------------
library(raster)
library(tmap)
library(sp)
# Read data -------------------------------------------------------------------
chem <- readRDS("./data/processed/Hollis_processed_2022_07_19.rds")
bio <- readRDS("./data/processed/bio_proxies_2022_08_08.RDS")
world <- read_sf(system.file("shapes/world.gpkg", package = "spData"))
# Data prep -------------------------------------------------------------------
chem <- chem[, c("longitude", "latitude", "proxy")]
bio <- bio[, c("lng", "lat", "type")]
# Rename data
bio[which(bio$type == "reef"), c("type")] <- "Coral reef"
bio[which(bio$type == "mangrove"), c("type")] <- "Mangrove"
# Rename cols 
names(chem) <- c("lng", "lat", "Proxy type")
names(bio) <- c("lng", "lat", "Proxy type")
# Bind data
locs <- rbind.data.frame(chem, bio)
locs <- unique(locs)
locs$uni <- paste(locs$lng, locs$lat)
filt <- names(which(table(locs$uni) > 1))
# Reformat data for plotting
locs <- st_as_sf(x = locs, coords = c("lng", "lat"), crs = 4326)
# Plot ------------------------------------------------------------------------
# Create base map
m <- tm_shape(world) +
  tm_fill()  +
  tm_borders()
  
m <- m + tm_shape(locs) +
  tm_dots(col = "Proxy type", shape = 21, border.col = "black", size = 0.1) +
  tm_legend(position = c(0.02, 0.1))

m

locs <- locs[which(locs$uni %in% filt),] 
locs[, c("Proxy type")] <- c("Various")

m <- m + tm_shape(locs) +
  tm_dots(col = "black", shape = 21, border.col = "black", size = 0.1) +
  tm_legend(position = c(0.02, 0.1))

m
