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
chem$shape <- 23
bio <- bio[, c("lng", "lat", "type")]
bio$shape <- 21
# Rename data
bio[which(bio$type == "reef"), c("type")] <- "Coral reef"
bio[which(bio$type == "mangrove"), c("type")] <- "Mangrove"
# Rename cols 
names(chem) <- c("lng", "lat", "Proxy type", "shape")
names(bio) <- c("lng", "lat", "Proxy type", "shape")
# Bind data
locs <- rbind.data.frame(chem, bio)
locs <- unique(locs)
locs$uni <- paste(locs$lng, locs$lat)
filt <- names(which(table(locs$uni) > 1))
locs[which(locs$uni %in% filt), c("Proxy type")] <- c("Various")

# Set factors
locs[, c("Proxy type")] <- factor(locs[, c("Proxy type")],
                                  levels = c("Coral reef", "Mangrove", "d18O",
                                             "d47", "TEX86","MgCa", "Various"))
# Shift point for visualisation
locs[which.max(locs$lat), c("lat")] <- 86

# Reformat data for plotting
locs <- st_as_sf(x = locs, coords = c("lng", "lat"), crs = 4326)
# Plot ------------------------------------------------------------------------
# Create base map
m <- tm_shape(world) +
  tm_fill()  +
  tm_borders()
  
m <- m + tm_shape(locs) +
  tm_dots(col = "Proxy type", palette = "Dark2",
          shape = "shape", shapes = c(21, 21, 23, 23, 23, 23, 23), 
          legend.shape.show = FALSE,
          legend.col.reverse = FALSE,
          border.col = "black", size = 0.1, alpha = 0.7,
          shapes.legend = c(21, 21, 23, 23, 23, 23, 23)) +
  tm_legend(position = c(0.02, 0.1)) +
  tm_layout(frame = FALSE)

tmap_save(m, "./figures/data_distribution.jpg",
          units = "mm",
          height = 100,
          width = 200,
          dpi = 600)


