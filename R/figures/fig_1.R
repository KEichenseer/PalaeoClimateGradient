# Load libraries --------------------------------------------------------------
library(raster)
library(tmap)
library(sf)
library(chronosphere)
# Read data -------------------------------------------------------------------
chem <- readRDS("./data/processed/Hollis_processed_EECO_2022_07_19.rds")
bio <- readRDS("./data/processed/bio_proxies_2022_08_08.RDS")
# PALEOMAP Palaeogeography
pgeog <- raster("./data/raw/DEM/early_eocene.tif")
shp <- fetch(dat = "paleomap", var = "paleocoastlines")
shp <- shp[11]
# Data prep -------------------------------------------------------------------
# Filter data
chem <- chem[, c("p_lng", "p_lat", "proxy")]
chem$shape <- 1
bio <- bio[, c("p_lng", "p_lat", "type")]
bio$shape <- 2
# Rename data
bio[which(bio$type == "reef"), c("type")] <- "Coral reef"
bio[which(bio$type == "mangrove"), c("type")] <- "Mangrove"
# Rename cols 
names(chem) <- c("p_lng", "p_lat", "Proxy type", "shape")
names(bio) <- c("p_lng", "p_lat", "Proxy type", "shape")
# Bind data
locs <- rbind.data.frame(chem, bio)
locs <- unique(locs)
locs$unique <- paste(locs$p_lng, locs$p_lat)
filt <- names(which(table(locs$uni) > 1))
locs[which(locs$uni %in% filt), c("Proxy type")] <- c("Various")
locs[which(locs$uni %in% filt), c("shape")] <- 3

# Set factors
locs[, c("Proxy type")] <- factor(locs[, c("Proxy type")],
                                  levels = c("Coral reef", "Mangrove", "d18O",
                                             "d47", "TEX86","MgCa", "Various"))
# Shift point for visualisation
locs[which.max(locs$p_lat), c("p_lat")] <- 86

# Reformat data for plotting
locs <- sf::st_as_sf(x = locs, coords = c("p_lng", "p_lat"), crs = 4326)
# Plot ------------------------------------------------------------------------
my_colours = rev(c(rgb(0.825,0.725,0.625), 
              rgb(0.9,0.8,0.7), 
              rgb(0.85,0.99,0.99), 
              rgb(0.6,0.85,0.2)))
# Create base map
m <- tm_shape(pgeog) +
  tm_raster(palette = "Blues",
            style = "pretty",
            legend.show = FALSE)

m <- m + tm_shape(shp) +
  tm_polygons(col = my_colours[4], alpha = 0.5)

m
m <- m + tm_shape(locs) +
  tm_dots(col = "Proxy type", palette = "-Dark2",
          shape = "Proxy type", shapes = c(21, 22, 23, 24, 25, 10, 8), 
          legend.shape.show = FALSE,
          legend.col.reverse = FALSE,
          border.col = "black", size = 0.15, alpha = 0.8,
          shapes.legend = c(21, 22, 23, 24, 25, 10, 8),
          legend.is.portrait = TRUE,
          title = "") +
  tm_legend(legend.position = c(0.9, 0.4),
            legend.frame = TRUE,
            legend.bg.color = "transparent") +
  tm_layout(frame = TRUE)

tmap_save(m, "./figures/fig_1.jpg",
          units = "mm",
          height = 110,
          width = 200,
          dpi = 600)


