# Load libraries --------------------------------------------------------------
library(raster)
library(tmap)
library(sp)
# Read data -------------------------------------------------------------------
chem <- readRDS("./data/processed/Hollis_processed_EECO_2022_07_19.rds")
bio <- readRDS("./data/processed/bio_proxies_2022_08_08.RDS")
# Cao et al 2017 Palaeogeog:
# https://www.earthbyte.org/improving-global-paleogeography-since-the-late-paleozoic-using-paleobiology/
pgeog <- raster("./data/raw/palaeogeog/PaleogeogRecon_Matthews2016_53Ma.tiff")
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
locs <- st_as_sf(x = locs, coords = c("p_lng", "p_lat"), crs = 4326)
# Plot ------------------------------------------------------------------------
my_colors = c("#9ecae1", "black", "#c6dbef", "#bababa")
# Create base map
m <- tm_shape(pgeog) +
      tm_raster(style = "cat", 
                palette = my_colors,
                legend.show = FALSE)

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
            legend.bg.color = "transparent")
  tm_layout(frame = TRUE)

tmap_save(m, "./figures/fig_1.jpg",
          units = "mm",
          height = 110,
          width = 200,
          dpi = 600)


