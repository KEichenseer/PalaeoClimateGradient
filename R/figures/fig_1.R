# Load libraries --------------------------------------------------------------
library(raster)
library(tmap)
library(sf)
library(broom)
# Read data -------------------------------------------------------------------
chem <- readRDS("./data/processed/Hollis_processed_EECO_2022_07_19.rds")
bio <- readRDS("./data/processed/bio_proxies_2022_08_08.RDS")
# Merdith2021 shapefiles
coastline <-
  "http://gws.gplates.org/reconstruct/coastlines/?time=51&model=MERDITH2021"
polygon <-
  "http://gws.gplates.org/reconstruct/static_polygons/?time=51&model=MERDITH2021"
coastline <- rgdal::readOGR(coastline)
coast_poly <- broom::tidy(coastline)
polygon <- rgdal::readOGR(polygon)
plate_poly <- broom::tidy(polygon)
#Data prep -------------------------------------------------------------------
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
#locs <- sf::st_as_sf(x = locs, coords = c("p_lng", "p_lat"), crs = 4326)
# Plot ------------------------------------------------------------------------
eocene_map <-
  ggplot() +
  geom_map(
    data = coast_poly, map = coast_poly,
    aes(x = long, y = lat, map_id = id),
    size = 0.15, fill = "grey80", colour = "grey80"
  ) +
  geom_rect(
    data = data.frame(xmin = -180, xmax = 180, ymin = -90, ymax = 90),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    color = 1, fill = NA, size = 0.3
  ) +
  coord_map("mollweide") +
  ggthemes::theme_map()
eocene_map +
  geom_point(
    data = locs,
    aes(x = p_lng,
        y = p_lat,
        shape = `Proxy type`,
        colour = `Proxy type`,
        fill = `Proxy type`)) + 
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 19, 20)) +
  scale_fill_manual(values = c("#33a02c",
                                 "#1f78b4",
                                 "#e31a1c",
                                 "#6a3d9a",
                                 "#b15928",
                                 "#ff7f00",
                                 "#fb9a99")) +
  scale_colour_manual(values = c("#33a02c",
                               "#1f78b4",
                               "#e31a1c",
                               "#6a3d9a",
                               "#b15928",
                               "#ff7f00",
                               "#fb9a99"))




my_colours = rev(c(rgb(0.825,0.725,0.625), 
              rgb(0.9,0.8,0.7), 
              rgb(0.85,0.99,0.99), 
              rgb(0.6,0.85,0.2)))
# Create base map
m <- tm_shape(coast_poly) +
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


