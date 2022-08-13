# Load libraries --------------------------------------------------------------
library(raster)
library(rgdal)
library(broom)
library(ggplot2)
# Read data -------------------------------------------------------------------
chem <- readRDS("./data/processed/Hollis_processed_EECO_2022_07_19.rds")
bio <- readRDS("./data/processed/bio_proxies_2022_08_08.RDS")
# Merdith2021 shapefiles
coastline <-
  "http://gws.gplates.org/reconstruct/coastlines/?time=51&model=MERDITH2021"
coastline <- rgdal::readOGR(coastline)
coast_poly <- broom::tidy(coastline)
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
names(chem) <- c("p_lng", "p_lat", "proxy", "shape")
names(bio) <- c("p_lng", "p_lat", "proxy", "shape")
# Unique data
bio <- unique(bio)
chem <- unique(chem)
# Bind data
locs <- rbind.data.frame(chem, bio)
locs$unique <- paste(locs$p_lng, locs$p_lat)
filt <- names(which(table(locs$uni) > 1))
locs[which(locs$uni %in% filt), c("proxy")] <- c("Various")
locs[which(locs$uni %in% filt), c("shape")] <- 3

# Set factors
locs[, c("proxy")] <- factor(locs[, c("proxy")],
                                  levels = c("Coral reef", "Mangrove", "d18O",
                                             "d47", "TEX86","MgCa", "Various"))
# Shift point for visualisation
locs[which.max(locs$p_lat), c("p_lat")] <- 86

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
  ggthemes::theme_map() +
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0, 0.8),
        legend.key.size = unit(x = 5, units = "mm"),
        legend.text.align = 0)

eocene_map <- eocene_map +
  geom_point(data = locs,
             aes(x = p_lng,
                 y = p_lat,
                 shape = proxy,
                 colour = proxy,
                 fill = proxy)) + 
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 19, 20),
                     labels = c("Coral reef", "Mangrove", expression(delta^18*"O"),
                                expression(Delta[47]), expression(TEX[86]), "Mg/Ca",
                                "Various")) +
  scale_fill_manual(values = c("#33a02c",
                                 "#1f78b4",
                                 "#e31a1c",
                                 "#6a3d9a",
                                 "#b15928",
                                 "#ff7f00",
                                 "#fb9a99"),
                    labels = c("Coral reef", "Mangrove", expression(delta^18*"O"),
                               expression(Delta[47]), expression(TEX[86]), "Mg/Ca",
                               "Various")) +
  scale_colour_manual(values = c("black",
                               "black",
                               "black",
                               "black",
                               "black",
                               "#ff7f00",
                               "black"),
                      labels = c("Coral reef", "Mangrove", expression(delta^18*"O"),
                                 expression(Delta[47]), expression(TEX[86]), "Mg/Ca",
                                 "Various")) +
  guides(shape=guide_legend(ncol=2))



ggsave("./figures/fig_1.png", plot = eocene_map, units = "mm",
       width = 200, height = 100, dpi = 600)

# display the figure in the R plot window
eocene_map

