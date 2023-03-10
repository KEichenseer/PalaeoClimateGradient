# Load libraries ---------------------------------------------------------
library(sf)
library(ggplot2)
library(ggthemes)
# Read data --------------------------------------------------------------
chem <- readRDS("./data/processed/Hollis_processed_EECO_2022_07_19.rds")
bio <- readRDS("./data/processed/bio_proxies_2022_08_08.RDS")
# Merdith2021 shapefiles
coastline <- "https://gws.gplates.org/reconstruct/coastlines/?time=51&model=MERDITH2021"
coastline <- sf::read_sf(coastline)
#Data prep ---------------------------------------------------------------
# Filter data
chem <- chem[, c("p_lng", "p_lat", "proxy")]
chem$shape <- 1
bio <- bio[, c("p_lng", "p_lat", "type")]
bio$shape <- 2
# Rename cols 
names(chem) <- c("p_lng", "p_lat", "proxy", "shape")
names(bio) <- c("p_lng", "p_lat", "proxy", "shape")
# Unique data
bio <- unique(bio)
chem <- unique(chem)
# Bind data
locs <- rbind.data.frame(chem, bio)

# Set factors
levels <- c("reef", "mangrove", "d18O", "d47", "TEX86", "MgCa")
locs$proxy <- factor(locs$proxy, levels = levels)

# Rename data for plotting
locs[which(locs$type == "reef"), c("type")] <- "Coral reef"
locs[which(locs$type == "mangrove"), c("type")] <- "Mangrove"
locs[which(locs$type == "d18O"), c("type")] <- expression(delta^18*"O")
locs[which(locs$type == "d47"), c("type")] <- expression(Delta[47])
locs[which(locs$type == "TEX86"), c("type")] <- expression(TEX[86])
locs[which(locs$type == "MgCa"), c("type")] <- "Mg/Ca"

locs <- sf::st_as_sf(x = locs, coords = c("p_lng", "p_lat"), crs = sf::st_crs(4326))

# Plot ------------------------------------------------------------------------
# Labels
labs <- c("Coral reef", "Mangrove", 
          expression(delta^18*"O"), 
          expression(Delta[47]), 
          expression(TEX[86]),
          "Mg/Ca")
# Shapes
#shps <- c(21, 22, 23, 24, 25, 19, 20)
# Plot
ggplot() +
  geom_sf(data = coastline, colour = "grey80", fill = "grey80") +
  geom_sf(data = locs, aes(colour = proxy, fill = proxy, shape = proxy),
          size = 2, alpha = 0.85) +
  coord_sf(crs = sf::st_crs("ESRI:54030")) +
  scale_shape_manual(values = 20:26, labels = labs) +
  scale_colour_manual(values = rep("black", 6), labels = labs) +
  scale_fill_manual(values = 1:6, labels = labs) +
  guides(shape = guide_legend(nrow = 1, override.aes = list(size = 3.5))) +
  theme_map() +
  theme(plot.margin = margin(0.5, 0.5, 1, 0.5, "cm"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.2, -0.075),
        legend.key.size = unit(x = 5, units = "mm"),
        legend.text = element_text(size = 10),
        legend.text.align = 0)
  
ggsave("./figures/fig_1.png", units = "mm",
       width = 200, height = 120, dpi = 600, bg = "white")
