# Load libraries --------------------------------------------------------------
library(raster)
library(ggplot2)
library(ggpubr)
# Load data -------------------------------------------------------------------
temp <- raster(
  "./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")
modern_sample <- readRDS("./results/modern/modern_sample_gradient.RDS")
# Data preparation ------------------------------------------------------------
# Extract temperature data
temp <- raster::as.data.frame(x = temp, xy = TRUE, centroids = TRUE)
colnames(temp) <- c("lng", "lat", "SST")
# Generate random sample
#temp <- temp[sample(1:nrow(temp), size = 1000),]
# Remove data from continents
temp <- na.omit(temp)
# Make lats absolute values
temp$lat <- abs(temp$lat)
# Add 1 deg latitudinal bins
temp$bin <- ceiling(temp$lat) - 0.5
# Calculate mean per lat bin
gradient <- tapply(temp$SST, temp$bin, median)
gradient <- data.frame(lat = as.numeric(names(gradient)), 
                       SST = as.vector(gradient))
# Plot ------------------------------------------------------------------------
# Modern gradient
p1 <- ggplot(data = temp, aes(x = lat, y = SST)) + 
  geom_point(shape = 21, 
             colour = "black",
             fill = "#02818a", 
             size = 1, 
             alpha = 0.4) +
  geom_line(data = gradient, aes(x = lat, y = SST),
            size = 1,
            colour = "black",
            alpha = 0.8) +
  xlab("Latitude (\u00B0)") +
  scale_x_continuous(limits = c(0, 90),
                     breaks = seq(0, 90, 30),
                     labels = seq(0, 90, 30)) +
  ylab("Temperature (\u00B0C)") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank(),
  )
# Modern model with random sample (n = 10,000)
p2 <- ggplot(data = gradient, aes(x = lat, y = SST)) + 
  geom_line(size = 1, colour = "black", alpha = 0.4) +
  geom_line(size = 1, 
            colour = "#1d91c0") +
  xlab("Latitude (\u00B0)") +
  scale_x_continuous(limits = c(0, 90),
                     breaks = seq(0, 90, 30),
                     labels = seq(0, 90, 30)) +
  ylab("Temperature (\u00B0C)") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank(),
  )
p2
# Modern model using Eocene sample distribution (n = XXX)
p3 <- ggplot(data = gradient, aes(x = lat, y = SST)) + 
  geom_line(size = 1, colour = "black", alpha = 0.4) +
  geom_line(size = 1, 
            colour = "#1d91c0") +
  xlab("Latitude (\u00B0)") +
  scale_x_continuous(limits = c(0, 90),
                     breaks = seq(0, 90, 30),
                     labels = seq(0, 90, 30)) +
  ylab("Temperature (\u00B0C)") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank(),
  )
p3

ggsave("./figures/fig_3.jpg",
       plot = p1,
       units = "mm", width = 150, height = 150,
       dpi = 600)
  
  