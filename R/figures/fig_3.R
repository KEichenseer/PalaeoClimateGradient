# Load libraries --------------------------------------------------------------
library(raster)
library(ggplot2)
library(ggpubr)
source("./R/functions/model_processing/temp_from_gradient.R")
source("./R/functions/model_processing/combine_posterior.R")
# Load data -------------------------------------------------------------------
temp <- raster(
  "./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")
modern_sample <- readRDS("./results/modern/modern_sample_gradient.RDS")
eocene_sample <- readRDS("./results/modern/modern_sample_eocene_p_lat_gradient.RDS")
eocene_reps <- readRDS("./results/modern/modern_sample_summary_stats.RDS")
# Data preparation ------------------------------------------------------------
# Reduce resolution of raster
r <- raster(res = 1)
temp <- resample(x = temp, y = r)
# Extract temperature data
temp <- raster::as.data.frame(x = temp, xy = TRUE, centroids = TRUE)
colnames(temp) <- c("lng", "lat", "SST")
# Remove data from continents
temp <- na.omit(temp)
# Make lats absolute values
temp$lat <- abs(temp$lat)
eocene_reps$p_lat <- abs(eocene_reps$p_lat)
# Add 1 deg latitudinal bins
temp$bin <- ceiling(temp$lat) - 0.5
# Calculate mean per lat bin
grad <- tapply(temp$SST, temp$bin, median)
grad <- data.frame(lat = as.numeric(names(grad)), 
                       SST = as.vector(grad))

# Calculate modern temperatures from sample
lat <- seq(from = 0, to = 90, by = 0.1)
modern_sample <- temp_from_gradient(lat = lat, model_out = modern_sample)
# Calculate modern temperatures from palaeodistribution of Eocene data
eocene_sample <- combine_posterior(mod = eocene_sample, burnin = 5000)
eocene_sample <- temp_from_gradient(lat = lat, model_out = eocene_sample)
# Plot ------------------------------------------------------------------------
# Modern gradient
p1 <- ggplot() + 
  geom_point(data = temp,
             aes(x = lat,
                     y = SST),
                     colour = "grey80",
             size = 1, 
             alpha = 0.7) +
  geom_line(data = grad,
            aes(x = lat,
                y = SST,
                colour = "Empirical gradient"),
            size = 1, 
            alpha = 1) +
  geom_line(data = modern_sample,
            aes(x = lat, 
                y = median,
                colour = "Modelled gradient"),
            size = 1.5, 
            alpha = 0.85) +
  xlab("Latitude (\u00B0)") +
  scale_x_continuous(limits = c(0, 90),
                     breaks = seq(0, 90, 30),
                     labels = seq(0, 90, 30)) +
  ylab("Sea surface temperature (\u00B0C)") +
  theme_bw() +
  theme(
    legend.position = c(0.8, 0.9),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank(),
  )
p1 <- p1 + scale_color_manual(values = c("black", "#02818a"))

# Modern model using Eocene sample distribution
p2 <- ggplot() + 
  geom_pointrange(data = eocene_reps,
                  aes(x = p_lat,
                      y = median,
                      ymin = min,
                      ymax = max),
                  colour = "grey80",
                  size = 1,
                  shape = 21,
                  fatten = 1,
                  alpha = 0.7) +
  geom_line(data = grad,
            aes(x = lat,
                y = SST,
                colour = "Empirical gradient"),
            size = 1, 
            alpha = 1) +
  geom_line(data = eocene_sample,
            aes(x = lat, 
                y = median,
                colour = "Modelled gradient"),
            size = 1.5, 
            alpha = 0.7) +
  xlab("Latitude (\u00B0)") +
  scale_x_continuous(limits = c(0, 90),
                     breaks = seq(0, 90, 30),
                     labels = seq(0, 90, 30)) +
  ylab("Sea surface temperature (\u00B0C)") +
  theme_bw() +
  theme(
    legend.position = c(0.8, 0.9),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank(),
  )
p2 <- p2 + scale_color_manual(values = c("black", "#54278f"))
# Combine plots and save ------------------------------------------------------
p <- ggarrange(p1, p2,
               ncol = 1,
               labels = "auto")

ggsave("./figures/fig_3.jpg",
       plot = p,
       units = "mm", width = 150, height = 250,
       dpi = 600)
  
  