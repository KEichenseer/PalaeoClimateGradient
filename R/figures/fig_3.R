# Load libraries ---------------------------------------------------------
library(terra)
library(ggplot2)
library(ggpubr)
source("./R/functions/model_processing/temp_from_gradient.R")
source("./R/functions/model_processing/combine_posterior.R")
# Load data --------------------------------------------------------------
temp <- rast("./data/raw/climate/BioOracle_20220711/Present.Surface.Temperature.Mean.asc")
modern_sample <- readRDS("./results/modern/modern_climate_model_output.rds")
eocene_sample <- readRDS("./results/modern/modern_sample_eocene_p_lat_gradient.RDS")
eocene_reps <- readRDS("./results/modern/modern_sample_summary_stats.RDS")
# Data preparation -------------------------------------------------------
# Reduce resolution of raster
r <- rast(res = 1)
temp <- resample(x = temp, y = r)
# Extract temperature data
temp <- terra::as.data.frame(x = temp, xy = TRUE)
colnames(temp) <- c("lng", "lat", "SST")
# Remove data from continents
temp <- na.omit(temp)
# Make lats absolute values
temp$lat <- abs(temp$lat)
eocene_reps$p_lat <- abs(eocene_reps$p_lat)
# Add 1 deg latitudinal bins
temp$bin <- ceiling(temp$lat) - 0.5
# Calculate median per lat bin
grad <- tapply(temp$SST, temp$bin, median)
grad <- data.frame(lat = as.numeric(names(grad)), SST = as.vector(grad))
# save for later use
saveRDS(grad,"results/modern/empirical_median_1-deg-lat.rds")
# Calculate modern temperatures from sample
lat <- seq(from = 0, to = 90, by = 0.1)
modern_sample <- temp_from_gradient(lat = lat, model_out = modern_sample)
# Calculate modern temperatures from palaeodistribution of Eocene data
#eocene_sample <- combine_posterior(mod = eocene_sample, burnin = 5000)
#eocene_sample <- temp_from_gradient(lat = lat, model_out = eocene_sample)
# Plot ------------------------------------------------------------------------
# Modern gradient
p1 <- ggplot() + 
   geom_point(data = temp,
              aes(x = lat, y = SST), 
              colour = "grey40", size = 0.7, alpha = 0.075) +
  # stat_density2d(data = temp, aes(x = lat, y = SST, fill = sqrt(..density..)), geom = "tile", contour = FALSE,
  #                show.legend=FALSE) +
  # scale_fill_gradient(low = "transparent", high = "grey20", guide = "none") + 
  # guides(fill = "none",  color = guide_legend(override.aes = list(fill = NA))) +  # Remove fill legend
  geom_line(data = grad,
            aes(x = lat, y = SST, colour = "Empirical gradient"),
            size = 1.15, alpha = 1) +
  geom_line(data = modern_sample,
            aes(x = lat, y = median, colour = "Modelled gradient \n(full distribution)"),
            size = 1.6, alpha = 0.85) +
  scale_color_manual(values = c("black", "#02818a")) +
  scale_x_continuous(limits = c(0, 90),
                     breaks = seq(0, 90, 15),
                     labels = seq(0, 90, 15)) +
  xlab("Absolute latitude (\u00B0)") +
  ylab("Sea surface temperature (\u00B0C)") +
  theme_bw() +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = c(0.8, 0.9),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(), # Remove the key background
        legend.key.size = unit(1, 'lines'),
        legend.box.background = element_blank(),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 12),
        legend.spacing.y = unit(0.5, "mm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank())
p1
# # Modern model using Eocene sample distribution
# p2 <- ggplot() + 
#   geom_pointrange(data = eocene_reps,
#                   aes(x = p_lat, y = median, ymin = min, ymax = max),
#                   colour = "#54278f", size = 0.7, fatten = 1, alpha = 0.35) +
#   geom_ribbon(data = eocene_sample,
#               aes(x = lat, ymin = l_ci_95, ymax = u_ci_95),
#               fill = "#54278f", colour = NA, alpha = 0.3) +
#   geom_line(data = grad,
#             aes(x = lat, y = SST, colour = "Empirical gradient"),
#             size = 1, alpha = 1) +
#   geom_line(data = eocene_sample,
#             aes(x = lat, y = median, colour = "Modelled gradient \n(sampled distribution)"),
#             size = 1.5, alpha = 0.9) +
#   xlab("Abs. latitude (\u00B0)") +
#   ylab("Sea surface temperature (\u00B0C)") +
#   scale_color_manual(values = c("black", "#54278f")) +
#   scale_x_continuous(limits = c(0, 90),
#                      breaks = seq(0, 90, 30),
#                      labels = seq(0, 90, 30)) +
#   theme_bw() +
#   theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
#         legend.position = c(0.8, 0.9),
#         legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.text = element_text(size = 12),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14, face = "bold"),
#         axis.ticks = element_line(colour = "grey70", size = 0.2),
#         panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
#         panel.grid.minor = element_blank())
# # Combine plots and save ------------------------------------------------------
# ggarrange(p1, p2, ncol = 1, labels = "AUTO", font.label = list(size = 18))
ggsave("./figures/fig_3_revised.png", p1,
       units = "mm", width = 150, height = 150, dpi = 600, scale = 1)

