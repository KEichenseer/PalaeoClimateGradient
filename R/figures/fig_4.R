# Load libraries ---------------------------------------------------------
library(ggplot2)
library(ggpubr)
source("./R/functions/model_processing/temp_from_gradient.R")
source("./R/functions/model_processing/combine_posterior.R")
# Load data --------------------------------------------------------------
eeco_mod <- readRDS("./results/eeco/eeco_climate_model_output.rds")
bio <- readRDS("./data/processed/distribution_matrix_for_EECO_model.rds")
chem <- readRDS("./data/processed/observation_matrix_for_EECO_model.rds")
modern_sample <- readRDS("./results/modern/modern_climate_model_output.rds")

# Data preparation -------------------------------------------------------
# Filter chem data
chem <- chem[, c("p_lat", "temperature", "sd", "proxy")]
chem$type <- "Geochemical"
# Filter eco data
bio <- bio[, c("p_lat", "mu", "scale")]
names(bio) <- c("p_lat", "temperature", "sd")
bio$type <- "Ecological"
bio$proxy <- "Mangrove"
# Update reef occurrences
bio$proxy[which(bio$sd == 2.125)] <- "Coral Reef"
# Bind data
locs <- rbind.data.frame(chem, bio)

# Set factors
levels <- c("Coral Reef", "Mangrove", "d18O", "d47", "TEX86", "MgCa")
locs$proxy <- factor(locs$proxy, levels = levels)

# Rename data for plotting
locs[which(locs$proxy == "d18O"), c("type")] <- expression(delta^18*"O")
locs[which(locs$proxy == "d47"), c("type")] <- expression(Delta[47])
locs[which(locs$proxy == "TEX86"), c("type")] <- expression(TEX[86])
locs[which(locs$proxy == "MgCa"), c("type")] <- "Mg/Ca"

# eeco model
eeco_mod <- combine_posterior(mod = eeco_mod, burnin = 100000)
# Calculate eeco temperatures
lat <- seq(from = 0, to = 90, by = 0.1)
eeco_temp <- temp_from_gradient(lat = lat, model_out = eeco_mod)
modern_sample <- temp_from_gradient(lat = lat, model_out = modern_sample)

# Plot data --------------------------------------------------------------
# Labels
labs <- c("Coral reef", "Mangrove", 
          expression(delta^18*"O"), 
          expression(Delta[47]), 
          expression(TEX[86]),
          "Mg/Ca")
# Modern gradient
ggplot() + 
  geom_point(data = locs,
             aes(x = p_lat, y = temperature, shape = proxy, fill = proxy),
             colour = "black", size = 2) +
  geom_line(data = modern_sample,
            aes(x = lat, y = median, colour = "Modern gradient"),
            size = 1.5, alpha = 0.85) +
  geom_ribbon(data = eeco_temp,
              aes(x = lat, ymin = l_ci_95, ymax = u_ci_95),
              colour = NA, fill = "#54278f", alpha = 0.35) +
  geom_line(data = eeco_temp,
             aes(x = lat, y = median, colour = "EECO gradient"),
             size = 1.5, alpha = 0.85) +
  xlab("Absolute Latitude (\u00B0)") +
  ylab("Sea surface temperature (\u00B0C)") +
  scale_x_continuous(limits = c(0, 90), labels = seq(0, 90, 15), breaks = seq(0, 90, 15)) +
  scale_colour_manual(values = c("#02818a", "#54278f"),
                     breaks = c("Modern gradient",
                                "EECO gradient")) + 
  scale_shape_manual(values = 20:26, labels = labs) +
  scale_fill_manual(values = 1:6, labels = labs) +
  theme_bw() +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = c(0.2, 0.25),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 12),
        legend.spacing.y = unit(0.5, "mm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold"),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank())

ggsave("./figures/fig_4_revised.png", units = "mm", width = 150, height = 150, dpi = 600)
