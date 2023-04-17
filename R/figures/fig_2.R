# Load libraries --------------------------------------------------------
library(ggplot2)
library(ggthemes)
# Read data -------------------------------------------------------------
df <- readRDS("./results/simulation/simulated_gradients_temperature_estimates.rds")
# Remove modern gradient
df <- subset(df, gradient_type != "modern")
# Remove 50 sample size
df <- subset(df, sample_size != 50)
# Rename sample size
df$sample_size[which(df$sample_size == "eocene")] <- "34 (Eocene samples)"
# Rename gradients
df$gradient_type[which(df$gradient_type == "extreme_greenhouse")] <- "Extreme Greenhouse"
df$gradient_type[which(df$gradient_type == "greenhouse")] <- "Greenhouse"
df$gradient_type[which(df$gradient_type == "icehouse")] <- "Icehouse"
df$gradient_type[which(df$gradient_type == "extreme_icehouse")] <- "Extreme Icehouse"
# Set factor levels
df$gradient_type <- factor(x = df$gradient_type, levels = c("Extreme Icehouse",
                                                            "Icehouse",
                                                            "Greenhouse",
                                                            "Extreme Greenhouse"))
df$sample_size <- factor(x = df$sample_size, levels = c(5, 10, 20, "34 (Eocene samples)"))
# Plot data -------------------------------------------------------------
col <- "#02818a"
ggplot(data = df, aes(x = lat, y = median)) +
  geom_ribbon(aes(x = lat, ymin = p_100, ymax = p_900),
              fill = col, alpha = 0.2) +
  geom_ribbon(aes(x = lat, ymin = p_005, ymax = p_995),
              fill = col, alpha = 0.2) +
  geom_ribbon(aes(x = lat, ymin = p_025, ymax = p_975),
              fill = col, alpha = 0.2) +
  geom_line(aes(x = lat, y = known), linetype = 1, colour = "black", size = 0.75, alpha = 1) +
  geom_line(colour = col, linetype = 1, size = 0.75, alpha = 0.9) +
  facet_grid(gradient_type~sample_size) +
  scale_x_continuous(limits = c(0, 90), labels = seq(0, 90, 15), breaks = seq(0, 90, 15)) +
  labs(y = "Sea surface temperature (\u00B0C)", x = "Absolute latitude (\u00B0)") +
  theme_bw() +
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position = c(0.8, 0.9),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12, face = "bold"),
        axis.ticks = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank())
# Save plot -------------------------------------------------------------
ggsave("./figures/fig_2.png",
       units = "mm", width = 200, height = 200, dpi = 600)
