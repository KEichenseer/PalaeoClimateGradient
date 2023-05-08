# Load libraries --------------------------------------------------------
library(ggplot2)
library(ggthemes)
library(ggtext)
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

# Get R2 values
r2 <- readRDS("./results/simulation/simulated_gradients_r_squared_gelman_unweighted.rds")
# Rename gradients
r2$gradient_type[which(r2$gradient_type == "extreme_greenhouse")] <- "Extreme Greenhouse"
r2$gradient_type[which(r2$gradient_type == "greenhouse")] <- "Greenhouse"
r2$gradient_type[which(r2$gradient_type == "icehouse")] <- "Icehouse"
r2$gradient_type[which(r2$gradient_type == "extreme_icehouse")] <- "Extreme Icehouse"
# Set factor levels
r2$gradient_type <- factor(x = r2$gradient_type, levels = c("Extreme Icehouse",
                                                            "Icehouse",
                                                            "Greenhouse",
                                                            "Extreme Greenhouse"))
# Rename sample size
r2$sample_size[which(r2$sample_size == "eocene")] <- "34 (Eocene samples)"
r2$sample_size <- factor(x = r2$sample_size, levels = c(5, 10, 20, "34 (Eocene samples)"))
# Round values
r2$median <- round(r2$median, 3)
r2$p_025 <- round(r2$p_025, 3)
r2$p_975 <- round(r2$p_975, 3)
r2$value <- paste0("*R*^2^ = ", r2$median, "  \n(95% CI: ", r2$p_025, "--", r2$p_975, ")") 
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
  geom_richtext(data = r2, aes(x = 90, y = 39.5, label = value), hjust = 1, fontface = "bold", size = 2.75,
                fill = NA, colour = NA, text.colour = "black") +
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
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
        panel.grid.minor = element_blank())
# Save plot -------------------------------------------------------------
ggsave("./figures/fig_2.png",
       units = "mm", width = 200, height = 200, dpi = 600)
