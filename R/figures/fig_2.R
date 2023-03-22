# Load libraries --------------------------------------------------------
library(ggplot2)
library(ggthemes)
# Read data -------------------------------------------------------------
df <- readRDS("./results/simulation/simulated_gradients_temperature_estimates.rds")
# Categorise gradients
df$gradient_type[which(df$gradient_type == "flat")] <- "Extreme Greenhouse"
df$gradient_type[which(df$gradient_type == "quadratic")] <- "Greenhouse"
df$gradient_type[which(df$gradient_type == "linear")] <- "Icehouse"
df$gradient_type[which(df$gradient_type == "steep")] <- "Extreme Icehouse"
# Set factor levels
df$gradient_type <- factor(x = df$gradient_type, levels = c("Extreme Icehouse",
                                                            "Icehouse",
                                                            "Greenhouse",
                                                            "Extreme Greenhouse"))
# Plot data -------------------------------------------------------------
col <- "#02818a"
ggplot(data = df, aes(x = lat, y = median)) +
  geom_ribbon(aes(x = lat, ymin = p_100, ymax = p_900),
              fill = col, alpha = 0.25) +
  geom_ribbon(aes(x = lat, ymin = p_005, ymax = p_995),
              fill = col, alpha = 0.25) +
  geom_ribbon(aes(x = lat, ymin = p_025, ymax = p_975),
              fill = col, alpha = 0.25) +
  geom_line(aes(x = lat, y = known), colour = "black", size = 0.75, alpha = 1) +
  geom_line(colour = col, size = 0.75, alpha = 1) +
  facet_grid(gradient_type~sample_size) +
  labs(y = "Sea surface temperature (\u00B0C)", x = "Abs. latitude (\u00B0)") +
  theme_bw()
# Save plot -------------------------------------------------------------
ggsave("./figures/fig_2.png",
       units = "mm", width = 200, height = 200, dpi = 600)
