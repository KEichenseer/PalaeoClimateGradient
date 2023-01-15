# Supplementary figure 1: Traceplots of Markov chains of the EECO model

# 1. Traceplot of A, K, M and log(B)
eeco_mod <- readRDS("./results/eeco/eeco_climate_model_output.rds")
source("R/functions/auxiliary_functions.R")
png("figures/SM/S1_traceplot_eeco.png", width = 4.5, height = 5, units = "in", res = 400)
plot_chains(eeco_mod, xlabel = expression("iteration (thinned to every 100"^"th"*")"))
dev.off()

# 2. Traceplot of some of the location mu !!! TODO !!!

# plot_chains_mu(eeco_mod, xlabel = expression("iteration (thinned to every 100"^"th"*")"))
