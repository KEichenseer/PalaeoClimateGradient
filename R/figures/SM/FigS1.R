# Supplementary figure 1: Traceplots of Markov chains of the EECO model
eeco_mod <- readRDS("./results/eeco/eeco_climate_model_output.rds")
source("R/functions/auxiliary_functions.R")

# 1. Traceplot of A, K, M and log(B)
png("figures/SM/S1_traceplot_eeco.png", width = 3.5, height = 4.9, units = "in", res = 400)
plot_chains(eeco_mod, select = list("params", 1:5), xlabel = expression("iteration"), ylabel = c("A", "K - A", "M", "B", expression(sigma)),
            xaxt = "n", xaxsat = c(1,c(1,2,3,4,5,6)*10^5), xaxslab = c(1,"100,000","200,000","300,000",
                                                                       "400,000","500,000","600,000"))
dev.off()

# 2. Traceplot of some of the location mu
png("figures/SM/S1_traceplot_eeco_muj.png", width = 3.5, height = 4.9, units = "in", res = 400)
plot_chains(eeco_mod, select = list("yestimate", c(1,8,10,25,30)), xlabel = expression("iteration"), ylabel = unlist(sapply(c(1,8,10,25,30), function(x) 
  bquote(mu[.(x)]))),
            xaxt = "n", xaxsat = c(1,c(1,2,3,4,5,6)*10^5), xaxslab = c(1,"100,000","200,000","300,000",
                                                                       "400,000","500,000","600,000"))
dev.off()

## comment: large values of sdyest occur when yestimate is far from the empirical mean
# 3. Traceplot of some of the location sd
png("figures/SM/S1_traceplot_eeco_sigmaj.png", width = 3.5, height = 4.9, units = "in", res = 400)
plot_chains(eeco_mod, select = list("sdyest", c(1,8,10,14,22)), xlabel = expression("iteration"), ylabel = unlist(sapply(c(1,8,10,14,22), function(x) 
  bquote(sigma[.(x)]))),
  xaxt = "n", xaxsat = c(1,c(1,2,3,4,5,6)*10^5), xaxslab = c(1,"100,000","200,000","300,000",
                                                             "400,000","500,000","600,000"))
dev.off()