source("R/options.R")
source("R/functions/model_components/dsnorm.R")
source("R/functions/model_components/dtnorm.R")

# Priors
n <- 1000
xval = list(A = seq(-6, 36, length.out = n),
            K = seq(-5, 65, length.out = n),
            M = seq(9, 80, length.out = n),
            B = seq(0, 0.4, length.out = n)
)
A <- priors$f1(xval[["A"]], log = FALSE)
A <- A / max(A)
K <- priors$f2(xval[["K"]], xval[["K"]][which.max(A)], log = FALSE)
K <- K / max(K)
M <- priors$f3(xval[["M"]], log = FALSE)
M <- M / max(M)
B <- priors$f4(xval[["B"]], log = FALSE)
B <- B / max(B)

# Set up dataframe
df <- data.frame(x = c(xval[["A"]], xval[["K"]], xval[["M"]], xval[["B"]]),
                 y = c(A, K, M, B),
                 prior = rep(c("A", "K", "M", "B"), each = n))

ggplot(data = df, aes(y = y)) +
  geom_density(fill="royalblue3") +
  facet_wrap(~prior, scales = "free")
