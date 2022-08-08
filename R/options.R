#define options

dat1 <- readRDS("data/processed/Hollis_processed_2022_07_19.rds")

temp_dist <- data.frame(name = c("Avicennia", "Avicennia-Rhizophoraceae", "Reefs"),
                        distribution = c("normal", "normal", "normal"),
                        mean = c(mean(c(15.6,22.5)), mean(c(20.7,29.5)), 27.6),
                        sd = c((22.5-15.6)/4, c(29.5-20.7)/4, (29.5-21)/4),
                        shape = rep(NA,3)
)

22.5-15.6

x <- seq(16,36,0.1)
mu = 27.6
mint <- 21
maxt <- 29.5
dens <- dnorm(x,mu,(maxt-mint)/4)
plot(x,dens,type = "l", lwd = 2)
abline(v = c(mu,mint,maxt), lty = c(1,2,2))

densm <- dnorm(x,mean(c(20.7,29.5)), c(29.5-20.7)/4)
points(x,densm,type = "l", col = "red")

densa <- dnorm(x,mean(c(15.6,22.5)), c(22.5-15.6)/4)
points(x,0.63*densa,type = "l", col = "blue")
