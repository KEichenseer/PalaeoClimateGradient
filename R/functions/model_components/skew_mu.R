# Gibbs sampling of mu with skew normal likelihood and normal prior
### this is for single observations (sample mean), so no mean(x) or mean(y)
skew_mu <- function(x, y, sigma, rho, mu_prior, sigma_prior) {
  n1 = 1
  rnorm(length(x),(n1*sigma_prior^2*(x-sigma*rho*y)+sigma^2*(1-rho^2)*mu_prior)/
          (n1*sigma_prior^2+sigma^2*(1-rho^2)),
        sqrt((sigma_prior^2*sigma^2*(1-rho^2))/(n1*sigma_prior^2+sigma^2*(1-rho^2))) )
}
