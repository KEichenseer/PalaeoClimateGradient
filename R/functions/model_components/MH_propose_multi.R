# create multivariate normal proposals
MH_propose_multi <- function(nprop,coeff,proposal_cov) {
  mvnfast::rmvn(n = nprop, mu = 0.95*c(coeff[1:3],log(coeff[4])),
                sigma = 2.4/sqrt(4)*proposal_cov)+
    rnorm(n = 4*nprop, mean = 0.05*c(coeff[1:3],log(coeff[4])),
          sd = 0.001)
}