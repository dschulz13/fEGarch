## Absolute moments of Sargan distr.

absMomentSM <- function(order, P) {

  force(order)

  FUN <- function(x, P) {
    x^order * c(pdf_sm_v1(x, P))
  }

  2 * integrate(FUN, lower = 0, upper = Inf, P = P,
                subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                stop.on.error = FALSE)[[1]]

}

## Central moments of skewed standardized Sargan distr.
central_mom_ssm <- function(skew, P) {
  M1 <- absMomentSM(order = 1, P = P)
  M2 <- absMomentSM(order = 2, P = P)
  Ez <- M1 * (skew - 1 / skew)
  skew2 <- skew^2
  M1_2 <- M1^2
  Varz <- (M2 - M1_2) * (skew^2 + 1 / skew2) + 2 * M1_2 - M2
  c("E" = Ez, "Var" = Varz)
}

# pdf of standardized skewed Sargan distr.
pdf_skew_ssm_s <- function(x, P, skew) {
  moms <- central_mom_ssm(skew, P)
  mom1 <- moms[[1]]
  mom2 <- sqrt(moms[[2]])
  fac <- mom2
  fac * c(pdf_skew_ssm(x * fac + mom1, P = P, skew = skew))
}

# pdf of standardized skewed Sargan distr. with adj. for mean and cond. SD
pdf_skew_ssm_final <- function(rt, mu, sigt, shape, skew) {
  pdf_skew_ssm_s((rt - mu) / sigt, P = P, skew = skew) / sigt
}
