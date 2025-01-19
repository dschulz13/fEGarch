## Absolute moments t-distr.
absMomentStd <- function(order, df) {

  force(order)

  FUN <- function(x, df) {
    x^order * pdf_std_v1(x, df)
  }

  2 * integrate(FUN, lower = 0, upper = Inf, df = df,
                subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                stop.on.error = FALSE)[[1]]

}

## Central moments of skewed standardized t
central_mom_sstd <- function(skew, df) {
  M1 <- absMomentStd(order = 1, df = df)
  M2 <- absMomentStd(order = 2, df = df)
  Ez <- M1 * (skew - 1 / skew)
  skew2 <- skew^2
  M1_2 <- M1^2
  Varz <- (M2 - M1_2) * (skew^2 + 1 / skew2) + 2 * M1_2 - M2
  c("E" = Ez, "Var" = Varz)
}

# pdf of standardized skewed t-distr.
pdf_skew_sstd_s <- function(x, df, skew) {
  moms <- central_mom_sstd(skew, df)
  mom1 <- moms[[1]]
  mom2 <- sqrt(moms[[2]])
  fac <- mom2
  fac * pdf_skew_sstd(x * fac + mom1, df = df, skew = skew)
}

# pdf of standardized skewed t-distr. with adj. for mean and cond. SD
pdf_skew_sstd_final <- function(rt, mu, sigt, df, skew) {
  pdf_skew_sstd_s((rt - mu) / sigt, df = df, skew = skew) / sigt
}
