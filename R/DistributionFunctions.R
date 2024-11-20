## Absolute moments normal distr.
absMomentNorm <- function(order) {

  force(order)

  FUN <- function(x) {
    x^order * pdf_norm_v1(x)
  }

  2 * integrate(FUN, lower = 0, upper = Inf,
                subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                stop.on.error = FALSE)[[1]]
}

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

## Absolute moments of GEDs.

absMomentGED <- function(order, shape) {

  force(order)

  FUN <- function(x, shape) {
    x^order * pdf_ged_v1(x, shape)
  }

  2 * integrate(FUN, lower = 0, upper = Inf, shape = shape,
                subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                stop.on.error = FALSE)[[1]]

}

## Central moments of skewed standard normal
central_mom_snorm <- function(skew) {
  M1 <- absMomentNorm(order = 1)
  M2 <- absMomentNorm(order = 2)
  Ez <- M1 * (skew - 1 / skew)
  skew2 <- skew^2
  M1_2 <- M1^2
  Varz <- (M2 - M1_2) * (skew^2 + 1 / skew2) + 2 * M1_2 - M2
  c("E" = Ez, "Var" = Varz)
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

## Central moments of skewed standardized GED
central_mom_sged <- function(skew, shape) {
  M1 <- absMomentGED(order = 1, shape = shape)
  M2 <- absMomentGED(order = 2, shape = shape)
  Ez <- M1 * (skew - 1 / skew)
  skew2 <- skew^2
  M1_2 <- M1^2
  Varz <- (M2 - M1_2) * (skew^2 + 1 / skew2) + 2 * M1_2 - M2
  c("E" = Ez, "Var" = Varz)
}

# pdf of standardized skewed normal
pdf_skew_norm_s <- function(x, skew) {
  moms <- central_mom_snorm(skew)
  mom1 <- moms[[1]]
  mom2 <- sqrt(moms[[2]])
  fac <- mom2
  fac * pdf_skew_snorm(x * fac + mom1, skew = skew)
}

# pdf of standardized skewed normal with adjustments for mean and cond. SD
pdf_skew_norm_final <- function(rt, mu, sigt, skew) {
  pdf_skew_norm_s((rt - mu) / sigt, skew = skew) / sigt
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

# pdf of standardized skewed GED
pdf_skew_sged_s <- function(x, shape, skew) {
  moms <- central_mom_sged(skew, shape)
  mom1 <- moms[[1]]
  mom2 <- sqrt(moms[[2]])
  fac <- mom2
  fac * pdf_skew_sged(x * fac + mom1, shape = shape, skew = skew)
}

# pdf of standardized skewed GED with adj. for mean and cond. SD
pdf_skew_sged_final <- function(rt, mu, sigt, shape, skew) {
  pdf_skew_sged_s((rt - mu) / sigt, shape = shape, skew = skew) / sigt
}
