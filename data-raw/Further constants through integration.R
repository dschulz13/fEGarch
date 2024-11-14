FUN1 <- function(x, expo, dfun, ...) {
  (sign(x) * (abs(x) + 1)^expo / expo) * dfun(x, ...)
}

dfun <- fEGarch:::pdf_skew_norm_s
expo <- 0.5
skew <- 1.3

E1 <- integrate(FUN1, -Inf, Inf, abs.tol = 1e-15, expo = expo,
                dfun = dfun, skew = skew)[[1]]

FUN <- function(x, expo, dfun, E1, ...) {
  (sign(x) * (abs(x) + 1)^expo / expo - E1) * dfun(x, ...)
}

integrate(FUN, -Inf, Inf, abs.tol = 1e-15, expo = expo,
          dfun = dfun, E1 = E1, skew = skew)[[1]]
