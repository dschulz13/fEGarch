### PDFs of different error distributions

### Standard normal distribution

pdf_norm <- function(x, ...) {
  dnorm(x)
}

### Standardized t-distribution

pdf_std <- function(x, df, ...) {
  sd <- sqrt(df / (df - 2))
  dt(x * sd, df = df) * sd
}


# Skew-normal distribution

pdf_snorm <- function(x, skew_par, ...) {

  delta <- skew_par / sqrt(1 + skew_par^2)
  scale <- 1 / sqrt(1 - (2 * delta^2) / pi)
  location <- -scale * delta * sqrt(2 / pi)

  skew_int <- function(t) {
    (1 / sqrt(2 * pi)) * exp(-t^2 / 2)
  }

  C <- integrate(skew_int, lower = -Inf, upper = skew_par * ((x - location) / scale))[[1]]

  (2 / (scale * sqrt(2 * pi))) * exp(-(x - location)^2 / (2 * scale^2)) * C


}

# Skew-t distribution

E_zt_ <- function(df, skew_par) {

}

# GED

pdf_ged <- function(x, shape, ...) {
  scale <- sqrt(1 / (gamma(3 / shape) / gamma(1 / shape)))
  shape / (2 * scale * gamma(1 / shape)) * exp(-(abs(x) / scale)^shape)
}

# Skew-GED distribution


### Obtain E(|z_t|^p), where z_t are the innovations with mean zero and variance
### one and where p is a power term

E_zt_f <- function(dist, power, ...) {
  f <- switch(dist,
    "norm" = pdf_norm,
    "std" = pdf_std,
    "ged" = pdf_ged,
    "snorm" = pdf_snorm
  )
  int_f <- function(x, skew_par, df, shape, power) {
    abs(x)^power * f(x, df = df, skew_par = skew_par, shape = shape)
  }
  integrate(int_f, -Inf, Inf, skew_par = skew_par, df = df, shape = shape, power = power)[[1]]
}

