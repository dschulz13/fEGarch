### PDFs of different error distributions

### Standard normal distribution

pdf_norm <- function(x, df, shape, skew_par) {
  dnorm(x)
}

### Standardized t-distribution

pdf_std <- function(x, df, shape, skew_par) {
  sd <- sqrt(df / (df - 2))
  dt(x * sd, df = df) * sd
}


# Skew-normal distribution

pdf_snorm <- function(x, df, shape, skew_par) {

  delta <- skew_par / sqrt(1 + skew_par^2)
  scale <- 1 / sqrt(1 - (2 * delta^2) / pi)
  location <- -scale * delta * sqrt(2 / pi)

  skew_int <- function(t) {
    (1 / sqrt(2 * pi)) * exp(-t^2 / 2)
  }

  C <- vapply(x,
    FUN = function(.x, skew_par, location, scale, skew_int) {
      integrate(skew_int, lower = -Inf, upper = skew_par * ((.x - location) / scale))[[1]]
    }, skew_par = skew_par, location = location, scale = scale, skew_int = skew_int,
    FUN.VALUE = numeric(1)
  )

  (2 / (scale * sqrt(2 * pi))) * exp(-(x - location)^2 / (2 * scale^2)) * C


}

# Skew-t distribution

pdf_sstd <- function(x, df, shape, skew_par) {

}

# GED

pdf_ged <- function(x, df, shape, skew_par) {
  scale <- sqrt(1 / (gamma(3 / shape) / gamma(1 / shape)))
  shape / (2 * scale * gamma(1 / shape)) * exp(-(abs(x) / scale)^shape)
}

# Skew-GED distribution

pdf_sged <- function(x, df, shape, skew_par) {
  scale <- sqrt(1 / (gamma(3 / shape) / gamma(1 / shape)))
  shape / (2 * scale * gamma(1 / shape)) * exp(-(abs(x) / scale)^shape)
}

### Obtain E(|z_t|^p), where z_t are the innovations with mean zero and variance
### one and where p is a power term

E_zt_f <- function(dist, power, df = 10, shape = 0, skew_par = 0) {
  f <- switch(dist,
    "norm" = pdf_norm,
    "std" = pdf_std,
    "ged" = pdf_ged,
    "snorm" = pdf_snorm,
    "sstd" = pdf_sstd,
    "sged" = pdf_sged
  )

  int_f <- function(x, power, df, shape, skew_par) {
    abs(x)^power * f(x = x, df = df, shape = shape, skew_par = skew_par)
  }
  integrate(int_f, -Inf, Inf, power = power, df = df, shape = shape, skew_par = skew_par)[[1]]
}


