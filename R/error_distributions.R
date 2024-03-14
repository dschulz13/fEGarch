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


# Skew-normal distribution (standardized)

pdf_snorm <- function(x, df, shape, skew_par) {

  delta <- skew_par / sqrt(1 + skew_par^2)
  scale <- 1 / sqrt(1 - (2 * delta^2) / pi)
  location <- -scale * delta * sqrt(2 / pi)

  skew_int <- function(t) {
    (1 / sqrt(2 * pi)) * exp(-t^2 / 2)
  }

  C_x <- integrate(skew_int, lower = -Inf, upper = skew_par * ((x - location) / scale), subdivisions = 1000L, abs.tol = 1e-20)[[1]]

  (2 / (scale * sqrt(2 * pi))) * exp(-(x - location)^2 / (2 * scale^2)) * C_x
}

# Skew-t distribution (standardized)

pdf_sstd <- function(x, df, shape, skew_par) {

}

# GED (standardized)

pdf_ged <- function(x, df, shape, skew_par) {
  scale <- sqrt(1 / (gamma(3 / shape) / gamma(1 / shape)))
  shape / (2 * scale * gamma(1 / shape)) * exp(-(abs(x) / scale)^shape)
}

# Skew-GED distribution (standardized)

pdf_sged <- function(x, df, shape, skew_par) {
  scale <- sqrt(1 / (gamma(3 / shape) / gamma(1 / shape)))
  shape / (2 * scale * gamma(1 / shape)) * exp(-(abs(x) / scale)^shape)
}

### Obtain E(|z_t|^p), where z_t are the innovations with mean zero and variance
### one and where p is a power term

E_zt_f <- function(dist, power, df = 10, shape = 0, skew_par = 0) {
  g <- switch(dist,
    "norm" = pdf_norm,
    "std" = pdf_std,
    "ged" = pdf_ged,
    "snorm" = pdf_snorm,
    "sstd" = pdf_sstd,
    "sged" = pdf_sged
  )

  if (dist %in% c("snorm", "sstd", "sged")) {
    f <- Vectorize(g, vectorize.args = c("x"))
  } else {
    f <- g
  }

  int_f <- function(x, power, df, shape, skew_par) {
    abs(x)^power * f(x = x, df = df, shape = shape, skew_par = skew_par)
  }
  integrate(int_f, -Inf, Inf, power = power, df = df, shape = shape, skew_par = skew_par, subdivisions = 1000L, abs.tol = 1e-20)[[1]]
}

E_zt_f("snorm", power = 1, skew_par = 30)
