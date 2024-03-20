### PDFs of different error distributions

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

# Skew-GED distribution (standardized)

pdf_sged <- function(x, df, shape, skew_par) {
  scale <- sqrt(1 / (gamma(3 / shape) / gamma(1 / shape)))
  shape / (2 * scale * gamma(1 / shape)) * exp(-(abs(x) / scale)^shape)
}

### Obtain E(|z_t|^p), where z_t are the innovations with mean zero and variance
### one and where p is a power term

E_zt_f <- function(dist, power, df, shape, skew_par) {
  g <- switch(dist,
    "norm" = pdfnormcpp,
    "std" = pdfstdcpp,
    "ged" = pdfgedcpp
  )

  # if (dist %in% c("snorm", "sstd", "sged")) {
  #   f <- Vectorize(g, vectorize.args = c("x"))
  # } else {
    func <- g
  #}

  if (dist == "norm") {
    int_f <- function(x, power, func) {
      abs(x)^power * func(x = x)
    }
    out <- integrate(int_f, -Inf, Inf, power = power, func = func, subdivisions = 1000L, abs.tol = 1e-20)[[1]]
  } else if (dist == "std") {
    int_f <- function(x, df, power, func) {
      abs(x)^power * func(x = x, df = df)
    }
    out <- integrate(int_f, -Inf, Inf, df = df, power = power, func = func, subdivisions = 1000L, abs.tol = 1e-20)[[1]]
  } else if (dist == "ged") {
    int_f <- function(x, shape, power, func) {
      abs(x)^power * func(x = x, shape = shape)
    }
    out <- integrate(int_f, -Inf, Inf, shape = shape, power = power, func = func, subdivisions = 1000L, abs.tol = 1e-20)[[1]]
  }

  out
}
