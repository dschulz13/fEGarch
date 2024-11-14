cdf_fun <- function(x, dfun, ...) {
  vapply(
    x,
    FUN = function(x, dfun, ...) {
      integrate(f = dfun, lower = -Inf, upper = x,
                ..., subdivisions = 1000L, abs.tol = 1e-15)[[1]]
    },
    FUN.VALUE = numeric(1), dfun = dfun, ...
  )
}

inv_cdf_fun_helper <- function(x, p, dfun, ...) {
  cdf_fun(x, dfun = dfun, ...) - p
}

inv_cdf_fun <- function(p, dfun, ...) {
  vapply(
    p,
    FUN = function(p, dfun, ...) {
       uniroot(inv_cdf_fun_helper, interval = c(-5, 5), dfun = dfun, p = p,
               ..., extendInt = "yes")[[1]]
    },
    FUN.VALUE = numeric(1), dfun = dfun, ...
  )

}

#'Sampling Functions for Innovations
#'
#'Draw random samples from a normal distribution, \eqn{t}-distribution,
#'generalized error distribution, or their skewed variants (all standardized
#'to have mean zero and variance one).
#'
#'@param n the number of observations to draw.
#'@param df the degrees of freedom for a (skewed) \eqn{t}-distribution.
#'@param shape the shape parameter for a (skewed) generalized error
#'distribution.
#'@param skew the skewness parameter in the skewed distributions.
#'
#'@details
#'Draw random samples from a normal distribution, \eqn{t}-distribution,
#'generalized error distribution, or their skewed variants (all standardized
#'to have mean zero and variance one).
#'
#'@return
#'These functions return a numeric vector of length \code{n}.
#'
#'@name sim_functions
#'
#'@export
#'
#'@examples
#'rnorm_s(10)
#'
#'rsstd_s(10, df = 7, skew = 0.9)

rnorm_s <- function(n) {
  rnorm(n = n, mean = 0, sd = 1)
}

#'@export
#'@rdname sim_functions

rstd_s <- function(n, df = 10000) {
  sdev <- sqrt(df / (df - 2))
  rt(n = n, df = df) / sdev
}
#'@export
#'@rdname sim_functions

rged_s <- function(n, shape = 2) {
  p <- runif(n)
  inv_cdf_fun(p = p, dfun = fEGarch:::pdf_ged_v1, shape = shape)
}

#'@export
#'@rdname sim_functions
rsnorm_s <- function(n, skew = 1) {
  p <- runif(n)
  inv_cdf_fun(p = p, dfun = fEGarch:::pdf_skew_norm_s, skew = skew)
}

#'@export
#'@rdname sim_functions
rsstd_s <- function(n, df = 10000, skew = 1) {
  p <- runif(n)
  inv_cdf_fun(p = p, dfun = fEGarch:::pdf_skew_sstd_s, df = df, skew = skew)
}

#'@export
#'@rdname sim_functions
rsged_s <- function(n, shape = 2, skew = 1) {
  p <- runif(n)
  inv_cdf_fun(p = p, dfun = fEGarch:::pdf_skew_sged_s, shape = shape, skew = skew)
}
