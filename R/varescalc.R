#Value-at-Risk Calculation
#
#Compute the value-at-risk (VaR) as a quantile of the
#standardized innovation distribution for all six
#distributions. Makes use of a numerical approach,
#where the uniroot() function looks for the input
#quantile of the distribution's corresponding CDF
#R function that matches the one minus the "level" argument.
#The result is thus with negative sign.
#This way, internally, there are inverse CDF functions
#defined, i.e. functions that return the quantile for
#an input cumulative probability.
#
#
#
baseVaR_calc <- function(level = 0.95, cond_dist, shape, skew) {
  dFUN <- lookup_table$fun1[[cond_dist]]
  vapply(
    X = 1 - level,
    FUN = function(.x, dfun, shape, skew) {
      inv_cdf_fun(p = .x, dfun = dfun, shape = shape, skew = skew)
    }, FUN.VALUE = numeric(1),
    dfun = dFUN, shape = shape, skew = skew
  )
}

#ES function does not work properly because uniroot fails at some point
#
#
baseES_calc <- function(level = 0.95, cond_dist, shape, skew) {

  q <- baseVaR_calc(level = level, cond_dist = cond_dist, shape = shape, skew = skew)
  dfun <- lookup_table[["fun1"]][[cond_dist]]
  int_fun <- function(x, dfun, shape, skew) {
    x * dfun(x = x, shape = shape, skew = skew)
  }

  integ <- vapply(
    X = q,
    FUN = function(.x, dfun, shape, skew, int_fun) {
      integrate(int_fun, lower = -Inf, upper = .x, dfun = dfun, shape = shape, skew = skew,
              subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6)[[1]]
    }, FUN.VALUE = numeric(1),
    shape = shape, skew = skew, dfun = dfun, int_fun = int_fun
  )

  integ / (1 - level)

}

#'VaR and ES Computation for Standardized Distributions
#'
#'Compute the value at risk (VaR) and the expected shortfall (ES) numerically
#'for the standardized distributions available in this package. These quantiles
#'can then be used to obtain the conditional VaR and ES following GARCH-type models.
#'
#'@param level a numeric vector with the confidence level(s) to calculate
#'the VaR or ES for; the quantiles are the VaR or ES
#'computed at one minus the input level; the result is thus with negative sign
#'for common \code{level} input such as \code{0.975} or \code{0.99}.
#'@param dist a single character value (or the default vector) that
#'specifies the distribution to consider (all distributions are
#'considered to be standardized with mean zero and variance one).
#'@param ... further arguments to consider for the distributions;
#'for \code{"std"} or \code{"sstd"}, specify the degrees of freedom
#'\code{df}, for \code{"ged"} or \code{"sged"}, give the shape
#'parameter \code{shape}, and for \code{"ald"} or \code{"sald"},
#'use the additional argument \code{P}; moreover, for the skewed
#'distributions \code{"snorm"}, \code{"sstd"}, \code{"sged"} and
#'\code{"sald"}, the skewness parameter \code{skew} must
#'be provided as well.
#'
#'@details
#'The VaR is found numerically using numerical root finding
#'via \code{uniroot} (of the \code{stats} package), whereas the ES is obtained through
#'numerical integration, where firstly the VaR at the corresponding
#'confidence level is computed using \code{VaR_calc}, and where
#'subsequently \code{integrate} (of the \code{stats} package) is
#'used at the tail of the distribution.
#'
#'@export
#'
#'@rdname VaR-ES-Calculation
#'
#'@return
#'Returns a numeric vector of the same length as the input
#'argument \code{level}.
#'
#'@examples
#'
#'# 99-percent VaR, normal distribution
#'VaR_calc(level = 0.99, dist = "norm")
#'
#'# 99-percent VaR, t-distribution with df = 10
#'VaR_calc(level = 0.99, dist = "std", df = 10)
#'
#'# 97.5-percent ES, normal distribution
#'ES_calc(level = 0.975, dist = "norm")
#'
#'# 97.5-percent ES, t-distribution with df = 10
#'ES_calc(level = 0.975, dist = "std", df = 10)
#'
#'

VaR_calc <- function(level = 0.99,
                     dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                     ...) {
  dist <- match.arg(dist)
  l <- list(...)
  shape <- if (dist %in% c("std", "sstd")) {
    l[["df"]]
  } else if (dist %in% c("ald", "sald")) {
    l[["P"]]
  } else {
    l[["shape"]]
  }
  if (!(dist %in% c("norm", "snorm"))) {
    stopifnot('Provide additional arguments "df", "shape" or "P" for the corresponding distributions "std", "ged", "ald", "sstd", "sged" and "sald".' = !is.null(shape))
  }
  skew <- l[["skew"]]
  if (dist %in% c("snorm", "sstd", "sged", "sald")) {
    stopifnot('For skewed distributions, please provide argument "skew".' = !is.null(skew))
  }
  baseVaR_calc(level = level, cond_dist = dist, shape = shape, skew = skew)
}

#'@export
#'@rdname VaR-ES-Calculation

ES_calc <- function(level = 0.975,
                     dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                     ...) {
  dist <- match.arg(dist)
  l <- list(...)
  shape <- if (dist %in% c("std", "sstd")) {
    l[["df"]]
  } else if (dist %in% c("ald", "sald")) {
    l[["P"]]
  } else {
    l[["shape"]]
  }
  if (!(dist %in% c("norm", "snorm"))) {
    stopifnot('Provide additional arguments "df", "shape" or "P" for the corresponding distributions "std", "ged", "ald", "sstd", "sged" and "sald".' = !is.null(shape))
  }
  skew <- l[["skew"]]
  if (dist %in% c("snorm", "sstd", "sged", "sald")) {
    stopifnot('For skewed distributions, please provide argument "skew".' = !is.null(skew))
  }

  baseES_calc(level = level, cond_dist = dist, shape = shape, skew = skew)
}
