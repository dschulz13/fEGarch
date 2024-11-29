#'Value-at-Risk Calculation
#'
#'Compute the value-at-risk (VaR) as a quantile of the
#'standardized innovation distribution for all six
#'distributions. Makes use of a numerical approach,
#'where the uniroot() function looks for the input
#'quantile of the distribution's corresponding CDF
#'R function that matches the "level" argument.
#'This way, internally, there are inverse CDF functions
#'defined, i.e. functions that return the quantile for
#'an input cumulative probability.
#'
baseVaR_calc <- function(level = 0.95, cond_dist, ...) {
  dFUN <- lookup_table$fun1[[cond_dist]]
  inv_cdf_fun(p = level, dfun = dFUN, ...)
}

#'ES function does not work properly because uniroot fails at some point
#'
#'
baseES_calc <- function(level = 0.95, cond_dist, ...) {

  integ <- integrate(baseVaR_calc, lower = level, upper = Inf, cond_dist = cond_dist, ...)[[1]]

  integ / (1 - level)

}
