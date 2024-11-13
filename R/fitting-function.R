#'Fitting Function for Models of the Broader EGARCH Family
#'
#'Use maximum-likelihood estimation to fit a model from the
#'broader EGARCH family to some observed time series.
#'
#'@param rt the input time series to fit the model to ordered from
#'past to present.
#'@param spec an S4 object of class \code{"egarch_spec"} that can
#'be created via the function \code{egarch_spec()}.
#'@param start_pars a vector with starting parameters for the
#'optimization; must be of the same length as the output vector of
#'parameters; the default NULL uses internally saved default sets
#'of starting values.
#'@param LB a vector with lower boundaries for parameters; must be
#'of the same length as the output vector of
#'parameters; the default NULL uses internally saved default sets
#'of lower boundary values.
#'@param UB a vector with upper boundaries for parameters; must be
#'of the same length as the output vector of
#'parameters; the default NULL uses internally saved default sets
#'of upper boundary values.
#'
#'@export
#'
#'@importFrom Rcpp evalCpp
#'@useDynLib fEGarch
#'
#'@details
#'To be filled soon.
#'


fEGarchFit <- function(rt,
    spec = egarch_spec(),
    start_pars = NULL, LB = NULL, UB = NULL) {

  model <- match.arg(model)
  cond_dist <- match.arg(cond_dist)

  fitFUN <- fitFUNS[[model]]

  fitFUN(
    rt = rt,
    order = order,
    include_mean = include_mean,
    cond_dist = cond_dist,
    start_pars = start_pars,
    LB = LB,
    UB = UB
  )

}


# fEGarchFit <- function(rt, spec) {
#   dist <- spec@distribution
#   orders <- spec@garch_order
#   p <- orders[[1]]
#   q <- orders[[2]]
#
#   llikhood <- llikhood
#   stopifnot("Test" = is.function(llikhood))
#
#   stats::nlminb(
#     start = c(0, 0, rep(0, p), rep(0, p), rep(0, q), 1, 1),
#     objective = llikhood,
#     rt = rt, p = p, q = q, dist = dist, include_mean = TRUE,
#     lower = c(-Inf, -Inf, rep(-Inf, 2 * p + q), 0, 0),
#     upper = c(Inf, Inf, rep(Inf, 2 * p + q), 2, 2),
#     control = list(trace = 0)
#
#
#   )
#
# }
