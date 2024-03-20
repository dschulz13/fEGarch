#'Fitting Function
#'
#'@export
#'
#'@importFrom Rcpp evalCpp
#'@useDynLib fEGarch
#'

egarchfit <- function(rt, spec) {
  dist <- spec@distribution
  orders <- spec@garch_order
  p <- orders[[1]]
  q <- orders[[2]]

  llikhood <- llikhood
  stopifnot("Test" = is.function(llikhood))

  stats::nlminb(
    start = c(0, 1e-4, rep(1e-1, p), rep(1e-1, p), rep(1e-1, q), 1, 1),
    objective = llikhood,
    rt = rt, p = p, q = q, dist = dist, include_mean = TRUE,
    lower = c(-Inf, 0, rep(-Inf, 2 * p + q), 0, 0),
    upper = c(Inf, Inf, rep(Inf, 2 * p + q), 2, 2),
    control = list(trace = 0)


  )

}
