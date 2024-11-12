#'Fitting Function
#'
#'@export
#'
#'@importFrom Rcpp evalCpp
#'@useDynLib fEGarch
#'


fEGarchFit <- function(rt,
    model = c("mloggarch", "megarch", "filoggarch", "fiegarch",
              "mafiloggarch"),
    order = c(1L, 1L),
    include_mean = TRUE,
    cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
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
