#'Fitting Function
#'
#'@export
#'

egarchfit <- function(rt, spec) {
  dist <- spec@distribution
  orders <- spec@garch_order
  p <- orders[[1]]
  q <- orders[[2]]

  llhood <- llhood
  stopifnot("Test" = is.function(llhood))

  Rsolnp::solnp(
    pars = c(0, 1e-4, rep(1e-1, p), rep(1e-1, p), rep(1e-1, q), 1, 1),
    fun = llhood,
    LB = c(-Inf, 0, rep(-Inf, 2 * p + q), 0, 0),
    UB = c(Inf, Inf, rep(Inf, 2 * p + q), 2, 2),
    control = list(trace = 0),
    rt = rt, p = p, q = q, dist = dist, include_mean = TRUE

  )

}
