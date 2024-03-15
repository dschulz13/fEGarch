#'Fitting Function
#'
#'@export
#'

egarchfit <- function(rt, spec) {
  dist <- spec@distribution
  orders <- spec@garch_order
  p <- orders[[1]]
  q <- orders[[2]]

  Rsolnp::solnp(
    pars = c(0, 1e-4, rep(1e-1, p), rep(1e-1, p), rep(1e-1, q), 1, 1),
    fun = fEGarch:::llh_norm,
    rt = rt, p = p, q = q, dist = dist,
    LB = rep(-Inf, 4 + 2 * p + q),
    UB = rep(Inf, 4 + 2 * p + q),
    control = list(trace = 0)
  )

}
