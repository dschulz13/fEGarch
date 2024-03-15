#'EGARCH Family Simulation Function
#'
#'@export
#'

sim_fun <- function(spec, n = 500, n_start = 1000) {
  dist <- "norm"      # possibly change

  mu <- 0
  omega <- -1.8
  kappa <- 0.1
  gamma <- 0.15
  beta <- 0.8
  power1 <- 1.5
  power2 <- 2

  l_k <- length(kappa)
  l_b <- length(beta)
  m <- max(l_k, l_b)

  r <- n + n_start + m
  z_t <- rnorm(r)

  eta_t <- (sign(z_t) * abs(z_t)^power1) / power1
  E_const <- E_zt_f(dist = dist, power = power2)
  u_t <- (abs(z_t)^power2 - E_const) / power2

  lsig2 <- rep(omega / 1 - sum(beta), times = r)

  p <- l_k
  q <- l_b

  for (i in (m + 1):r) {
    lsig2[[i]] <- omega + kappa %*% eta_t[(i - 1):(i - l_k)] +
      gamma %*% u_t[(i - 1):(i - p)] +
      beta %*% lsig2[(i - 1):(i - q)]
  }

  sig <- tail(exp(lsig2 / 2), n)
  zt <- tail(z_t, n)
  rt <- sig * zt + mu

  list(series = rt, sigma = sig, innov = zt)

}
