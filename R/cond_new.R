cond_check <- function(beta, phi, d, trunc) {

  beta_n <- -beta

  c2 <- abs(polyroot(c(1, beta_n))) - 1 # check that invertible

  coef_inf <- ar_infty(ar = phi, ma = beta_n, d = d, max_i = trunc)[-1]

  j <- 1:trunc

  fs <- (j - 1 - d) / j

  if (beta >= 0) {
    r1 <- coef_inf[[1]]
    r2 <- fs[[2]] - phi

    c1_1 <- min(r1, r2)    # must be >= 0

    # or

    idx <- fs[2:(trunc - 1)] < phi & phi <= fs[3:trunc]
    c1_2 <- min(coef_inf[2:(trunc - 1)][idx])  # must be >= 0

  } else if (beta < 0) {

    r1 <- coef_inf[[1]]
    r2 <- coef_inf[[2]]
    r3 <- fs[[2]] * (beta + fs[[3]]) / (beta + fs[[2]]) - phi

    c2_1 <- min(r1, r2, r3)   # must be >= 0

    # or

    idx <- fs[2:(trunc - 2)] * (beta * fs[3:(trunc - 1)]) / (beta + fs[2:(trunc - 2)]) < phi &
      phi <= fs[3:(trunc - 1)] * (beta + fs[4:trunc]) / (beta + fs[3:(trunc - 1)])

    psi_km2 <- psi[2:(trunc - 2)][idx]
    psi_km1 <- psi[3:(trunc - 1)][idx]
    c2_2 <- min(psi_km1, psi_km2)    # must be >= 0

  }

}
