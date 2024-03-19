g_power1 <- function(innov, kappa, power1, E1) {
  kappa * (sign(innov) * abs(innov)^power1 - E1) / power1
}

g_power2 <- function(innov, gamma, power2, E2) {
  gamma * (abs(innov)^power2 - E2) / power2
}

ged_dens <- function(x, shape = 2) {
  scale <- sqrt(gamma(1 / shape) / gamma(3 / shape))
  log(x^2) * shape / (2 * scale * gamma(1 / shape)) * exp(-(abs(x) / scale)^shape)
}

# Expectation of ln(zt^2) where zt are the innovations with mean zero and variance one

Elnz2 <- function(dist, df, shape) {
  switch(dist,
    "norm" = digamma(0.5) - log(0.5),
    "std" = log(df - 2) + digamma(0.5) - digamma(df / 2),
    "ged" = integrate(ged_dens, -Inf, Inf, shape = shape, abs.tol = 1e-20)[[1]]
  )
}


llhood <- function(parameters, rt, p, q, include_mean, dist) {

  if (include_mean) {
    i_p1 <- 2 + 2 * p + q + 1
  } else {
    i_p1 <- 2 + 2 * p + q
  }
  i_p2 <- i_p1 + 1
  power1 <- parameters[[i_p1]]
  power2 <- parameters[[i_p2]]

  df <- 0
  shape <- 0
  skew <- 0

  if (dist == "std") {
    df <- parameters[[i_p2 + 1]]
  } else if (dist == "ged") {
    shape <- parameters[[i_p2 + 1]]
  }

  # Obtain E(|z_t|^p), where z_t are the innovations with mean zero and variance
  # one and where p is a power term
  E1 <- 0     # adjust later on for non-symmetric densities
  E2 <- E_zt_f(dist, power2, df = df, shape = shape, skew_par = skew)

  # Select the log-likelihood function according to the conditional distribution
  llh <- switch(dist,
    "norm" = llhNORM,
    "std" = llhSTD,
    "ged" = llhGED
  )

  # Call to the log-likelihood C++ function
  llh(parameters, rt, p, q, include_mean, E1, E2)

}


# Under standard normal innovations

llh_norm <- function(parameters, rt, p, q, dist) {

  mu <- parameters[[1]]
  omega <- parameters[[2]]
  kappa <- gamma <- beta <- 0    # adjust in future version

  ind <- 3

  if (p > 0) {
    begin_kappa <- ind
    end_kappa <- 3 + p - 1
    begin_gamma <- end_kappa + 1
    end_gamma <- begin_gamma + p - 1
    kappa <- parameters[begin_kappa:end_kappa]
    gamma <- parameters[begin_gamma:end_gamma]
    ind <- end_gamma + 1
  }
  if (q > 0) {
    begin_beta <- ind
    end_beta <- end_gamma + p - 1
    beta <- parameters[begin_beta:end_beta]
    ind <- end_beta + 1
  }
  power1 <- parameters[[ind]]
  power2 <- parameters[[ind + 1]]
  df <- shape <- skew <- 0     # not required for standard normal innovations
  # if (dist == "std") {
  #   df <- parameters[[end_beta + 5]]
  # }
  # if (dist == "ged") {
  #   shape <- parameters[[end_beta + 5]]
  # }    # adjust later on

  E1 <- 0     # adjust later on for non-symmetric densities
  E2 <- E_zt_f(dist, power2, df = df, shape = shape, skew_par = skew)

  n <- length(rt)

  lsigt  <- zt <- rep(0, n)

  m <- max(p, q)
  lsigt[1:m] <- mean(log(rt^2)) - Elnz2("norm", df = 0)
  zt[1:(m - 1)] <- rt[1:(m - 1)] / exp(lsigt[1:(m - 1)] / 2)


  for (i in (m + 1):n) {
    zt[[i - 1]] <- rt[[i - 1]] / exp(lsigt[[i - 1]] / 2)
    lsigt[[i]] <- omega + sum(g_power1(zt[(i - 1):(i - p)], kappa, power1, E1) + g_power2(zt[(i - 1):(i - p)], gamma, power2, E2)) + beta %*% lsigt[(i - 1):(i - q)]
  }

  sigt2 <- exp(lsigt)

  sum(0.5 * (log(2 * pi) + log(sigt2) + (rt - mu)^2 / sigt2))

}
