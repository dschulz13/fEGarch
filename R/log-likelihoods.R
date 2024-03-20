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

llhood <- function(pars, rt, p, q, include_mean, dist) {

  if (include_mean) {
    i_p1 <- 2 + 2 * p + q + 1
  } else {
    i_p1 <- 2 + 2 * p + q
  }
  i_p2 <- i_p1 + 1
  power1 <- pars[[i_p1]]
  power2 <- pars[[i_p2]]

  df <- 0
  shape <- 0
  skew <- 0

  if (dist == "std") {
    df <- pars[[i_p2 + 1]]
  } else if (dist == "ged") {
    shape <- pars[[i_p2 + 1]]
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
  llh(pars, rt, p, q, include_mean, E1, E2)

}


