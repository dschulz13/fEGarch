cdf_fun <- function(x, dfun, ...) {
  vapply(
    x,
    FUN = function(x, dfun, ...) {
      integrate(f = dfun, lower = -Inf, upper = x,
                ..., subdivisions = 1000L, abs.tol = 1e-15)[[1]]
    },
    FUN.VALUE = numeric(1), dfun = dfun, ...
  )
}

inv_cdf_fun_helper <- function(x, p, dfun, ...) {
  cdf_fun(x, dfun = dfun, ...) - p
}

inv_cdf_fun <- function(p, dfun, ...) {
  vapply(
    p,
    FUN = function(p, dfun, ...) {
       uniroot(inv_cdf_fun_helper, interval = c(-5, 5), dfun = dfun, p = p,
               ..., extendInt = "yes")[[1]]
    },
    FUN.VALUE = numeric(1), dfun = dfun, ...
  )

}

#'Sampling Functions for Innovations
#'
#'Draw random samples from a normal distribution, \eqn{t}-distribution,
#'generalized error distribution, or their skewed variants (all standardized
#'to have mean zero and variance one).
#'
#'@param n the number of observations to draw.
#'@param df the degrees of freedom for a (skewed) \eqn{t}-distribution.
#'@param shape the shape parameter for a (skewed) generalized error
#'distribution.
#'@param skew the skewness parameter in the skewed distributions.
#'
#'@details
#'Draw random samples from a normal distribution, \eqn{t}-distribution,
#'generalized error distribution, or their skewed variants (all standardized
#'to have mean zero and variance one).
#'
#'@return
#'These functions return a numeric vector of length \code{n}.
#'
#'@name sim_functions
#'
#'@export
#'
#'@examples
#'rnorm_s(10)
#'
#'rsstd_s(10, df = 7, skew = 0.9)

rnorm_s <- function(n) {
  rnorm(n = n, mean = 0, sd = 1)
}

#'@export
#'@rdname sim_functions

rstd_s <- function(n, df = 10000) {
  sdev <- sqrt(df / (df - 2))
  rt(n = n, df = df) / sdev
}
#'@export
#'@rdname sim_functions

rged_s <- function(n, shape = 2) {
  p <- runif(n)
  inv_cdf_fun(p = p, dfun = pdf_ged_v1, shape = shape)
}

#'@export
#'@rdname sim_functions
rsnorm_s <- function(n, skew = 1) {
  p <- runif(n)
  inv_cdf_fun(p = p, dfun = pdf_skew_norm_s, skew = skew)
}

#'@export
#'@rdname sim_functions
rsstd_s <- function(n, df = 10000, skew = 1) {
  p <- runif(n)
  inv_cdf_fun(p = p, dfun = pdf_skew_sstd_s, df = df, skew = skew)
}

#'@export
#'@rdname sim_functions
rsged_s <- function(n, shape = 2, skew = 1) {
  p <- runif(n)
  inv_cdf_fun(p = p, dfun = pdf_skew_sged_s, shape = shape, skew = skew)
}

#'@export
fiegarch_sim <- function(n = 750, nstart = 1000,
                         mu = 0, omega_sig = -9, phi = 0.6, psi = numeric(0),
                         kappa = -0.2, gamma = 0.3, d = 0.25, shape = 2,
                         df = 10,
                         skew = 1,
                         cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged")) {

  cond_dist <- match.arg(cond_dist)
  nn <- n + nstart          # pre-sample plus sample number of observations
  if ( cond_dist == "norm") {
    innov <- rnorm(nn + 1)    # we need innovations for our pre-sample and sample time points plus for one time point before the pre-sample
  } else if (cond_dist == "std") {
    innov <- rt(nn + 1, df = df) / sqrt(df / (df + 2))
  }


  # Obtain the probability density function for the innov. distribution
  # from the internal lookup table
  PDFunc <- lookup_table[["fun1"]][[cond_dist]]
  if (cond_dist %in% c("std", "sstd")) {
    shape <- df    # all pdfs are defined with "shape" argument (even if not used then within the function; for t-distr. types, we need to make shape our df)
  }

  # Define a function that follows the most basic definition of the expectation of
  # a continuous random variable (here: for the absolute innovations) but just for
  # the part behind the integration symbol
  Efun <- function(x, PDFunc, shape, skew) {
    abs(x) * PDFunc(x, shape = shape, skew = skew)
  }

  # Then we can integrate numerically over that function and thus
  # follow the definition of an expectation
  Econst <- stats::integrate(
    f = Efun, lower = -Inf, upper = Inf,
    subdivisions = 1000L, rel.tol = 1e-6, abs.tol = 1e-6,
    stop.on.error = FALSE, PDFunc = PDFunc, shape = shape,
    skew = skew
  )[[1]]

  # we need coefficients for the 1000 pre-pre-sample time points, for the pre-sample,
  # and for the sample with the exception of the very last sample time point
  m <- M <- nn + 1000 - 1
  if (m < 5000) {
    M <- 5000
  }
  coef_inf <- esemifar::farima_to_ar(ar = phi, ma = psi, d = d,
                                     max_i = M)[-1]

  C_omega <- 1 - sum(coef_inf)

  # problem: uncond. variance seems to be continuously (but slowly) decreasing over time

  sim <- fiegarch_sim_ar(
    innov = innov, omega_sig = omega_sig, coef_inf = head(coef_inf, m),
    kappa = kappa, gamma = gamma, n_out = n, mu = mu,
    Econst = Econst, C_omega = C_omega
  )

  out <- list(
    "rt" = ts(c(sim[["rt"]])),
    "sigt" = ts(c(sim[["sigt"]])),
    "etat" = ts(c(sim[["etat"]]))
  )

}
