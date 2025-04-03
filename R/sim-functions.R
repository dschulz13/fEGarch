cdf_fun <- function(x, dfun, ...) {
  vapply(
    x,
    FUN = function(x, dfun, ...) {
      integrate(f = dfun, lower = -Inf, upper = x,
                ..., subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6)[[1]]
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
#'@param P the number of Laplace distributions (minus 1) to derive the arithmetic
#'mean from as the basis for a (skewed) average Laplace (AL) distribution
#'distribution.
#'@param skew the skewness parameter in the skewed distributions.
#'
#'@details
#'Draw random samples from a normal distribution, \eqn{t}-distribution,
#'generalized error distribution, an average Laplace distribution,
#'or their skewed variants (all standardized to have mean zero and
#'variance one).
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
rald_s <- function(n, P = 8) {
  p <- runif(n)
  inv_cdf_fun(p = p, dfun = pdf_ald_v1, P = P)
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
#'@rdname sim_functions
rsald_s <- function(n, P = 8, skew = 1) {
  p <- runif(n)
  inv_cdf_fun(p = p, dfun = pdf_skew_sald_s, P = P, skew = skew)
}

#
# fiegarch_sim <- function(n = 750, nstart = 1000,
#                          mu = 0, omega_sig = -9, phi = 0.6, psi = numeric(0),
#                          kappa = -0.2, gamma = 0.3, d = 0.25, shape = 2,
#                          df = 10,
#                          skew = 1,
#                          cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged")) {
#
#   cond_dist <- match.arg(cond_dist)
#   nn <- n + nstart          # pre-sample plus sample number of observations
#   if ( cond_dist == "norm") {
#     innov <- rnorm(nn + 1)    # we need innovations for our pre-sample and sample time points plus for one time point before the pre-sample
#   } else if (cond_dist == "std") {
#     innov <- rt(nn + 1, df = df) / sqrt(df / (df + 2))
#   }
#
#
#   # Obtain the probability density function for the innov. distribution
#   # from the internal lookup table
#   PDFunc <- lookup_table[["fun1"]][[cond_dist]]
#   if (cond_dist %in% c("std", "sstd")) {
#     shape <- df    # all pdfs are defined with "shape" argument (even if not used then within the function; for t-distr. types, we need to make shape our df)
#   }
#
#   # Define a function that follows the most basic definition of the expectation of
#   # a continuous random variable (here: for the absolute innovations) but just for
#   # the part behind the integration symbol
#   Efun <- function(x, PDFunc, shape, skew) {
#     abs(x) * PDFunc(x, shape = shape, skew = skew)
#   }
#
#   # Then we can integrate numerically over that function and thus
#   # follow the definition of an expectation
#   Econst <- stats::integrate(
#     f = Efun, lower = -Inf, upper = Inf,
#     subdivisions = 1000L, rel.tol = 1e-6, abs.tol = 1e-6,
#     stop.on.error = FALSE, PDFunc = PDFunc, shape = shape,
#     skew = skew
#   )[[1]]
#
#   # we need coefficients for the 1000 pre-pre-sample time points, for the pre-sample,
#   # and for the sample with the exception of the very last sample time point
#   m <- M <- nn + 1000 - 1
#   if (m < 5000) {
#     M <- 5000
#   }
#   coef_inf <- esemifar::farima_to_ar(ar = phi, ma = psi, d = d,
#                                      max_i = M)[-1]
#
#   C_omega <- 1 - sum(coef_inf)
#
#   # problem: uncond. variance seems to be continuously (but slowly) decreasing over time
#
#   sim <- fiegarch_sim_ar(
#     innov = innov, omega_sig = omega_sig, coef_inf = head(coef_inf, m),
#     kappa = kappa, gamma = gamma, n_out = n, mu = mu,
#     Econst = Econst, C_omega = C_omega
#   )
#
#   out <- list(
#     "rt" = ts(c(sim[["rt"]])),
#     "sigt" = ts(c(sim[["sigt"]])),
#     "etat" = ts(c(sim[["etat"]]))
#   )
#
# }

#'Simulate From Models of the Broader EGARCH Family
#'
#'A streamlined simulation function to simulate from models
#'that are part of the broader EGARCH family specifiable
#'through \code{\link{fEGarch_spec}}.
#'
#'@param spec an object of class \code{"egarch_type_spec"} or
#'\code{"loggarch_type_spec"} returned by \code{\link{fEGarch_spec}}
#'or related wrapper functions; note that the model orders in the
#'conditional mean and in the
#'conditional variance as well as the long-memory settings are not
#'controlled through the element
#'\code{orders} and \code{long_memo} of \code{spec} but
#'solely through the parameter settings in \code{pars}.
#'@param pars a named list with the parameter specifications; the user
#'can provide a named list with only the settings they would like to adjust
#'relative to the default settings.
#'@param n the number of observations to return.
#'@param nstart the number of burn-in observations to simulate before
#'the final \code{n} values to keep; the first \code{nstart} values
#'are not returned; if a dual model, i.e. with model in the conditional
#'mean and in the conditional variance, is considered, two times \code{nstart}
#'is considered in the first simulation step in the conditional variance,
#'so that \code{n + nstart} values can be fed into the second simulation
#'step for the conditional mean.
#'@param trunc a truncation for the finite-order coefficient series
#'in long-memory models; can either be the character \code{"none"} for truncation
#'back to the very first observation at each time point, or to any positive integer
#'for setting the corresponding truncation length of the infinite-order representation
#'polynomial.
#'
#'@details
#'See the documentation on \code{\link{fEGarch_spec}} for information
#'on the models of the broader EGARCH family. This function provides
#'an easy way to simulate from these models.
#'
#'@return
#'A list with four elements is returned: \code{rt} are the simulated
#'observations, \code{etat} are the underlying innovations,
#'\code{sigt} are the correspondingly simulated conditional
#'standard deviations, and \code{cmeans} are the simulated
#'conditional means. These four elements are formatted as
#'\code{"ts"} class time series objects.
#'
#'@export
#'
#'@examples
#'spec <- megarch_spec()
#'sim <- fEGarch_sim(spec = spec)
#'mat <- do.call(cbind, sim)
#'plot(mat, main = "")
#'
fEGarch_sim <- function(spec = egarch_spec(),
                        pars = list(
                           mu = 0,
                           ar = numeric(0),
                           ma = numeric(0),
                           D = 0,
                           omega_sig = -9,
                           phi = 0.8,
                           psi = numeric(0),
                           kappa = -0.2,
                           gamma = 0.3,
                           d = 0,
                           df = 10,
                           shape = 2,
                           P = 3,
                           skew = 1
                         ), n = 1000, nstart = 5e3, trunc = "none") {

  nn <- n + nstart          # pre-sample plus sample number of observations
  N <- n

  trunc_vol <- trunc
  trunc_mean <- trunc

  # Parameter defaults (relevant if a truncated parameter list is provided)
  pars_default <- lookup_table$sim_pars_default
  names_pars <- names(pars)
  for (nam in names_pars) {
    pars_default[[nam]] <- pars[[nam]]
  }

  # Extract parameter settings
  mu <- pars_default$mu
  omega_sig <- pars_default$omega_sig
  phi <- pars_default$phi
  psi <- pars_default$psi
  kappa <- pars_default$kappa
  gamma <- pars_default$gamma
  d <- pars_default$d
  df <- pars_default$df
  shape <- pars_default$shape
  P <- pars_default$P
  skew <- pars_default$skew
  cond_dist <- spec@cond_dist
  lm <- d != 0

  # ARMA / FARIMA parameters
  ar <- pars_default$ar
  ma <- pars_default$ma
  D <- pars_default$D
  lm_arma <- D != 0

  # Simulate even more observations if ARMA / FARIMA required
  check <- length(ar) + length(ma) > 0 || lm_arma
  nn <- if (check) {
    N <- N + nstart
    nn + nstart
  } else {
    nn
  }

  if ((lm || lm_arma) && ((is.character(trunc) && trunc == "none") || (is.numeric(trunc) && trunc >= nn))) {
    trunc_vol <- nn - 1
    trunc_mean <- N - 1
  }

  if (cond_dist %in% c("std", "sstd")) {
    shape <- df    # all functions are defined with "shape" argument (even if not used then within the function; for t-distr. types, we need to make shape our df)
  } else if (cond_dist %in% c("ald", "sald")) {
    shape <- P
  }

  # Simulate innovations following selected distribution
  sim_fun <- lookup_table[["sim_funs"]][[cond_dist]]
  innov <- sim_fun(nn, shape = shape, skew = skew)

  PDF <- lookup_table$fun1[[cond_dist]]
  model_type <- c("logGARCH", "eGARCH")[[inherits(spec, "egarch_type_spec") + 1]]

  if (model_type == "eGARCH") {

  modulus <- spec@modulus
  powers <- spec@powers

  mod_a <- modulus[[1]]
  mod_m <- modulus[[2]]
  pow_a <- powers[[1]]
  pow_m <- powers[[2]]
  log_a <- (pow_a > 0) + 1
  log_m <- (pow_m > 0) + 1

  # Obtain the various functions to integrate over for the expectations
  # in the magnitude and asymmetry terms
  Func_a <- exp_funs[["asy"]][[mod_a + 1]][[log_a]]
  Func_m <- exp_funs[["mag"]][[mod_m + 1]][[log_m]]

  # Then we can integrate numerically over those function and thus
  # follow the definition of an expectation
  E_asy <- stats::integrate(
    f = Func_a, lower = -Inf, upper = Inf,
    subdivisions = 1000L, rel.tol = 1e-6, abs.tol = 1e-6,
    stop.on.error = FALSE, pow = pow_a, dfun = PDF, shape = shape,
    skew = skew
  )[[1]]

  E_mag <- stats::integrate(
    f = Func_m, lower = -Inf, upper = Inf,
    subdivisions = 1000L, rel.tol = 1e-6, abs.tol = 1e-6,
    stop.on.error = FALSE, pow = pow_m, dfun = PDF, shape = shape,
    skew = skew
  )[[1]]

  mode <- mode_mat[log_a, log_m]
  if (lm) {

  # we need coefficients for the 1 pre-pre-sample time points, for the pre-sample,
  # and for the sample with the exception of the very last sample time point

  coef_inf <- c(0, ma_infty(ar = phi, ma = psi, d = d,
                       max_i = trunc_vol - 1))

  np2 <- nextpow2(2 * nn - 1)
  # Compute the function g for the simulated innovations
  sim <- egarch_longmemo_type_sim(
    innov = innov, omega_sig = omega_sig, coef_inf = coef_inf,
    kappa = kappa, gamma = gamma, n_out = N, mu = mu,
    E_mag = E_mag, E_asy = E_asy, powers = powers,
    modulus = modulus, mode = mode, np2)

  } else {

  sim <- egarch_shortmemo_type_sim(
    innov = innov, omega_sig = omega_sig, ar_coef = phi,
    ma_coef = c(1, psi),
    kappa = kappa, gamma = gamma, n_out = N, mu = mu,
    E_mag = E_mag, E_asy = E_asy, powers = powers,
    modulus = modulus, mode = mode)

  }

  } else if (model_type == "logGARCH") {

    Func <- function(x, shape, skew, dfun) {
      log(x^2) * dfun(x, shape = shape, skew = skew)
    }

    E_const <- stats::integrate(
      f = Func, lower = -Inf, upper = Inf,
      subdivisions = 1000L, rel.tol = 1e-6, abs.tol = 1e-6,
      stop.on.error = FALSE, dfun = PDF, shape = shape,
      skew = skew
    )[[1]]

    if (lm) {

      np2 <- nextpow2(2 * nn - 1)

      coef_inf <- ma_infty(
        ar = phi, ma = psi, d = d, max_i = trunc_vol
      )
      coef_inf[[1]] <- 0

      sim <- loggarch_longmemo_type_sim(
        innov, omega_sig,
        coef_inf = coef_inf,
        n_out = N,
        mu = mu,
        Econst = E_const, np2 = np2)

    } else {

      n_ar <- length(phi)
      n_ma <- length(psi)
      if (n_ma > n_ar) {
        phi <- c(phi, rep(0, n_ma - n_ar))
      } else if (n_ar > n_ma) {
        psi <- c(psi, rep(0, n_ar - n_ma))
      }

      sim <- loggarch_shortmemo_type_sim(
        innov, omega_sig,
        ar_coef = phi, ma_coef = psi + phi,
        n_out = N,
        mu = mu,
        Econst = E_const)

    }

  }

  out <- list(
    "rt" = ts(c(sim[["rt"]])),
    "sigt" = ts(c(sim[["sigt"]])),
    "etat" = ts(c(sim[["etat"]])),
    "cmeans" = ts(rep(mu, n))
  )

  if (check) {
    rt <- out$rt - mu        # subtract mean again and then feed these as
                             # the errors into an ARMA / FARIMA model
    if (!lm_arma) {

      yt_data <- arma_sim_Cpp(innov = rt, ar = ar, ma = ma, nout = n, mu = mu)

    } else {

      coef_inf <- ar_infty(ar = ar, ma = ma, d = D, max_i = trunc_mean)[-1]
      yt_data <- farima_sim_Cpp(innov = rt, coef_inf = coef_inf, nout = n, mu = mu)

    }
    out[["rt"]] <- ts(c(yt_data[["rt"]]))
    out[["cmeans"]] <- ts(c(yt_data[["cmeans"]]))
    out[["etat"]] <- ts(tail(out[["etat"]], n))
    out[["sigt"]] <- ts(tail(out[["sigt"]], n))
  }

  out

}

#'Simulate From FIAPARCH Models
#'
#'A streamlined simulation function to simulate from
#'fractionally integrated asymmetric power autoregressive
#'conditional heteroskedasticity (FIAPARCH) models.
#'
#'@param pars a named list with the parameter specifications; the user
#'can provide a named list with only the settings they would like to adjust
#'relative to the default settings.
#'@param cond_dist a one-element character vector specifying
#'the conditional distribution to consider.
#'@param n the number of observations to return.
#'@param nstart the number of burn-in observations to simulate before
#'the final \code{n} values to keep; the first \code{nstart} values
#'are not returned; if a dual model, i.e. with model in the conditional
#'mean and in the conditional variance, is considered, two times \code{nstart}
#'is considered in the first simulation step in the conditional variance,
#'so that \code{n + nstart} values can be fed into the second simulation
#'step for the conditional mean.
#'
#'@details
#'See the documentation on \code{\link{fiaparch}} for information
#'on the FIAPARCH model. This function provides
#'an easy way to simulate from these models.
#'
#'@return
#'A list with four elements is returned: \code{rt} are the simulated
#'observations, \code{etat} are the underlying innovations,
#'\code{sigt} are the correspondingly simulated conditional
#'standard deviations, and \code{cmeans} are the simulated
#'conditional means. These four elements are formatted as
#'\code{"ts"} class time series objects.
#'
#'@export
#'
#'@examples
#'sim <- fiaparch_sim(n = 1000)
#'mat <- do.call(cbind, sim)
#'plot(mat, main = "")
#'

fiaparch_sim <- function(pars = list(
                           mu = 0,
                           ar = numeric(0),
                           ma = numeric(0),
                           D = 0,
                           omega = 0.0004,
                           phi = 0.2,
                           beta = 0.4,
                           gamma = 0.1,
                           delta = 2,
                           d = 0.25,
                           df = 10,
                           shape = 2,
                           P = 3,
                           skew = 1
                         ), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                         n = 1000, nstart = 5e3, trunc = "none") {

  cond_dist <- match.arg(cond_dist)

  nn <- n + nstart          # pre-sample plus sample number of observations
  N <- n

  trunc_vol <- trunc
  trunc_mean <- trunc

  # Parameter defaults (relevant if a truncated parameter list is provided)
  pars_default <- lookup_table$sim_pars_default_fiaparch
  names_pars <- names(pars)
  for (nam in names_pars) {
    pars_default[[nam]] <- pars[[nam]]
  }

  # Extract parameter settings
  mu <- pars_default$mu
  omega <- pars_default$omega
  phi <- pars_default$phi
  beta <- -pars_default$beta      # "-" required to create a polynomial with positive signs
  gamma <- pars_default$gamma
  delta <- pars_default$delta
  d <- pars_default$d
  df <- pars_default$df
  P <- pars_default$P
  shape <- pars_default$shape
  skew <- pars_default$skew

  # ARMA / FARIMA parameters
  ar <- pars_default$ar
  ma <- pars_default$ma
  D <- pars_default$D
  lm_arma <- D != 0

  # Simulate even more observations if ARMA / FARIMA required
  check <- length(ar) + length(ma) > 0 || lm_arma
  nn <- if (check) {
    N <- N + nstart
    nn + nstart
  } else {
    nn
  }

  lm <- TRUE

  if ((lm || lm_arma) && ((is.character(trunc) && trunc == "none") || (is.numeric(trunc) && trunc >= nn))) {
    trunc_vol <- nn - 1
    trunc_mean <- N - 1
  }

  if (cond_dist %in% c("std", "sstd")) {
    shape <- df    # all functions are defined with "shape" argument (even if not used then within the function; for t-distr. types, we need to make shape our df)
  } else if (cond_dist %in% c("ald", "sald")) {
    shape <- P
  }

  # Simulate innovations following selected distribution
  sim_fun <- lookup_table[["sim_funs"]][[cond_dist]]
  innov <- sim_fun(nn, shape = shape, skew = skew)

  # we need coefficients for the 1 pre-pre-sample time points, for the pre-sample,
  # and for the sample with the exception of the very last sample time point

  coef_inf <- ar_infty(ar = phi, ma = beta, d = d,
                                     max_i = trunc_vol)
  coef_inf[[1]] <- 0

  sim <- fiaparch_sim_Cpp(
    innov = innov, omega = omega, coef_inf = coef_inf,
    gamma = gamma, delta = delta, n_out = N, mu = mu)

  out <- list(
    "rt" = ts(c(sim[["rt"]])),
    "sigt" = ts(c(sim[["sigt"]])),
    "etat" = ts(c(sim[["etat"]])),
    "cmeans" = ts(rep(mu, n))
  )

  if (check) {
    rt <- out$rt - mu        # subtract mean again and then feed these as
                             # the errors into an ARMA / FARIMA model
    if (!lm_arma) {

      yt_data <- arma_sim_Cpp(innov = rt, ar = ar, ma = ma, nout = n, mu = mu)

    } else {

      coef_inf <- ar_infty(ar = ar, ma = ma, d = D, max_i = trunc_mean)[-1]
      yt_data <- farima_sim_Cpp(innov = rt, coef_inf = coef_inf, nout = n, mu = mu)

    }
    out[["rt"]] <- ts(c(yt_data[["rt"]]))
    out[["cmeans"]] <- ts(c(yt_data[["cmeans"]]))
    out[["etat"]] <- ts(tail(out[["etat"]], n))
    out[["sigt"]] <- ts(tail(out[["sigt"]], n))
  }

  out

}

#'Simulate From FIGJRGARCH Models
#'
#'A streamlined simulation function to simulate from
#'FIGJR-GARCH models.
#'
#'@param pars a named list with the parameter specifications; the user
#'can provide a named list with only the settings they would like to adjust
#'relative to the default settings.
#'@param cond_dist a one-element character vector specifying
#'the conditional distribution to consider.
#'@param n the number of observations to return.
#'@param nstart the number of burn-in observations to simulate before
#'the final \code{n} values to keep; the first \code{nstart} values
#'are not returned; if a dual model, i.e. with model in the conditional
#'mean and in the conditional variance, is considered, two times \code{nstart}
#'is considered in the first simulation step in the conditional variance,
#'so that \code{n + nstart} values can be fed into the second simulation
#'step for the conditional mean.
#'
#'@details
#'See the documentation on \code{\link{figjrgarch}} for information
#'on the FIGJRGARCH model. This function provides
#'an easy way to simulate from these models.
#'
#'@return
#'A list with four elements is returned: \code{rt} are the simulated
#'observations, \code{etat} are the underlying innovations,
#'\code{sigt} are the correspondingly simulated conditional
#'standard deviations, and \code{cmeans} are the simulated
#'conditional means. These four elements are formatted as
#'\code{"ts"} class time series objects.
#'
#'@export
#'
#'@examples
#'sim <- figjrgarch_sim(n = 1000)
#'mat <- do.call(cbind, sim)
#'plot(mat, main = "")
#'

figjrgarch_sim <- function(pars = list(
                           mu = 0,
                           ar = numeric(0),
                           ma = numeric(0),
                           D = 0,
                           omega = 0.0004,
                           phi = 0.2,
                           beta = 0.4,
                           gamma = 0.1,
                           d = 0.25,
                           df = 10,
                           shape = 2,
                           P = 3,
                           skew = 1
                         ), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                         n = 1000, nstart = 5e3, trunc = "none") {

  pars[["delta"]] <- 2
  fiaparch_sim(pars = pars, cond_dist = cond_dist, n = n, nstart = nstart, trunc = trunc)

}

#'Simulate From FIGARCH Models
#'
#'A streamlined simulation function to simulate from
#'fractionally integrated generalized autoregressive conditional
#'heteroskedasticity (FIGARCH) models.
#'
#'@param pars a named list with the parameter specifications; the user
#'can provide a named list with only the settings they would like to adjust
#'relative to the default settings.
#'@param cond_dist a one-element character vector specifying
#'the conditional distribution to consider.
#'@param n the number of observations to return.
#'@param nstart the number of burn-in observations to simulate before
#'the final \code{n} values to keep; the first \code{nstart} values
#'are not returned; if a dual model, i.e. with model in the conditional
#'mean and in the conditional variance, is considered, two times \code{nstart}
#'is considered in the first simulation step in the conditional variance,
#'so that \code{n + nstart} values can be fed into the second simulation
#'step for the conditional mean.
#'@param trunc a truncation for the finite-order coefficient series
#'in long-memory models; can either be the character \code{"none"} for truncation
#'back to the very first observation at each time point, or to any positive integer
#'for setting the corresponding truncation length of the infinite-order representation
#'polynomial.
#'
#'@details
#'See the documentation on \code{\link{figarch}} for information
#'on the FIGARCH model. This function provides
#'an easy way to simulate from these models.
#'
#'@return
#'A list with four elements is returned: \code{rt} are the simulated
#'observations, \code{etat} are the underlying innovations,
#'\code{sigt} are the correspondingly simulated conditional
#'standard deviations, and \code{cmeans} are the simulated
#'conditional means. These four elements are formatted as
#'\code{"ts"} class time series objects.
#'
#'@export
#'
#'@examples
#'sim <- figarch_sim(n = 1000)
#'mat <- do.call(cbind, sim)
#'plot(mat, main = "")
#'

figarch_sim <- function(pars = list(
                           mu = 0,
                           ar = numeric(0),
                           ma = numeric(0),
                           D = 0,
                           omega = 0.0004,
                           phi = 0.2,
                           beta = 0.4,
                           d = 0.25,
                           df = 10,
                           shape = 2,
                           P = 3,
                           skew = 1
                         ), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                         n = 1000, nstart = 5e3, trunc = "none") {

  cond_dist <- match.arg(cond_dist)

  nn <- n + nstart          # pre-sample plus sample number of observations
  N <- n

  trunc_vol <- trunc
  trunc_mean <- trunc

  # Parameter defaults (relevant if a truncated parameter list is provided)
  pars_default <- lookup_table$sim_pars_default_figarch
  names_pars <- names(pars)
  for (nam in names_pars) {
    pars_default[[nam]] <- pars[[nam]]
  }

  # Extract parameter settings
  mu <- pars_default$mu
  omega <- pars_default$omega
  phi <- pars_default$phi
  beta <- -pars_default$beta      # "-" required to create a polynomial with positive signs
  d <- pars_default$d
  df <- pars_default$df
  P <- pars_default$P
  shape <- pars_default$shape
  skew <- pars_default$skew

  # ARMA / FARIMA parameters
  ar <- pars_default$ar
  ma <- pars_default$ma
  D <- pars_default$D
  lm_arma <- D != 0

  # Simulate even more observations if ARMA / FARIMA required
  check <- length(ar) + length(ma) > 0 || lm_arma
  nn <- if (check) {
    N <- N + nstart
    nn + nstart
  } else {
    nn
  }

  lm <- TRUE

  if ((lm || lm_arma) && ((is.character(trunc) && trunc == "none") || (is.numeric(trunc) && trunc >= nn))) {
    trunc_vol <- nn - 1
    trunc_mean <- N - 1
  }

  if (cond_dist %in% c("std", "sstd")) {
    shape <- df    # all functions are defined with "shape" argument (even if not used then within the function; for t-distr. types, we need to make shape our df)
  } else if (cond_dist %in% c("ald", "sald")) {
    shape <- P
  }

  # Simulate innovations following selected distribution
  sim_fun <- lookup_table[["sim_funs"]][[cond_dist]]
  innov <- sim_fun(nn, shape = shape, skew = skew)

  # we need coefficients for the 1 pre-pre-sample time points, for the pre-sample,
  # and for the sample with the exception of the very last sample time point

  coef_inf <- ar_infty(ar = phi, ma = beta, d = d,
                                     max_i = trunc_vol)
  coef_inf[[1]] <- 0

  sim <- figarch_sim_Cpp(
    innov = innov, omega = omega, coef_inf = coef_inf,
    n_out = N, mu = mu)

  out <- list(
    "rt" = ts(c(sim[["rt"]])),
    "sigt" = ts(c(sim[["sigt"]])),
    "etat" = ts(c(sim[["etat"]])),
    "cmeans" = ts(rep(mu, n))
  )

  if (check) {
    rt <- out$rt - mu        # subtract mean again and then feed these as
                             # the errors into an ARMA / FARIMA model
    if (!lm_arma) {

      yt_data <- arma_sim_Cpp(innov = rt, ar = ar, ma = ma, nout = n, mu = mu)

    } else {

      coef_inf <- ar_infty(ar = ar, ma = ma, d = D, max_i = trunc_mean)[-1]
      yt_data <- farima_sim_Cpp(innov = rt, coef_inf = coef_inf, nout = n, mu = mu)

    }
    out[["rt"]] <- ts(c(yt_data[["rt"]]))
    out[["cmeans"]] <- ts(c(yt_data[["cmeans"]]))
    out[["etat"]] <- ts(tail(out[["etat"]], n))
    out[["sigt"]] <- ts(tail(out[["sigt"]], n))
  }

  out

}


#'Simulate From APARCH Models
#'
#'A streamlined simulation function to simulate from
#'asymmetric power autoregressive
#'conditional heteroskedasticity (APARCH) models.
#'
#'@param pars a named list with the parameter specifications; the user
#'can provide a named list with only the settings they would like to adjust
#'relative to the default settings.
#'@param cond_dist a one-element character vector specifying
#'the conditional distribution to consider.
#'@param n the number of observations to return.
#'@param nstart the number of burn-in observations to simulate before
#'the final \code{n} values to keep; the first \code{nstart} values
#'are not returned; if a dual model, i.e. with model in the conditional
#'mean and in the conditional variance, is considered, two times \code{nstart}
#'is considered in the first simulation step in the conditional variance,
#'so that \code{n + nstart} values can be fed into the second simulation
#'step for the conditional mean.
#'@param trunc a truncation for the finite-order coefficient series
#'in long-memory models; can either be the character \code{"none"} for truncation
#'back to the very first observation at each time point, or to any positive integer
#'for setting the corresponding truncation length of the infinite-order representation
#'polynomial.
#'
#'@details
#'See the documentation on \code{\link{aparch}} for information
#'on the APARCH model. This function provides
#'an easy way to simulate from these models.
#'
#'@return
#'A list with four elements is returned: \code{rt} are the simulated
#'observations, \code{etat} are the underlying innovations,
#'\code{sigt} are the correspondingly simulated conditional
#'standard deviations, and \code{cmeans} are the simulated
#'conditional means. These four elements are formatted as
#'\code{"ts"} class time series objects.
#'
#'@export
#'
#'@examples
#'sim <- aparch_sim(n = 1000)
#'mat <- do.call(cbind, sim)
#'plot(mat, main = "")
#'

aparch_sim <- function(pars = list(
                           mu = 0,
                           ar = numeric(0),
                           ma = numeric(0),
                           D = 0,
                           omega = 0.0004,
                           phi = 0.05,
                           beta = 0.8,
                           gamma = 0.1,
                           delta = 2,
                           df = 10,
                           shape = 2,
                           P = 3,
                           skew = 1
                         ), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                         n = 1000, nstart = 5e3, trunc = "none") {

  cond_dist <- match.arg(cond_dist)

  nn <- n + nstart          # pre-sample plus sample number of observations
  N <- n

  trunc_mean <- trunc

  # Parameter defaults (relevant if a truncated parameter list is provided)
  pars_default <- lookup_table$sim_pars_default_aparch
  names_pars <- names(pars)
  for (nam in names_pars) {
    pars_default[[nam]] <- pars[[nam]]
  }

  # Extract parameter settings
  mu <- pars_default$mu
  omega <- pars_default$omega
  phi <- pars_default$phi
  beta <- pars_default$beta
  gamma <- pars_default$gamma
  delta <- pars_default$delta
  df <- pars_default$df
  P <- pars_default$P
  shape <- pars_default$shape
  skew <- pars_default$skew

  l_phi <- length(phi)
  l_gamma <- length(gamma)

  stopifnot("Input parameter vectors phi and gamma must be of the same length." = l_phi == l_gamma)

  # ARMA / FARIMA parameters
  ar <- pars_default$ar
  ma <- pars_default$ma
  D <- pars_default$D
  lm_arma <- D != 0

  # Simulate even more observations if ARMA / FARIMA required
  check <- length(ar) + length(ma) > 0 || lm_arma
  nn <- if (check) {
    N <- N + nstart
    nn + nstart
  } else {
    nn
  }

  if (lm_arma && ((is.character(trunc) && trunc == "none") || (is.numeric(trunc) && trunc >= nn))) {
    trunc_mean <- N - 1
  }

  if (cond_dist %in% c("std", "sstd")) {
    shape <- df    # all functions are defined with "shape" argument (even if not used then within the function; for t-distr. types, we need to make shape our df)
  } else if (cond_dist %in% c("ald", "sald")) {
    shape <- P
  }

  # Simulate innovations following selected distribution
  sim_fun <- lookup_table[["sim_funs"]][[cond_dist]]
  innov <- sim_fun(nn, shape = shape, skew = skew)

  dfun <- lookup_table[["fun1"]][[cond_dist]]
  int_fun <- function(x, shape, skew, gamma_i, delta, dfun) {
    (abs(x) - gamma_i * x)^delta * dfun(x, shape = shape, skew = skew)
  }
  E_g <- rep(NA, length(gamma))
  for (i in seq_along(gamma)) {
    E_g <- stats::integrate(
      f = int_fun, lower = -Inf, upper = Inf, subdivisions = 1000L,
      rel.tol = 1e-6, abs.tol = 1e-6,
      shape = shape, skew = skew, gamma_i = gamma[[i]], delta = delta,
      dfun = dfun
    )[[1]]
  }
  E_sigd <- omega / (1 - sum(phi * E_g) - sum(beta))

  sim <- aparch_sim_Cpp(
    innov = innov, omega = omega, phi = phi, beta = beta,
    gamma = gamma, delta = delta, n_out = N, mu = mu,
    E_sigd = E_sigd)

  out <- list(
    "rt" = ts(c(sim[["rt"]])),
    "sigt" = ts(c(sim[["sigt"]])),
    "etat" = ts(c(sim[["etat"]])),
    "cmeans" = ts(rep(mu, n))
  )

  if (check) {
    rt <- out$rt - mu        # subtract mean again and then feed these as
                             # the errors into an ARMA / FARIMA model
    if (!lm_arma) {

      yt_data <- arma_sim_Cpp(innov = rt, ar = ar, ma = ma, nout = n, mu = mu)

    } else {

      coef_inf <- ar_infty(ar = ar, ma = ma, d = D, max_i = N)[-1]
      yt_data <- farima_sim_Cpp(innov = rt, coef_inf = coef_inf, nout = n, mu = mu)

    }
    out[["rt"]] <- ts(c(yt_data[["rt"]]))
    out[["cmeans"]] <- ts(c(yt_data[["cmeans"]]))
    out[["etat"]] <- ts(tail(out[["etat"]], n))
    out[["sigt"]] <- ts(tail(out[["sigt"]], n))
  }

  out

}

#'Simulate From GJR-GARCH Models
#'
#'A streamlined simulation function to simulate from
#'GJR-GARCH models.
#'
#'@param pars a named list with the parameter specifications; the user
#'can provide a named list with only the settings they would like to adjust
#'relative to the default settings.
#'@param cond_dist a one-element character vector specifying
#'the conditional distribution to consider.
#'@param n the number of observations to return.
#'@param nstart the number of burn-in observations to simulate before
#'the final \code{n} values to keep; the first \code{nstart} values
#'are not returned; if a dual model, i.e. with model in the conditional
#'mean and in the conditional variance, is considered, two times \code{nstart}
#'is considered in the first simulation step in the conditional variance,
#'so that \code{n + nstart} values can be fed into the second simulation
#'step for the conditional mean.
#'@param trunc a truncation for the finite-order coefficient series
#'in long-memory models; can either be the character \code{"none"} for truncation
#'back to the very first observation at each time point, or to any positive integer
#'for setting the corresponding truncation length of the infinite-order representation
#'polynomial.
#'
#'@details
#'See the documentation on \code{\link{gjrgarch}} for information
#'on the GJR-GARCH model. This function provides
#'an easy way to simulate from these models.
#'
#'@return
#'A list with four elements is returned: \code{rt} are the simulated
#'observations, \code{etat} are the underlying innovations,
#'\code{sigt} are the correspondingly simulated conditional
#'standard deviations, and \code{cmeans} are the simulated
#'conditional means. These four elements are formatted as
#'\code{"ts"} class time series objects.
#'
#'@export
#'
#'@examples
#'sim <- gjrgarch_sim(n = 1000)
#'mat <- do.call(cbind, sim)
#'plot(mat, main = "")
#'

gjrgarch_sim <- function(pars = list(
                           mu = 0,
                           ar = numeric(0),
                           ma = numeric(0),
                           D = 0,
                           omega = 0.0004,
                           phi = 0.05,
                           beta = 0.8,
                           gamma = 0.1,
                           df = 10,
                           shape = 2,
                           P = 3,
                           skew = 1
                         ), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                         n = 1000, nstart = 5e3, trunc = "none") {


  pars[["delta"]] <- 2   # Fix delta at 2 and use APARCH functions
  aparch_sim(
    pars = pars,
    cond_dist = cond_dist,
    n = n, nstart = nstart, trunc = trunc
  )


}


#'Simulate From GARCH Models
#'
#'A streamlined simulation function to simulate from
#'generalized autoregressive conditional heteroskedasticity
#'(GARCH) models.
#'
#'@param pars a named list with the parameter specifications; the user
#'can provide a named list with only the settings they would like to adjust
#'relative to the default settings.
#'@param cond_dist a one-element character vector specifying
#'the conditional distribution to consider.
#'@param n the number of observations to return.
#'@param nstart the number of burn-in observations to simulate before
#'the final \code{n} values to keep; the first \code{nstart} values
#'are not returned; if a dual model, i.e. with model in the conditional
#'mean and in the conditional variance, is considered, two times \code{nstart}
#'is considered in the first simulation step in the conditional variance,
#'so that \code{n + nstart} values can be fed into the second simulation
#'step for the conditional mean.
#'@param trunc a truncation for the finite-order coefficient series
#'in long-memory models; can either be the character \code{"none"} for truncation
#'back to the very first observation at each time point, or to any positive integer
#'for setting the corresponding truncation length of the infinite-order representation
#'polynomial.
#'
#'@details
#'See the documentation on \code{\link{garch}} for information
#'on the GARCH model. This function provides
#'an easy way to simulate from these models.
#'
#'@return
#'A list with four elements is returned: \code{rt} are the simulated
#'observations, \code{etat} are the underlying innovations,
#'\code{sigt} are the correspondingly simulated conditional
#'standard deviations, and \code{cmeans} are the simulated
#'conditional means. These four elements are formatted as
#'\code{"ts"} class time series objects.
#'
#'@export
#'
#'@examples
#'sim <- garch_sim(n = 1000)
#'mat <- do.call(cbind, sim)
#'plot(mat, main = "")
#'

garch_sim <- function(pars = list(
                           mu = 0,
                           ar = numeric(0),
                           ma = numeric(0),
                           D = 0,
                           omega = 0.0004,
                           phi = 0.05,
                           beta = 0.8,
                           df = 10,
                           shape = 2,
                           P = 3,
                           skew = 1
                         ), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                         n = 1000, nstart = 5e3, trunc = "none") {

  cond_dist <- match.arg(cond_dist)

  nn <- n + nstart          # pre-sample plus sample number of observations
  N <- n

  trunc_mean <- trunc

  # Parameter defaults (relevant if a truncated parameter list is provided)
  pars_default <- lookup_table$sim_pars_default_aparch
  names_pars <- names(pars)
  for (nam in names_pars) {
    pars_default[[nam]] <- pars[[nam]]
  }

  # Extract parameter settings
  mu <- pars_default$mu
  omega <- pars_default$omega
  phi <- pars_default$phi
  beta <- pars_default$beta
  df <- pars_default$df
  P <- pars_default$P
  shape <- pars_default$shape
  skew <- pars_default$skew

  # ARMA / FARIMA parameters
  ar <- pars_default$ar
  ma <- pars_default$ma
  D <- pars_default$D
  lm_arma <- D != 0

  # Simulate even more observations if ARMA / FARIMA required
  check <- length(ar) + length(ma) > 0 || lm_arma
  nn <- if (check) {
    N <- N + nstart
    nn + nstart
  } else {
    nn
  }

  if (lm_arma && ((is.character(trunc) && trunc == "none") || (is.numeric(trunc) && trunc >= nn))) {
    trunc_mean <- N - 1
  }

  if (cond_dist %in% c("std", "sstd")) {
    shape <- df    # all functions are defined with "shape" argument (even if not used then within the function; for t-distr. types, we need to make shape our df)
  } else if (cond_dist %in% c("ald", "sald")) {
    shape <- P
  }

  # Simulate innovations following selected distribution
  sim_fun <- lookup_table[["sim_funs"]][[cond_dist]]
  innov <- sim_fun(nn, shape = shape, skew = skew)

  EY2 <- omega / (1 - sum(phi) - sum(beta))

  sim <- garch_sim_Cpp(
    innov = innov, omega = omega, phi = phi, beta = beta,
    n_out = N, mu = mu, EY2 = EY2)

  out <- list(
    "rt" = ts(c(sim[["rt"]])),
    "sigt" = ts(c(sim[["sigt"]])),
    "etat" = ts(c(sim[["etat"]])),
    "cmeans" = ts(rep(mu, n))
  )

  if (check) {
    rt <- out$rt - mu        # subtract mean again and then feed these as
                             # the errors into an ARMA / FARIMA model
    if (!lm_arma) {

      yt_data <- arma_sim_Cpp(innov = rt, ar = ar, ma = ma, nout = n, mu = mu)

    } else {

      coef_inf <- ar_infty(ar = ar, ma = ma, d = D, max_i = trunc_mean)[-1]
      yt_data <- farima_sim_Cpp(innov = rt, coef_inf = coef_inf, nout = n, mu = mu)

    }
    out[["rt"]] <- ts(c(yt_data[["rt"]]))
    out[["cmeans"]] <- ts(c(yt_data[["cmeans"]]))
    out[["etat"]] <- ts(tail(out[["etat"]], n))
    out[["sigt"]] <- ts(tail(out[["sigt"]], n))
  }

  out

}
