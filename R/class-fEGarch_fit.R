
setClass("fEGarch_fit",
  slots = c(
    pars = "numeric",
    se = "numeric",
    vcov_mat = "matrix",
    rt = "ANY",
    cmeans = "ANY",
    sigt = "ANY",
    etat = "ANY",
    llhood = "numeric",
    inf_criteria = "numeric",
    cond_dist = "character",
    orders = "numeric",
    long_memo = "logical",
    armaspec = "arma_spec"
  ),
  prototype = list(
    pars = numeric(1),
    se = numeric(1),
    vcov_mat = matrix(numeric(1)),
    rt = numeric(1),
    cmeans = numeric(1),
    sigt = numeric(1),
    etat = numeric(1),
    llhood = numeric(1),
    inf_criteria = numeric(1),
    cond_dist = character(1),
    orders = numeric(1),
    long_memo = logical(1),
    armaspec = arma_spec()
  )
)

setClass("fEGarch_fit_egarch",
  contains = "fEGarch_fit",
  slots = c(
    powers = "numeric",
    modulus = "logical"
  ),
  prototype = list(
    orders = numeric(1),
    modulus = logical(1)
  )
)

setClass("fEGarch_fit_loggarch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_fiaparch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_lmgarch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_rugarch_wrapper",
  contains = "fEGarch_fit",
  slots = c(
    rugarch_model = "ANY",
    model = "character"
  ),
  prototype = list(
    rugarch_model = numeric(1),
    model = character(1)
  )
)

fEGarch_fit_egarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    cmeans,
    sigt,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    powers,
    modulus,
    long_memo,
    armaspec
) {

  new("fEGarch_fit_egarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    cmeans = cmeans,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    powers = powers,
    modulus = modulus,
    long_memo = long_memo,
    armaspec = armaspec
  )

}

fEGarch_fit_loggarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    armaspec
) {

  new("fEGarch_fit_loggarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    armaspec = armaspec
  )

}


fEGarch_fit_fiaparch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    armaspec
) {

  new("fEGarch_fit_fiaparch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    armaspec = armaspec
  )

}

fEGarch_fit_lmgarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    armaspec
) {

  new("fEGarch_fit_lmgarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    armaspec = armaspec
  )

}

fEGarch_fit_rugarch_wrapper <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    rugarch_model,
    model,
    armaspec
) {

  new("fEGarch_fit_rugarch_wrapper",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    rugarch_model = rugarch_model,
    model = model,
    armaspec = armaspec
  )

}

#================================================#

fetch_arma_pars_fit <- function(object) {

  lm_arma <-  object@armaspec@long_memo
  arma_orders <- object@armaspec@orders

  p_ar <- arma_orders[[1]]
  q_ma <- arma_orders[[2]]

  pg0 <- p_ar > 0
  qg0 <- q_ma > 0
  pqg0 <- p_ar + q_ma > 0

  pars <- object@pars

  # Get all estimated parameters from output

  ar <- if (pg0) {
    names_ar <- paste0("ar", 1:p_ar)
    pars[names_ar]
  } else {
    numeric(0)
  }

  ma <- if (qg0) {
    names_ma <- paste0("ma", 1:q_ma)
    pars[names_ma]
  } else {
    numeric(0)
  }

  D <- if (lm_arma) pars[["D"]] else numeric(0)

  list(
   ar = ar,
   ma = ma,
   D = D,
   pqg0 = pqg0,
   lm_arma = lm_arma
  )

}

arma_farima_forecasts <- function(ar, ma, D, pqg0, lm_arma, mu, x, innov, n.ahead) {

  mean_fc <- if (pqg0 || lm_arma) {

    if (!lm_arma) {
      forecast_arma_Cpp(
        x = x,
        et = innov,
        mu = mu,
        ma_e = ma,
        ar_e = ar,
        horizon = n.ahead
      )
    } else {
      coef_inf_arma <- esemifar::farima_to_ar(ar = ar, ma = ma, d = D, max_i = length(x) + n.ahead - 1)[-1]
      forecast_farima_Cpp(
        x = x,
        mu = mu,
        horizon = n.ahead,
        coef_inf = coef_inf_arma
      )
    }

  } else {

    rep(mu, n.ahead)

  }

  mean_fc
}

# Forecasting methods

# For EGARCH-type
methods::setMethod("fEGarch_predict", "fEGarch_fit_egarch",
  function(object, n.ahead = 1, ...) {

  lm_garch <- object@long_memo
  garch_orders <- object@orders
  p <- garch_orders[[1]]
  q <- garch_orders[[2]]

  names_phi <- paste0("phi", 1:p)
  phi <- pars[names_phi]

  psi <- if (q > 1) {
    names_psi <- paste0("psi", 1:(q - 1))
    pars[names_psi]
  } else {
    numeric(0)
  }
  psi_2 <- psi
  psi <- c(1, psi)
  kappa <- pars[["kappa"]]
  gamma <- pars[["gamma"]]

  d <- if (lm_garch) pars[["d"]] else numeric(0)

  cond_dist <- object@cond_dist

  skew <- if (cond_dist %in% c("snorm", "sstd", "sged", "sald")) {
    pars[["skew"]]
  } else {
    numeric(0)
  }

  shape <- if (cond_dist %in% c("std", "sstd", "ged", "sged", "ald", "sald")) {
    lp <- length(pars)
    pars[[lp - 1]]
  } else  {
    numeric(0)
  }

  mu <- if (exists("mu", where = as.list(pars))) {
    pars[["mu"]]
  } else {
    0
  }

  # ARMA / FARIMA specifications
  arma_specs <- fetch_arma_pars_fit(object)
  ar <- arma_specs$ar
  ma <- arma_specs$ma
  D <- arma_specs$D
  pqg0 <- arma_specs$pqg0
  lm_arma <- arma_specs$lm_arma

  # EGARCH settings
  powers <- object@powers
  modulus <- object@modulus

  dfun <- lookup_table[["fun1"]][[cond_dist]]

  pow1 <- powers[[1]]
  pow1g0 <- pow1 > 0
  mod1 <- modulus[[1]]
  fun_asy <- exp_funs[["asy"]][[mod1 + 1]][[pow1g0 + 1]]
  E_asy <- stats::integrate(fun_asy, lower = -Inf, upper = Inf,
                            subdivisions = 1000L, rel.tol = 1e-6,
                            abs.tol = 1e-6,
                            pow = pow1, shape = shape, skew = skew, dfun = dfun)[[1]]

  pow2 <- powers[[2]]
  pow2g0 <- pow2 > 0
  mod2 <- modulus[[2]]
  fun_mag <- exp_funs[["mag"]][[mod2 + 1]][[pow2g0 + 1]]
  E_mag <- stats::integrate(fun_mag, lower = -Inf, upper = Inf,
                            subdivisions = 1000L, rel.tol = 1e-6,
                            abs.tol = 1e-6,
                            pow = pow1, shape = shape, skew = skew, dfun = dfun)[[1]]

  Elnsig2 <- pars[["omega_sig"]]

  # ARMA / FARIMA point forecasts of cond. mean
  mean_fc <- arma_farima_forecasts(ar = ar, ma = ma, D = D, pqg0 = pqg0,
                                   lm_arma = lm_arma, mu = mu,
                                   x = zoo::coredata(object@rt),
                                   innov = zoo::coredata(object@rt - object@cmeans),
                                   n.ahead = n.ahead)

  # EGARCH point forecasts of cond. volatility

  mode <- mode_mat[pow1g0 + 1, pow2g0 + 1]

  sigt_fc <- if (!lm_garch) {
    alpha <- psi * kappa
    beta <- psi * gamma
    omega <- Elnsig2 * (1 - sum(phi))
    sigt_egarch_forecast_shortCpp(
      et = zoo::coredata(object@etat),
      sigt = zoo::coredata(object@sigt),
      omega = omega,
      phi = phi,
      alpha = alpha,
      beta = beta,
      E_asy = E_asy,
      E_mag = E_mag,
      powers = powers,
      modulus = modulus,
      mode = mode,
      horizon = n.ahead
    )
  } else {
    coef_inf_garch <- esemifar::farima_to_ma(
      ar = phi, ma = psi_2, d = d, max_i = length(object@rt) + horizon - 2
    )
    sigt_egarch_forecast_longCpp(
      et = zoo::coredata(object@etat),
      coef_inf = coef_inf_garch,
      kappa = kappa,
      gamma = gamma,
      E_asy = E_asy,
      E_mag = E_mag,
      Elnsig2 = Elnsig2,
      powers = powers,
      modulus = modulus,
      mode = mode,
      horizon = horizon
    )
  }

  list(
    cmeans = mean_fc,
    sigt = sigt_fc
  )

  }
)

# For Log-GARCH-type
methods::setMethod("fEGarch_predict", "fEGarch_fit_loggarch",
  function(object, n.ahead = 1, ...) {

  lm_garch <- object@long_memo
  garch_orders <- object@orders
  p <- garch_orders[[1]]
  q <- garch_orders[[2]]

  names_phi <- paste0("phi", 1:p)
  phi <- pars[names_phi]

  names_psi <- paste0("psi", 1:q)
  psi <- pars[names_psi]

  d <- if (lm_garch) pars[["d"]] else numeric(0)

  cond_dist <- object@cond_dist

  skew <- if (cond_dist %in% c("snorm", "sstd", "sged", "sald")) {
    pars[["skew"]]
  } else {
    numeric(0)
  }

  shape <- if (cond_dist %in% c("std", "sstd", "ged", "sged", "ald", "sald")) {
    lp <- length(pars)
    pars[[lp - 1]]
  } else  {
    numeric(0)
  }

  mu <- if (exists("mu", where = as.list(pars))) {
    pars[["mu"]]
  } else {
    0
  }

  # ARMA / FARIMA specifications
  arma_specs <- fetch_arma_pars_fit(object)
  ar <- arma_specs$ar
  ma <- arma_specs$ma
  D <- arma_specs$D
  pqg0 <- arma_specs$pqg0
  lm_arma <- arma_specs$lm_arma

  dfun <- lookup_table[["fun1"]][[cond_dist]]

  int_fun <- function(x, shape, skew, dfun) {
    log(x^2) * dfun(x, shape = shape, skew = skew)
  }

  Elneta2 <- stats::integrate(
    int_fun, lower = -Inf, upper = Inf,
    subdivisions = 1000L, rel.tol = 1e-6, abs.tol = 1e-6,
    shape = shape, skew = skew, dfun = dfun
  )[[1]]

  Elnsig2 <- pars[["omega_sig"]]

  # ARMA / FARIMA point forecasts of cond. mean
  mean_fc <- arma_farima_forecasts(ar = ar, ma = ma, D = D, pqg0 = pqg0,
                                   lm_arma = lm_arma, mu = mu,
                                   x = zoo::coredata(object@rt),
                                   innov = zoo::coredata(object@rt - object@cmeans),
                                   n.ahead = n.ahead)

  # EGARCH point forecasts of cond. volatility

  mode <- mode_mat[pow1g0 + 1, pow2g0 + 1]

  sigt_fc <- if (!lm_garch) {
    omega <- Elnsig2 * (1 - sum(phi))
    sigt_loggarch_forecast_short(
      et = zoo::coredata(object@etat),
      sigt = zoo::coredata(object@sigt),
      omega = omega,
      phi = phi,
      psi = psi,
      Elneta2 = Elneta2,
      horizon = n.ahead
    )
  } else {
    coef_inf_garch <- esemifar::farima_to_ma(ar = phi, ma = psi, d = d,
                                     max_i = n + n.ahead - 1)[-1]
    sigt_egarch_forecast_longCpp(
      et = zoo::coredata(object@etat),
      coef_inf = coef_inf_garch,
      kappa = kappa,
      gamma = gamma,
      E_asy = E_asy,
      E_mag = E_mag,
      Elnsig2 = Elnsig2,
      powers = powers,
      modulus = modulus,
      mode = mode,
      horizon = horizon
    )
  }

  list(
    cmeans = mean_fc,
    sigt = sigt_fc
  )

  }
)

methods::setMethod("predict", "fEGarch_fit",
  function(object, n.ahead = 1, ...) {

    fcasts <- fEGarch_predict(object = object, n.ahead = n.ahead, ...)

  }
)

