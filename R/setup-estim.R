# Split series into training and test set
# and keep and omit potential time series formatting
initial_split <- function(rt, n_test) {

  split_rt <- split_ts(rt, n_test)
  rt_core_o <- zoo::coredata(split_rt$train)

  scale_const <- sd(rt_core_o)         # for rescaling later
  rt_core <- rt_core_o / scale_const   # now has sample variance 1

  list(
    scale_const = scale_const,
    rt_core = rt_core,
    rt_core_o = rt_core_o,
    test_obs = split_rt$test,
    train_obs = split_rt$train
  )

}

# Compute the negative log-likelihood from
# estimated input series, observations and
# the density function
nllhood_calc <- function(input_list, x, dfun) {

  sigt <- input_list$sigt

  cmeans <- input_list$cmeans
  skew <- input_list$skew
  shape <- input_list$shape

  -sum(log(dfun(x, cmeans, sigt, shape, skew)))

}

# Calculate AIC and BIC from fitted parameter vector length,
# observation number and log-likelihood value
crit_calc <- function(pars, rt, llhood) {

  n <- length(rt)
  k <- length(pars)
  m2llhood <- -2 * llhood

  list(
    aic = (m2llhood + 2 * k) / n,
    bic = (m2llhood + k * log(n)) / n
  )

}

identify_arma_settings <- function(meanspec, Drange) {

  # ARMA orders
  order_arma <- orders(meanspec)
  p_ar <- order_arma[[1]]
  q_ma <- order_arma[[2]]

  # Include unconditional mean parameter?
  incl_mean <- meanspec@include_mean

  # Long-memory setting (for ARMA-part)
  lm_arma <- long_memo(meanspec)

  # Assign ARMA / FARIMA function
  # and parameter grabbing function for ARMA / FARIMA
  arma_fun <- if (lm_arma) {
    grab_fun <- farima_grab_pars
    D_par_name <- "D"
    D_start <- 0.25
    D_low <- Drange[[1]]
    D_up <- Drange[[2]]
    farima_fit
  } else {
    grab_fun <- arma_grab_pars
    D_par_name <- character(0)
    D_start <- D_low <- D_up <- numeric(0)
    arma_fit
  }

  # Adjust (if ARMA relevant)
  pg0 <- p_ar > 0
  qg0 <- q_ma > 0
  # Check if adjustment of sigt-function required due to
  # conditional mean part
  mean_adj_required <- pg0 || qg0 || lm_arma
  # Number of parameters for the mean specification
  n_arma_pars <- incl_mean + p_ar + q_ma + lm_arma

  # Positions of AR-parameters in parameter vector
  low_ar <- 1 + incl_mean
  up_ar <- incl_mean + p_ar

  # Return all
  list(
    p_ar = p_ar,
    q_ma = q_ma,
    incl_mean = incl_mean,
    lm_arma = lm_arma,
    arma_grab_fun = grab_fun,
    cmean_fun = arma_fun,
    pg0 = pg0,
    qg0 = qg0,
    mean_adj_required = mean_adj_required,
    n_arma_pars = n_arma_pars,
    D_par_name = D_par_name,
    D_start = D_start,
    D_low = D_low,
    D_up = D_up,
    low_ar = low_ar,
    up_ar = up_ar,
    ar_start = rep(0.05 / p_ar, p_ar),
    ma_start = rep(0.05 / q_ma, q_ma),
    ar_lb = rep(-2, p_ar),
    ma_lb = rep(-2, q_ma),
    ar_ub = rep(2, p_ar),
    ma_ub = rep(2, q_ma)
  )

}

# Function to adjust the calculations for the conditional means and conditional
# standard deviations for an ARMA or FARIMA extended volatility model;
# returns the base function for each model, if no additional ARMA / FARIMA
# settings are given
adj_sigt_fun <- function(p_ar, q_ma, pg0, qg0, lm_arma, incl_mean, n_arma_pars,
                         grab_fun, arma_fun,
            model_type = c("garch", "figarch", "aparch", "fiaparch", "gjrgarch",
                      "figjrgarch", "loggarch", "filoggarch", "egarch", "fiegarch")) {

  check <- pg0 || qg0 || lm_arma

  model_type <- match.arg(model_type)

  # Select the correct base function for the model type
  sigt_fun <- switch(
    model_type,
    "garch" = sigt_garch_R,
    "figarch" = sigt_figarch_R,
    "aparch" = sigt_aparch_R,
    "fiaparch" = sigt_fiaparch_R,
    "gjrgarch" = sigt_gjrgarch_R,
    "figjrgarch" = sigt_figjrgarch_R,
    "loggarch" = sigt_loggarch_short_R,
    "filoggarch" = sigt_loggarch_long_R,
    "egarch" = sigt_egarch_short_R,
    "fiegarch" = sigt_egarch_long_R
  )

  sigt_fun_adj <- if (check) {
    function(theta, x, p, q, p_ar, q_ma,
         incl_mean, extra_par, skew_par, trunc, ...) {

      # Collect input arguments (including also the ellipsis ...)
      # as a named list
      args <- mget(names(formals()), environment())
      args[["..."]] <- NULL
      dot_args <- list(...)
      all_args <- c(args, dot_args)
      rm(dot_args)
      rm(args)

      # Obtain ARMA / FARIMA parameters from theta
      arma_pars <- grab_fun(
        theta = theta,
        incl_mean = incl_mean,
        pg0 = pg0,
        qg0 = qg0,
        p = p_ar,
        q = q_ma
      )
      mu <- arma_pars$mu
      ar <- arma_pars$ar
      ma <- arma_pars$ma

      coef_inf <- if (lm_arma) {
        coefs <- ar_infty(ar = ar, ma = ma, d = arma_pars$D, max_i = trunc)
        coefs[[1]] <- 0
        coefs
      } else {
        0
      }

      # Calculate fitted values following these parameters...
      fitted_vs <- arma_fun(x = x, mu = mu, ma = ma, ar = ar, coef_inf = coef_inf)

      # ... adjust input arguments from intermediate results ...
      all_args[["theta"]] <- theta[-(1:n_arma_pars)]     # Remove the first few elements from theta (those for the parameters relevant for ARMA / FARIMA)
      all_args[["x"]] <- x - fitted_vs                   # Compute conditional mean adjusted values
      all_args[["incl_mean"]] <- FALSE                   # Remove unconditional mean for further steps

      #... and calculate conditional standard deviations from the residuals
      # following the specified GARCH-type
      out <- do.call(sigt_fun, args = all_args)
      out$cmeans <- fitted_vs
      out
    }

  } else {
    sigt_fun             # ...else keep function for volatility only
  }

  sigt_fun_adj

}

# Function to grab conditional distribution parameter restrictions from
# internally saved tables
cond_d_par_restr <- function(cond_d) {

  vals <- lookup_table[["vals"]][[cond_d]]

  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  list(
    sval = vals[[1]], # starting value
    lval = vals[[2]], # lower value
    uval = vals[[2]]  # upper value
  )
}

set_start_LB_UB <- function(
    mean_val,
    ar_start, ma_start, D_start,
    omega_val,
    phi_start, psi_start,
    add_pars_start,
    d_start,
    sval, skew_start,
    mean_low,
    ar_low, ma_low, D_low,
    omega_low,
    phi_low,
    psi_low,
    add_low,
    d_low,
    lval, skew_low,
    mean_up,
    ar_up, ma_up, D_up,
    omega_up,
    phi_up, psi_up, add_up, d_up,
    uval, skew_up,
    start_pars, LB, UB
) {

  if (is.null(start_pars)) {
    start_pars <- c(
      mean_val,
      ar_start,
      ma_start,
      D_start,
      omega_val,
      phi_start,
      psi_start,
      add_pars_start,
      d_start,
      sval, skew_start
    )
  }

  if (is.null(LB)) {
    LB <- c(
      mean_low,
      ar_low,
      ma_low,
      D_low,
      omega_low,
      phi_low,
      psi_low,
      add_low,
      d_low,
      lval, skew_low
    )
  }

  if (is.null(UB)) {
    UB <- c(
      mean_up,
      ar_up,
      ma_up,
      D_up,
      omega_up,
      phi_up,
      psi_up,
      add_up,
      d_up,
      uval, skew_up
    )
  }

  list(
    start_pars = start_pars,
    LB = LB,
    UB = UB
  )
}
