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

nllhood_calc <- function(input_list, x, dfun) {

  sigt <- input_list$sigt

  cmeans <- input_list$cmeans
  skew <- input_list$skew
  shape <- input_list$shape

  -sum(log(dfun(x, cmeans, sigt, shape, skew)))

}

crit_calc <- function(pars, rt, llhood) {

  n <- length(rt)
  k <- length(pars)
  m2llhood <- -2 * llhood

  list(
    aic = (m2llhood + 2 * k) / n,
    bic = (m2llhood + k * log(n)) / n
  )

}

# Function to adjust the calculations for the conditional means and conditional
# standard deviations for an ARMA or FARIMA extended volatility model;
# returns the base function for each model, if no additional ARMA / FARIMA
# settings are given
adj_sigt_fun <- function(meanspec,
            model_type = c("garch", "figarch", "aparch", "fiaparch", "gjrgarch",
                      "figjrgarch", "loggarch", "filoggarch", "egarch", "fiegarch")) {

  model_type <- match.arg(model_type)

  # # ARMA orders
  order_arma <- orders(meanspec)
  p_ar <- order_arma[[1]]
  q_ma <- order_arma[[2]]

  incl_mean <- meanspec@include_mean

  # Long-memory setting (for ARMA-part)
  lm_arma <- long_memo(meanspec)
  # Assign ARMA / FARIMA function
  arma_fun <- if (lm_arma) {
    grab_fun <- farima_grab_pars
    farima_fit
  } else {
    grab_fun <- arma_grab_pars
    arma_fit
  }

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

  # Adjust (if ARMA relevant)
  pg0 <- p_ar > 0
  qg0 <- q_ma > 0
  check <- pg0 || qg0 || lm_arma
  n_arma_pars <- incl_mean + p_ar + q_ma + lm_arma

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
