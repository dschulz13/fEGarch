initial_split <- function(rt, n_test) {

  split_rt <- split_ts(rt, n_test)
  rt_core_o <- zoo::coredata(split_rt$train)

  scale_const <- sd(rt_core_o)         # for rescaling later
  rt_core <- rt_core_o / scale_const   # now has sample variance 1

  list(
    scale_const = scale_const,
    rt_core = rt_core,
    rt_core_o = rt_core_o
  )

}

adj_sigt_fun <- function(meanspec) {

  # p <- orders[[1]]
  # q <- orders[[2]]
  #
  # cond_d <- match.arg(cond_dist)
  #
  # # ARMA orders
  order_arma <- orders(meanspec)
  p_ar <- order_arma[[1]]
  q_ma <- order_arma[[2]]

  # pdf_fun1 <- lookup_table$fun1[[cond_d]]
  # pdf_fun2 <- lookup_table$fun2[[cond_d]]
  #
  # incl_mean <- meanspec@include_mean
  #
  # skew_par <- lookup_table[["skew_par"]][[cond_d]]
  # extra_par <- lookup_table[["extra_par"]][[cond_d]]

  # Long-memory setting (for ARMA-part)
  lm_arma <- long_memo(meanspec)
  # Assign ARMA / FARIMA function
  arma_fun <- if (lm_arma) {
    grab_fun <- farima_grab_pars
    # D_par_name <- "D"
    # D_start <- 0.25 #mean(Drange)
    # D_low <- Drange[[1]]
    # D_up <- Drange[[2]]
    farima_fit
  } else {
    grab_fun <- arma_grab_pars
    # D_par_name <- character(0)
    # D_start <- D_low <- D_up <- numeric(0)
    arma_fit
  }


  # n <- length(rt_core)
  # if (Prange[[1]] < 1) {Prange[[1]] <- 1}
  # trunc_i <- trunc
  # if (is.character(trunc) && trunc == "none") {
  #   trunc <- n - 1
  # } else if (is.numeric(trunc) && trunc >= n) {
  #   trunc <- n - 1
  #   message("trunc reduced to n - 1")
  # }
  sigt_fun <- sigt_gjrgarch_R
  # Adjust (if ARMA relevant)
  pg0 <- p_ar > 0
  qg0 <- q_ma > 0
  check <- pg0 || qg0 || lm_arma
  n_arma_pars <- incl_mean + p_ar + q_ma + lm_arma

  sigt_fun_adj <- if (check) {
    function(theta, x, p, q, p_ar, q_ma,
         incl_mean, extra_par, skew_par, dfun, sigd_init, etransf_init, trunc) {

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
      # ... and calculate the residuals
      res <- x - fitted_vs

      # Remove the first few elements from theta (those for the parameters
      # relevant for ARMA / FARIMA) ...
      theta2 <- theta[-(1:n_arma_pars)]
      #... and calculate conditional standard deviations from the residuals
      # following the specified GARCH-type
      out <- sigt_fun(theta2, res, p, q, p_ar, q_ma, FALSE, extra_par, skew_par, sigd_init, etransf_init, trunc)
      out$cmeans <- fitted_vs
      out
    }

  } else {
    sigt_fun             # ...else keep function for EGARCH-type
  }

  sigt_fun_adj

}
