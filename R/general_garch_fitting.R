general_garch_fitting <- function(rt, model_type, orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                        drange = c(0, 1),
                        meanspec = mean_spec(), Drange = c(0, 1),
                        nonparspec = locpol_spec(), use_nonpar = FALSE,
                        n_test = 0,
                        start_pars = NULL, LB = NULL, UB = NULL, control = list(),
                        control_nonpar = list(), mean_after_nonpar = FALSE,
                        parallel = TRUE, ncores = max(1, future::availableCores() - 1),
                        trunc = "none", presample = 50, Prange = c(1, 5), skip_vcov = FALSE) {


  p <- orders[[1]]
  q <- orders[[2]]

  cond_d <- match.arg(cond_dist)
  lm <- (model_type %in% c("figarch", "figjrgarch", "fiaparch", "fitgarch", "afiaparch", "afigjrgarch", "afitgarch", "fingarch", "afingarch"))


  if (lm && !all(orders == c(1, 1))) {
    orders <- c(1, 1)
    message("Only orders p=1 with q=1 are currently supported for FIGARCH, FITGARCH, FIGJR-GARCH, FIAPARCH, FINGARCH, AFIAPARCH, AFIGJR-GARCH, AFITGARCH and AFINGARCH models.")
  }


  if (use_nonpar) {

    nonpar_out <- run_nonpar(nonparspec = nonparspec,
               control_nonpar = control_nonpar, rt = rt, n_test = n_test,
               mean_after_nonpar = mean_after_nonpar, lm = lm)
    rt <- nonpar_out$rt
    meanspec <- nonpar_out$meanspec
    nonpar_result <- nonpar_out$nonpar_result
    n_test <- nonpar_out$n_test

  }

  split_data <- initial_split(rt, n_test)

  n <- length(split_data$rt_core)
  if (Prange[[1]] < 1) {Prange[[1]] <- 1}
  trunc_i <- trunc

  # Settings for ARMA / FARIMA part in the mean; makes objects available in this
  # environment
  arma_sets <- identify_arma_settings(meanspec = meanspec, Drange = Drange)

  if ((lm || arma_sets$lm_arma) && is.character(trunc) && trunc == "none") {
    trunc <- n + presample - 1
  } else if ((lm || arma_sets$lm_arma) && is.numeric(trunc) && trunc >= n + presample) {
    trunc <- n + presample - 1
    message("trunc reduced to n + ", presample, " - 1")
  }

  dfun1 <- fun1_selector(cond_d)
  dfun2 <- fun2_selector(cond_d)

  skew_par <- lookup_table[["skew_par"]][[cond_d]]
  extra_par <- lookup_table[["extra_par"]][[cond_d]]

  sigt_fun_adj <- adj_sigt_fun(
    p_ar = arma_sets$p_ar, q_ma = arma_sets$q_ma,
    pg0 = arma_sets$pg0, qg0 = arma_sets$qg0,
    lm_arma = arma_sets$lm_arma, incl_mean = arma_sets$incl_mean,
    n_arma_pars = arma_sets$n_arma_pars,
    grab_fun = arma_sets$arma_grab_fun,
    arma_fun = arma_sets$cmean_fun,
    model_type = model_type
  )

  creator_fun <- switch(
    model_type,
    "garch" = goal_fun_creator_garch,
    "gjrgarch" = goal_fun_creator_gjrgarch,
    "tgarch" = goal_fun_creator_tgarch,
    "aparch" = goal_fun_creator_aparch,
    "figarch" = goal_fun_creator_figarch,
    "figjrgarch" = goal_fun_creator_figjrgarch,
    "fitgarch" = goal_fun_creator_fitgarch,
    "fiaparch" = goal_fun_creator_fiaparch,
    "afiaparch" = goal_fun_creator_afiaparch,
    "afigjrgarch" = goal_fun_creator_afigjrgarch,
    "afitgarch" = goal_fun_creator_afitgarch,
    "ngarch" = goal_fun_creator_ngarch,
    "afingarch" = goal_fun_creator_afingarch,
    "fingarch" = goal_fun_creator_fingarch
  )

  base_args <- list(
    p = p,
    q = q,
    p_ar = arma_sets$p_ar,
    q_ma = arma_sets$q_ma,
    incl_mean = arma_sets$incl_mean,
    extra_par = extra_par,
    skew_par = skew_par,
    sig_fun = sigt_fun_adj,
    trunc = trunc,
    presample = presample
  )

  mean_req_1 <- (!arma_sets$lm_arma && (meanspec@orders[[1]] > 0 || meanspec@orders[[2]]))
  mean_req_2 <- arma_sets$lm_arma
  mean_res <- if (mean_req_1) {   # run an initial ARMA to get better pre-sample values for the volatility step
    init_arma <- suppressWarnings(stats::arima(split_data$rt_core, order = c(meanspec@orders[[1]], 0, meanspec@orders[[2]]), include.mean = meanspec@include_mean, method = "ML"))
    init_arma$residuals
  } else if (mean_req_2) {
    init_farima <- tryCatch({
      suppressWarnings(prelim_farima(x = split_data$rt_core, p = meanspec@orders[[1]], q = meanspec@orders[[2]], include_mean = meanspec@include_mean, presample = presample, trunc = trunc, Drange = Drange))
    },
    error = {
      suppressWarnings(stats::arima(split_data$rt_core, order = c(meanspec@orders[[1]], 0, meanspec@orders[[2]]), include.mean = meanspec@include_mean, method = "ML"))
    })
    init_farima$residuals
  } else {
    split_data$rt_core
  }

  extra_args <- switch(
    model_type,
    "garch" = list(sig2_init = stats::var(mean_res)),
    "gjrgarch" = {
      delta0 <- 2    # use delta0 = 2
      gamma0 <- 0  # use gamma0 = 0   as suitable initial guesses
      sigd_init <- stats::sd(mean_res)^delta0
      rt_dm <-  mean_res - mean(mean_res)
      etransf_init <- mean((abs(rt_dm ) - gamma0 *  rt_dm)^delta0)
      list(sigd_init = sigd_init, etransf_init = etransf_init)
    },
    "tgarch" = {
      delta0 <- 1    # use delta0 = 1
      gamma0 <- 0  # use gamma0 = 0   as suitable initial guesses
      sigd_init <- stats::sd(mean_res)^delta0
      rt_dm <- mean_res - mean(mean_res)
      etransf_init <- mean((abs(rt_dm) - gamma0 * rt_dm)^delta0)
      list(sigd_init = sigd_init, etransf_init = etransf_init)
    },
    "ngarch" = {
      delta0 <- 2    # use delta0 = 2
      gamma0 <- 0  # use gamma0 = 0   as suitable initial guesses
      sigd_init <- stats::sd(mean_res)^delta0
      rt_dm <- mean_res - mean(mean_res)
      etransf_init <- mean((abs(rt_dm) - gamma0 * rt_dm)^delta0)
      list(sigd_init = sigd_init, etransf_init = etransf_init)
    },
    "aparch" = {
      delta0 <- 2    # use delta0 = 2
      gamma0 <- 0  # use gamma0 = 0   as suitable initial guesses
      sigd_init <- stats::sd(mean_res)^delta0
      rt_dm <- mean_res - mean(mean_res)
      etransf_init <- mean((abs(rt_dm) - gamma0 * rt_dm)^delta0)
      list(sigd_init = sigd_init, etransf_init = etransf_init)
    },
    "figarch" = {
      list(presample_val = var(mean_res))
    },
    "figjrgarch" = {
      est_s <- suppressWarnings(gjrgarch(
        rt = rt,
        orders = orders,
        cond_dist = cond_dist,
        meanspec = meanspec,
        use_nonpar = use_nonpar,
        nonparspec = nonparspec,
        parallel = parallel,
        ncores = ncores,
        Prange = Prange,
        n_test = n_test,
        skip_vcov = TRUE
      ))
      delta0 <- 2
      gamma0 <- est_s@pars[["gamma1"]]
      rt_dm <- mean_res - mean(mean_res)
      list(presample_val = mean((abs(rt_dm) - gamma0 * rt_dm)^delta0))
    },
    "fingarch" = {
      est_s <- suppressWarnings(ngarch(
        rt = rt,
        orders = orders,
        cond_dist = cond_dist,
        meanspec = meanspec,
        use_nonpar = use_nonpar,
        nonparspec = nonparspec,
        parallel = parallel,
        ncores = ncores,
        Prange = Prange,
        n_test = n_test,
        skip_vcov = TRUE
      ))
      delta0 <- est_s@pars[["delta"]]
      gamma0 <- 0  # fixed for NGARCH-type
      rt_dm <- mean_res - mean(mean_res)
      list(presample_val = mean((abs(rt_dm) - gamma0 * rt_dm)^delta0))
    },
    "fitgarch" = {
      est_s <- suppressWarnings(tgarch(
        rt = rt,
        orders = orders,
        cond_dist = cond_dist,
        meanspec = meanspec,
        use_nonpar = use_nonpar,
        nonparspec = nonparspec,
        parallel = parallel,
        ncores = ncores,
        Prange = Prange,
        n_test = n_test,
        skip_vcov = TRUE
      ))
      delta0 <- 1
      gamma0 <- est_s@pars[["gamma1"]]
      rt_dm <- mean_res - mean(mean_res)
      list(presample_val = mean((abs(rt_dm) - gamma0 * rt_dm)^delta0))
    },
    "fiaparch" = {
      est_s <- suppressWarnings(aparch(
        rt = rt,
        orders = orders,
        cond_dist = cond_dist,
        meanspec = meanspec,
        use_nonpar = use_nonpar,
        nonparspec = nonparspec,
        parallel = parallel,
        ncores = ncores,
        Prange = Prange,
        n_test = n_test,
        skip_vcov = TRUE
      ))
      delta0 <- est_s@pars[["delta"]]
      gamma0 <- est_s@pars[["gamma1"]]
      rt_dm <- mean_res - mean(mean_res)
      list(presample_val = mean((abs(rt_dm) - gamma0 * rt_dm)^delta0))
    },
    "afiaparch" = {
      # Initial APARCH fit to obtain guesses of relevant parameters
      est_s <- suppressWarnings(aparch(
        rt = rt,
        orders = orders,
        cond_dist = cond_dist,
        meanspec = meanspec,
        use_nonpar = use_nonpar,
        nonparspec = nonparspec,
        parallel = parallel,
        ncores = ncores,
        Prange = Prange,
        n_test = n_test,
        skip_vcov = TRUE
      ))
      # For AFIAPARCH we need all distribution parameters to obtain
      # presample values due to the rescaling
      delta0 <- est_s@pars[["delta"]]
      gamma0 <- est_s@pars[["gamma1"]]
      shape_nam <- c("norm" = "NULL", "snorm" = "NULL", "std" = "df", "sstd" = "df",
                     "ged" = "shape", "sged" = "shape" , "ald" = "P", "sald" = "P")[[cond_d]]
      shape0 <- ifelse(cond_d %in% c("std", "sstd", "ged", "sged", "ald", "sald"),
                       est_s@pars[[shape_nam]], 0)
      skew0 <- ifelse(cond_d %in% c("snorm", "sstd", "sged", "sald"), est_s@pars[["skew"]], 0)
      C0 <- comp_C_zd(shape = shape0, skew = skew0, gamma = gamma0, delta = delta0, dens_fun = dfun1)
      rt_dm <- afiaparch_adjust_returns(mean_res - mean(mean_res), gamma = gamma0,
                                        delta = delta0, Czd = C0)
      # Additionally, also pass the innovation density to the goal function
      # for AFIAPARCH, as it is required for the computation of the constant
      list(presample_val = mean(rt_dm))
    },
    "afitgarch" = {
      # Initial APARCH fit to obtain guesses of relevant parameters
      est_s <- suppressWarnings(tgarch(
        rt = rt,
        orders = orders,
        cond_dist = cond_dist,
        meanspec = meanspec,
        use_nonpar = use_nonpar,
        nonparspec = nonparspec,
        parallel = parallel,
        ncores = ncores,
        Prange = Prange,
        n_test = n_test,
        skip_vcov = TRUE
      ))
      # For AFITGARCH we need all distribution parameters to obtain
      # presample values due to the rescaling
      delta0 <- 1  # fixed
      gamma0 <- est_s@pars[["gamma1"]]
      shape_nam <- c("norm" = "NULL", "snorm" = "NULL", "std" = "df", "sstd" = "df",
                     "ged" = "shape", "sged" = "shape" , "ald" = "P", "sald" = "P")[[cond_d]]
      shape0 <- ifelse(cond_d %in% c("std", "sstd", "ged", "sged", "ald", "sald"),
                       est_s@pars[[shape_nam]], 0)
      skew0 <- ifelse(cond_d %in% c("snorm", "sstd", "sged", "sald"), est_s@pars[["skew"]], 0)
      C0 <- comp_C_zd(shape = shape0, skew = skew0, gamma = gamma0, delta = delta0, dens_fun = dfun1)
      rt_dm <- afiaparch_adjust_returns(mean_res - mean(mean_res), gamma = gamma0,
                                        delta = delta0, Czd = C0)
      # Additionally, also pass the innovation density to the goal function
      # for AFIAPARCH, as it is required for the computation of the constant
      list(presample_val = mean(rt_dm))
    },
    "afingarch" = {
      # Initial APARCH fit to obtain guesses of relevant parameters
      est_s <- suppressWarnings(ngarch(
        rt = rt,
        orders = orders,
        cond_dist = cond_dist,
        meanspec = meanspec,
        use_nonpar = use_nonpar,
        nonparspec = nonparspec,
        parallel = parallel,
        ncores = ncores,
        Prange = Prange,
        n_test = n_test,
        skip_vcov = TRUE
      ))
      delta0 <- est_s@pars[["delta"]]
      gamma0 <- 0    # fixed
      shape_nam <- c("norm" = "NULL", "snorm" = "NULL", "std" = "df", "sstd" = "df",
                     "ged" = "shape", "sged" = "shape" , "ald" = "P", "sald" = "P")[[cond_d]]
      shape0 <- ifelse(cond_d %in% c("std", "sstd", "ged", "sged", "ald", "sald"),
                       est_s@pars[[shape_nam]], 0)
      skew0 <- ifelse(cond_d %in% c("snorm", "sstd", "sged", "sald"), est_s@pars[["skew"]], 0)
      C0 <- comp_C_zd(shape = shape0, skew = skew0, gamma = gamma0, delta = delta0, dens_fun = dfun1)
      rt_dm <- afiaparch_adjust_returns(mean_res - mean(mean_res), gamma = gamma0,
                                        delta = delta0, Czd = C0)
      # Additionally, also pass the innovation density to the goal function
      # for AFIAPARCH, as it is required for the computation of the constant
      list(presample_val = mean(rt_dm))
    },
    "afigjrgarch" = {
      # Initial APARCH fit to obtain guesses of relevant parameters
      est_s <- suppressWarnings(gjrgarch(
        rt = rt,
        orders = orders,
        cond_dist = cond_dist,
        meanspec = meanspec,
        use_nonpar = use_nonpar,
        nonparspec = nonparspec,
        parallel = parallel,
        ncores = ncores,
        Prange = Prange,
        n_test = n_test,
        skip_vcov = TRUE
      ))
      # For AFIGJRGARCH we need all distribution parameters to obtain
      # presample values due to the rescaling
      delta0 <- 2   # Fixed for GJR-GARCH type
      gamma0 <- est_s@pars[["gamma1"]]
      shape_nam <- c("norm" = "NULL", "snorm" = "NULL", "std" = "df", "sstd" = "df",
                     "ged" = "shape", "sged" = "shape" , "ald" = "P", "sald" = "P")[[cond_d]]
      shape0 <- ifelse(cond_d %in% c("std", "sstd", "ged", "sged", "ald", "sald"),
                       est_s@pars[[shape_nam]], 0)
      skew0 <- ifelse(cond_d %in% c("snorm", "sstd", "sged", "sald"), est_s@pars[["skew"]], 0)
      C0 <- comp_C_zd(shape = shape0, skew = skew0, gamma = gamma0, delta = delta0, dens_fun = dfun1)
      rt_dm <- afiaparch_adjust_returns(mean_res - mean(mean_res), gamma = gamma0,
                                        delta = delta0, Czd = C0)
      # Additionally, also pass the innovation density to the goal function
      # for AFIAPARCH, as it is required for the computation of the constant
      list(presample_val = mean(rt_dm))
    }
  )

  adjusted_model_check <- model_type %in% c("afiaparch", "afingarch", "afigjrgarch", "afitgarch")

  goal_fun_s <- if (adjusted_model_check) {
      function(theta, rt, dfun1, dfun2) {

        do.call(
          what = creator_fun,
          args = c(
            list(x = rt, theta = theta, pdf_fun2 = dfun2, dens_fun = dfun1),
            base_args,
            extra_args
          )
        )

      }
    } else {
      function(theta, rt, dfun1, dfun2) {

        do.call(
          what = creator_fun,
          args = c(
            list(x = rt, theta = theta, pdf_fun2 = dfun2),
            base_args,
            extra_args
          )
        )

      }
    }


  # Get cond. distribution restrictions and starting values
  # as objects sval, lval and uval in the environment
  cond_d_vals <- cond_d_par_restr(cond_d)

  # Make further starting parameter values and restrictions available
  # using internally saved table function
  pars_and_restr <- start_pars_and_restr[[model_type]](drange, p, q)

  if (is.null(start_pars) || is.null(LB) || is.null(UB)) {
    mean_rt <- list(NULL, mean(split_data$rt_core))[[arma_sets$incl_mean + 1]]
    min_rt <- list(NULL, suppressWarnings(0.2 * min(split_data$rt_core - mean_rt) + mean_rt))[[arma_sets$incl_mean + 1]]
    max_rt <- list(NULL, suppressWarnings(0.2 * max(split_data$rt_core - mean_rt) + mean_rt))[[arma_sets$incl_mean + 1]]
    skew_start <- list(NULL, 0.98)[[skew_par + 1]]
    skew_low <- list(NULL, 1e-15)[[skew_par + 1]]
    skew_up <- list(NULL, Inf)[[skew_par + 1]]

    start_LB_UB <- set_start_LB_UB(
      mean_val = mean_rt,
      ar_start = arma_sets$ar_start,
      ma_start = arma_sets$ma_start,
      D_start = arma_sets$D_start,
      omega_val = pars_and_restr$omega_start,
      phi_start = pars_and_restr$phi_start,
      psi_start = pars_and_restr$psi_start,
      add_pars_start = pars_and_restr$add_pars_start,
      d_start = pars_and_restr$d_spar,
      sval = cond_d_vals$sval,
      skew_start = skew_start,
      mean_low = min_rt,
      ar_low = arma_sets$ar_lb,
      ma_low = arma_sets$ma_lb,
      D_low = arma_sets$D_low,
      omega_low = pars_and_restr$omega_low,
      phi_low = pars_and_restr$phi_low,
      psi_low = pars_and_restr$psi_low,
      add_low = pars_and_restr$add_low,
      d_low = pars_and_restr$d_low,
      lval = cond_d_vals$lval, skew_low = skew_low,
      mean_up = max_rt,
      ar_up = arma_sets$ar_ub, ma_up = arma_sets$ma_ub,
      D_up = arma_sets$D_up,
      omega_up = pars_and_restr$omega_up,
      phi_up = pars_and_restr$phi_up, psi_up = pars_and_restr$psi_up,
      add_up = pars_and_restr$add_up,
      d_up = pars_and_restr$d_up,
      uval = cond_d_vals$uval, skew_up = skew_up,
      start_pars = start_pars, LB = LB, UB = UB
    )


  }

  constr_mean <- if (arma_sets$pg0) {

    ineqLB_mean <- rep(1e-6, arma_sets$p_ar)
    ineqUB_mean <- rep(Inf, arma_sets$p_ar)

    # Stationarity constraint
    function(theta) {
      ar <- c(1, -theta[arma_sets$ar_pos_low:arma_sets$ar_pos_up])
      abs(polyroot(ar)) - 1
    }

  } else {

    ineqLB_mean <- NULL
    ineqUB_mean <- NULL

    # Not required if no AR-part
    function(theta) {
      NULL
    }

  }

  constr_vol <- if (lm) {

    low_phi <- 2 + arma_sets$n_arma_pars
    up_phi <- low_phi
    low_beta <- low_phi + 1
    up_beta <- low_beta

    add_s <- switch(
      model_type,
      "figarch" = 1,
      "figjrgarch" = 2,
      "afigjrgarch" = 2,
      "fitgarch" = 2,
      "afitgarch" = 2,
      "afingarch" = 2,
      "fingarch" = 2,
      "fiaparch" = 3,
      "afiaparch" = 3
    )

    d_idx <- up_beta + add_s


    ineqLB_vol <- c(rep(0, 50), 1e-6)
    ineqUB_vol <- rep(Inf, 51)

    function(theta) {

      phi <- theta[low_phi:up_phi]
      beta <- theta[low_beta:up_beta]
      d <- theta[[d_idx]]

      coef_inf <- ar_infty(ar = phi, ma = -beta, d = d, max_i = 50)[-1]

      c2 <- abs(polyroot(c(1, -beta))) - 1
      c(coef_inf, c2)

    }
  } else {

      if (model_type == "garch") {
        low_phi <- 2 + arma_sets$n_arma_pars
        up_phi <- low_phi - 1 + p
        low_beta <- low_phi + 1
        up_beta <- low_beta - 1 + q

        ineqLB_vol <- c(0)
        ineqUB_vol <- 1 - 1e-6

        function(theta) {
          phi <- theta[low_phi:up_phi]
          beta <- theta[low_beta:up_beta]

          sum(phi) + sum(beta)
        }
      } else if (model_type %in% c("gjrgarch", "tgarch", "aparch", "ngarch")) {
        ineqLB_vol <- NULL
        ineqUB_vol <- NULL

        function(theta) {
          NULL
        }
      }

    }

  # Combined constraint bounds
  ineqLB <- c(ineqLB_mean, ineqLB_vol)
  ineqUB <- c(ineqUB_mean, ineqUB_vol)

  # Combined constraint function
  constr_fun <- if (lm || arma_sets$pg0) {
    function(theta, rt) {

      constr1 <- constr_mean(theta)
      constr2 <- constr_vol(theta)

      c(constr1, constr2)

    }

  } else {

    NULL

  }

  message_occurred <- FALSE
  if (is.null(control$trace)) {control$trace <- FALSE}

  if (cond_d %in% c("ald", "sald")) {

    # Makes "P_sel" and "result" available
    out <- withCallingHandlers({ald_fit(
      Prange = Prange,
      parallel = parallel,
      ncores = ncores,
      dfun1 = dfun1,
      dfun2 = dfun2,
      goal_fun_s = goal_fun_s,
      start_pars = start_LB_UB$start_pars,
      LB = start_LB_UB$LB, UB = start_LB_UB$UB,
      ineqfun = constr_fun,
      ineqLB = ineqLB,
      ineqUB = ineqUB,
      rt_core = split_data$rt_core,
      control = control
    )}, message = function(m1) {
      if (grepl("Czd is NA or NaN on some iterations", conditionMessage(m1), fixed = TRUE)) {
        message_occurred <<- TRUE
        invokeRestart("muffleMessage")
      }
    })

    result <- out$result
    P_sel <- out$P_sel

  } else {

    P_sel <- NULL

    goal_fun <- if (model_type %in% c("afiaparch", "afigjrgarch", "afitgarch", "afingarch")) {
      function(theta, rt) {
        goal_fun_s(theta = theta, rt = rt, dfun1 = dfun1, dfun2 = dfun2)
      }
    } else {
      function(theta, rt) {
        goal_fun_s(theta = theta, rt = rt, dfun1 = NULL, dfun2 = dfun2)
      }
    }

    result <- withCallingHandlers({tryCatch(
      expr = {
        suppressWarnings(Rsolnp::solnp(
          pars = start_LB_UB$start_pars,
          fun = goal_fun,
          LB = start_LB_UB$LB,
          UB = start_LB_UB$UB,
          control = control,
          ineqfun = constr_fun,
          ineqLB = ineqLB,
          ineqUB = ineqUB,
          rt = split_data$rt_core
        ))
      },
      error = function(e1) {
        stop("Error during optimization. You may want to try different starting parameter and / or solver settings.", call. = FALSE)
      })},
      message = function(m1) {
        if (grepl("Czd is NA or NaN on some iterations", conditionMessage(m1), fixed = TRUE)) {
          message_occurred <<- TRUE
          invokeRestart("muffleMessage")
        }
      })

    # result <- Rsolnp::solnp(
    #       pars = start_pars,
    #       fun = goal_fun,
    #       LB = LB,
    #       UB = UB,
    #       control = control,
    #       ineqfun = constr_fun,
    #       ineqLB = ineqLB,
    #       ineqUB = ineqUB,
    #       rt = rt_core
    #     )

  }

  pf <- parent.frame()
  if (!exists(x = ".debug_mode_fEGarch", envir = pf, inherits = TRUE)) {.debug_mode_fEGarch <- FALSE}
  if (.debug_mode_fEGarch && message_occurred) {warning("Note: Czd is NA or NaN on some iterations.")}
  if (result$convergence %in% c(1, 2)) {
    warning("Convergence failed. You may want to try different starting parameter and / or solver settings.", call. = FALSE)
  }

  pars <- result$pars

  # Compute log-likelihood which is shifted from the obtained negative
  # log-likelihood at the optimum; "-tail(result$values, 1)" is the
  # log-likelihood at the optimum for the scaled data, which is too large by
  # "n * log(scale_const)" in comparison to original data
  llhood <- -tail(result$values, 1) - n * log(split_data$scale_const)

  beta_names <- paste0("beta", 1:q)

  extra_par_name <- switch(
    cond_d,
    "std" = "df",
    "sstd" = "df",
    "ged" = "shape",
    "sged" = "shape",
    "ald" = "P",
    "sald" = "P",
    character(0)
  )


  delta_name <- if (model_type %in% c("aparch", "fiaparch", "afiaparch", "ngarch", "afingarch", "fingarch")) {
    "delta"
  } else {
    character(0)
  }

  gamma_name <- if (model_type %in% c("fiaparch", "figjrgarch", "fitgarch", "afiaparch", "afigjrgarch", "afitgarch")) {
    "gamma"
  } else if (model_type %in% c("aparch", "gjrgarch", "tgarch")) {
    paste0("gamma", 1:p)
  } else {
    character(0)
  }

  par_names <- c(
    list(NULL, "mu")[[arma_sets$incl_mean + 1]],
    arma_sets$names_ar,
    arma_sets$names_ma,
    arma_sets$D_par_name,
    "omega",
    paste0("phi", 1:p),
    beta_names,
    gamma_name,
    delta_name,
    pars_and_restr$d_name,
    extra_par_name,
    list(NULL, "skew")[[skew_par + 1]]
  )

  if (cond_d %in% c("ald", "sald")) {

    # Fix dfun1
    dfun1_P <- function(x, mu, sigt, shape, skew) {
      dfun1(x = x, shape = P_sel, skew = skew)
    }

    # Fix dfun2
    dfun2_P <- function(x, mu, sigt, shape, skew) {
      dfun2(x = x, mu = mu, sigt = sigt, shape = P_sel, skew = skew)
    }

    # Fix the final goal function to optimize over
    # using dfun1_P and dfun2_P
    goal_fun <- function(theta, rt) {
      goal_fun_s(theta = theta, rt = rt, dfun1 = dfun1_P, dfun2 = dfun2_P)
    }

  }

  fargs <- c(list(
      theta = unname(pars),
      x = split_data$rt_core,
      p = p,
      q = q,
      p_ar = arma_sets$p_ar,
      q_ma = arma_sets$q_ma,
      incl_mean = arma_sets$incl_mean,
      extra_par = extra_par,
      skew_par = skew_par,
      trunc = trunc,
      presample = presample
    ), extra_args)

  if (cond_d %in% c("ald", "sald")) {
    fargs[["dens_fun"]] <- dfun1_P
  } else {
    fargs[["dens_fun"]] <- dfun1
  }

  # Compute the conditional standard deviations at the optimum
  final_series <- do.call(
    sigt_fun_adj,
    args = fargs
  )

  delta <- switch(
    model_type,
    "garch" = 2,
    "gjrgarch" = 2,
    "tgarch" = 1,
    "aparch" = pars[[arma_sets$n_arma_pars + 2 * p + q + 2]],
    "figarch" = 2,
    "figjrgarch" = 2,
    "afigjrgarch" = 2,
    "fitgarch" = 1,
    "afitgarch" = 1,
    "fiaparch" = pars[[arma_sets$n_arma_pars + 2 * p + q + 2]],
    "afiaparch" = pars[[arma_sets$n_arma_pars + 2 * p + q + 2]],
    "ngarch" = pars[[arma_sets$n_arma_pars + 2 + p + q]],
    "afingarch" = pars[[arma_sets$n_arma_pars + 2 + p + q]],
    "fingarch" = pars[[arma_sets$n_arma_pars + 2 + p + q]],
  )

  sc_delta <- split_data$scale_const^delta

  type_inp <- c("aparch", "fiaparch")[[lm + 1]]
  if (model_type %in% c("aparch", "ngarch", "fiaparch", "fingarch", "afiaparch", "afingarch")) {
    delta_pos <- which(par_names == "delta")
    delta_inp <- NULL
    adj_fun_delta <- function(theta, delta_pos, scale_const, delta_inp) {
      scale_const^theta[[delta_pos]]
    }
  } else {
    delta_pos <- NULL
    delta_inp <- switch(
      model_type,
      "garch" = 2,
      "figarch" = 2,
      "tgarch" = 1,
      "fitgarch" = 1,
      "afitgarch" = 1,
      "gjrgarch" = 2,
      "figjrgarch" = 2,
      "afigjrgarch" = 2
    )
    adj_fun_delta <- function(theta, delta_pos, scale_const, delta_inp) {
      scale_const^delta_inp
    }
  }

  vcov_output <- compute_vcov(
    goal_fun = goal_fun,
    pars = pars,
    rt_core = split_data$rt_core,
    scale_const = split_data$scale_const,
    incl_mean = arma_sets$incl_mean,
    P_sel = P_sel,
    cond_d = cond_d,
    par_names = par_names,
    sc_delta = sc_delta,
    delta_pos = delta_pos,
    delta_inp = delta_inp,
    adj_fun_delta = adj_fun_delta,
    n_arma_pars = arma_sets$n_arma_pars,
    model_type = type_inp,
    skip_vcov = skip_vcov
  )

  rescale_output <- rescale_estim(
    pars = vcov_output$pars,
    rt_core_o = split_data$rt_core_o,
    sigt = final_series$sigt,
    cmeans = final_series$cmeans,
    scale_const = split_data$scale_const,
    n_arma_pars = arma_sets$n_arma_pars,
    incl_mean = arma_sets$incl_mean,
    model_type = type_inp,
    sc_delta = sc_delta
  )

  names(rescale_output$pars) <- par_names

  # Makes "aic" and "bic" available
  sel_criteria <- crit_calc(pars = rescale_output$pars, rt = split_data$rt_core_o,
                            llhood = llhood)

  # Apply time series formatting (if relevant) to all
  # estimated series
  series_out <- format_applier_ts(
    rt = split_data$train_obs,
    list_of_ts = list(
      "cmeans" = rescale_output$cmeans,
      "sigt" = rescale_output$sigt,
      "etat" = rescale_output$etat
    )
  )

  cmeans <- series_out$cmeans
  sigt <- series_out$sigt
  etat <- series_out$etat

  construct_fun <- switch(
    model_type,
    "garch" = fEGarch_fit_garch,
    "gjrgarch" = fEGarch_fit_gjrgarch,
    "tgarch" = fEGarch_fit_tgarch,
    "aparch" = fEGarch_fit_aparch,
    "figarch" = fEGarch_fit_figarch,
    "figjrgarch" = fEGarch_fit_figjrgarch,
    "fitgarch" = fEGarch_fit_fitgarch,
    "fiaparch" = fEGarch_fit_fiaparch,
    "afiaparch" = fEGarch_fit_afiaparch,
    "afigjrgarch" = fEGarch_fit_afigjrgarch,
    "afitgarch" = fEGarch_fit_afitgarch,
    "ngarch" = fEGarch_fit_ngarch,
    "afingarch" = fEGarch_fit_afingarch,
    "fingarch" = fEGarch_fit_fingarch
  )

  final_Cz <- if (model_type %in% c("afiaparch", "afigjrgarch", "afitgarch", "afingarch")) {
    final_series$Czd
  } else if (model_type %in% c("fiaparch", "figjrgarch", "fitgarch", "fingarch")) {
    shape_o <- if (cond_d %in% c("std", "sstd", "ged", "sged", "ald", "sald")) {rescale_output$pars[[extra_par_name]]} else {0}
    skew_o <- if (cond_d %in% c("snorm", "sstd", "sged", "sald")) {rescale_output$pars[["skew"]]} else {0}
    delta_o <- if (model_type %in% c("fiaparch", "fingarch")) {
      rescale_output$pars[["delta"]]
    } else if (model_type == "fitgarch") {
      1
    } else if (model_type == "figjrgarch") {
      2
    }
    gamma_o <- if (model_type %in% c("fiaparch", "fitgarch", "figjrgarch")) {rescale_output$pars[["gamma"]]} else {0}

    comp_C_zd(shape = shape_o, skew = skew_o, gamma = gamma_o, delta = delta_o, dens_fun = dfun1)
  } else {
    0
  }

  list_out <- if (model_type %in% c("fiaparch", "fitgarch", "figjrgarch", "fingarch", "afiaparch", "afitgarch", "afigjrgarch", "afingarch")) {

    max_c_len <- if (is.null(trunc) || trunc == "none") {
      presample + length(split_data$train_obs) - 1
    } else {
      min(trunc, presample + length(split_data$train_obs) - 1)
    }
    d_coef <- unname(rescale_output$pars[["d"]])
    phi_coef <- unname(rescale_output$pars[["phi1"]])
    beta_coef <- -unname(rescale_output$pars[["beta1"]])
    coef_inf <- ar_infty(ar = phi_coef, ma = beta_coef, d = d_coef, max_i = max_c_len)[-1]
    C_em <- if (model_type %in% c("fiaparch", "fitgarch", "figjrgarch", "fingarch")) {
      final_Cz * sum(coef_inf)
    } else {
      sum(coef_inf)
    }


    construct_fun(
      pars = rescale_output$pars,
      se = vcov_output$serrors,
      vcov_mat = vcov_output$vcov_mat,
      rt = split_data$train_obs,
      sigt = sigt,
      cmeans = cmeans,
      scale_fun = NULL,
      etat = etat,
      llhood = llhood,
      inf_criteria = c("aic" = sel_criteria$aic, "bic" = sel_criteria$bic),
      cond_dist = cond_d,
      orders = orders,
      long_memo = lm,
      meanspec = meanspec,
      test_obs = split_data$test_obs,
      nonpar_model = NULL,
      trunc = trunc_i,
      Cz = final_Cz,
      Cem = C_em
    )
  } else {
    construct_fun(
      pars = rescale_output$pars,
      se = vcov_output$serrors,
      vcov_mat = vcov_output$vcov_mat,
      rt = split_data$train_obs,
      sigt = sigt,
      cmeans = cmeans,
      scale_fun = NULL,
      etat = etat,
      llhood = llhood,
      inf_criteria = c("aic" = sel_criteria$aic, "bic" = sel_criteria$bic),
      cond_dist = cond_d,
      orders = orders,
      long_memo = lm,
      meanspec = meanspec,
      test_obs = split_data$test_obs,
      nonpar_model = NULL,
      trunc = trunc_i
    )
  }


  if (use_nonpar) {
    mu <- rep(nonpar_result$mu, length(rt))
    list_out@rt <- nonpar_result$rt_train
    format_list <- format_applier_ts(
      rt = nonpar_result$rt_train,
      list_of_ts = list(
        scale_fun = nonpar_result$scale_fun,
        cmeans = mu
      )
    )
    list_out@scale_fun <- format_list$scale_fun
    list_out@cmeans <- format_list$cmeans
    list_out@test_obs <- nonpar_result$test_obs
    list_out@sigt <- list_out@sigt * format_list$scale_fun
    list_out@nonpar_model <- nonpar_result$est_nonpar
  }

  list_out

}
