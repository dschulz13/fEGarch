garch_grab_pars <- function(theta, p, q, incl_mean, extra_par, skew_par) {

  imean <- !incl_mean
  mu <- if (incl_mean) theta[[1]] else 0
  step0 <- 2 - imean
  omega <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q
  phi <- theta[(step0 + 1):step1]
  beta <- theta[(step1 + 1):step2]
  shape <- if (extra_par) theta[[step2 + 1]] else 0
  skew <- if (skew_par) theta[[step2 + 1 + extra_par]] else 0

  list(
    mu = mu,
    omega = omega,
    phi = phi,
    beta = beta,
    shape = shape,
    skew = skew
  )

}

sigt_garch_R <- function(theta, x, p, q, p_ar, q_ma,
              incl_mean, extra_par, skew_par, trunc, ...) {

  list2env(list(...), envir = environment())

  all_pars <- garch_grab_pars(
    theta = theta,
    p = p,
    q = q,
    incl_mean = incl_mean,
    extra_par = extra_par,
    skew_par = skew_par
  )
  mu <- all_pars$mu
  omega <- all_pars$omega
  phi <- all_pars$phi
  beta <- all_pars$beta
  shape <- all_pars$shape
  skew <- all_pars$skew

  sigt <- c(sigt_garch_Cpp(
    x = x,
    mu = mu,
    phi = phi,
    beta = beta,
    omega = omega,
    sig2_init
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape, cmeans = rep(mu, length(sigt)))

}

goal_fun_creator_garch <- function(
    theta,
    x,
    pdf_fun2,
    p,
    q,
    p_ar,
    q_ma,
    incl_mean,
    extra_par,
    skew_par,
    sig_fun,
    sig2_init,
    trunc
) {

  res_list <- sig_fun(
    theta = theta,
    x = x,
    p = p,
    q = q,
    p_ar = p_ar,
    q_ma = q_ma,
    incl_mean = incl_mean,
    extra_par = extra_par,
    skew_par = skew_par,
    trunc = trunc,
    sig2_init = sig2_init)

  nllhood_calc(res_list, x, pdf_fun2)

}



#'GARCH Model Fitting
#'
#'Fit a standard GARCH model under the six
#'most common conditional distributions (and more) to observed data
#'using quasi maximum-likelihood estimation.
#'
#'@param rt the observed series ordered from past to present; can be
#'a numeric vector, a \code{"zoo"} class time series object, or a
#'\code{"ts"} class time series object.
#'@param orders a two-element numeric vector containing the two model
#'orders \eqn{p} and \eqn{q} (see Details for more information); currently,
#'only the default \code{orders = c(1, 1)} is supported; other specifications
#'of a two-element numeric vector will lead to \code{orders = c(1, 1)} being
#'run and a warning message being returned.
#'@param cond_dist the conditional distribution to consider as a
#'character object; the default is a conditional normal distribution
#'\code{"norm"}; available are also, however, a \eqn{t}-distribution
#'(\code{"std"}), a generalized error distribution (\code{"ged"}), an
#'average Laplace distribution (\code{"ald"}),
#'and their four skewed variants (\code{"snorm"}, \code{"sstd"},
#'\code{"sged"}, \code{"sald"}).
#'@param start_pars the starting parameters for the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used.
#'@param meanspec an object of class "mean_spec"; indicates the
#'specifications for the model in the conditional mean.
#'@param Drange a two-element numeric vector that indicates the boundaries
#'of the interval over which to search for the fractional differencing
#'parameter \eqn{D} in a long-memory ARMA-type model in the conditional
#'mean model part; by default,
#'\eqn{D} being searched for on the
#'interval from 0 to \eqn{0.5 - 1\times 10^{-6}}; note that specific
#'settings in the arguments
#'\code{LB} and \code{UB} overwrite this argument.
#'@param nonparspec an object of class \code{"locpol_spec"} returned
#'by \code{\link{locpol_spec()}}; defines the settings of the nonparametric
#'smoothing technique for \code{use_nonpar = TRUE}.
#'@param use_nonpar a logical indicating whether or not to implement a
#'semiparametric extension of the volatility model defined through \code{spec};
#'see "Details" for more information.
#'@param n_test a single numerical value indicating, how many observations
#'at the end of \code{rt} not to include in the fitting process and to
#'reserve for backtesting.
#'@param LB the lower boundaries of the parameters in the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used.
#'@param UB the upper boundaries of the parameters in the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used.
#'@param control a list that is passed to \code{control} of the
#'function \code{solnp} of the package \code{Rsolnp}.
#'@param control_nonpar a list containing changes to the arguments
#'for the hyperparameter estimation algorithm in the nonparametric
#'scale function estimation for
#'\code{use_nonpar = TRUE}; see "Details" for more information.
#'@param parallel only relevant for a (skewed) average Laplace (AL)
#'distribution, i.e.
#'if \code{cond_dist} in \code{spec} is set to \code{cond_dist = "ald"} or
#'\code{cond_dist = "sald"}; \code{parallel} is a logical value indicating whether
#'or not the slices for the positive integer-valued parameter of the SM
#'distribution should be fitted in parallel for a speed boost.
#'@param ncores only relevant for a (skewed) average Laplace (AL)
#'distribution, i.e.
#'if \code{cond_dist} in \code{spec} is set to \code{cond_dist = "ald"} or
#'\code{cond_dist = "sald"}, and if simultaneously \code{parallel = TRUE};
#'\code{ncores} is a single numeric value indicating the number of cores to
#'use for parallel computations.
#'@param Prange a two-element vector that indicates the search boundaries for
#'the parameter \eqn{P} in a (skewed) average Laplace distribution.
#'
#'@details
#'Let \eqn{\left\{r_t\right\}}, with \eqn{t \in \mathbb{Z}} as the
#'time index, be a theoretical time series that follows
#'\deqn{r_t=\mu+\varepsilon_t \text{ with } \varepsilon_t=\sigma_t \eta_t \text{ and } \eta_t \sim \text{IID}(0,1), \text{ where}}
#'\deqn{\sigma_t^{2}=\omega+\sum_{i=1}^{p}\phi_i \varepsilon_{t-i}^2 + \sum_{j=1}^{q}\beta_j \sigma_{t-j}^2.}
#'Here, \eqn{\eta_t\sim\text{IID}(0,1)} means that the innovations
#'\eqn{\eta_t} are independent but identically distributed (iid) with mean zero
#'and variance one, whereas \eqn{\sigma_t > 0} are the conditional standard
#'deviations in \eqn{r_t}.
#'\eqn{\phi_i}, \eqn{i=1,2,\dots, p}, and \eqn{\beta_j}, \eqn{j = 1, \dots, q},
#'are non-negative coefficients. \eqn{p}
#'and \eqn{q} are the model orders definable through the argument \code{orders},
#'where \eqn{p} is the first element and \eqn{q} is the second element in the
#'argument. In addition, we have \eqn{\mu = E\left(r_t\right)} as a
#'real-valued parameter. \eqn{\omega > 0} is the intercept.
#'
#'For a conditional average Laplace distribution, an optimal model for each
#'distribution parameter \eqn{P} from 1 to 8 is estimated (assuming that
#'\eqn{P} is then fixed to the corresponding value). Afterwards, \eqn{P} is then
#'estimated by selecting the estimated model among the eight fitted models that
#'has the largest log-likelihood. The eight models are, by default, fitted
#'simultaneously using parallel programming techniques (see also the arguments
#'\code{parallel} and \code{ncores}, which are only relevant for a conditional
#'average Laplace distribution). After the optimal model (including
#'the estimate of \eqn{P} called \eqn{\hat{P}}) has been determined, \eqn{P=\hat{P}}
#'is seen as fixed to obtain the standard errors via the Hessian matrix for the
#'estimates of the continuous parameters. A standard error for \eqn{\hat{P}} is therefore
#'not obtained and the ones obtained for the remaining estimates do not account
#'for \eqn{\hat{P}}.
#'
#'An ARMA-GARCH or a FARIMA-GARCH can be fitted by adjusting the
#'argument \code{meanspec} correspondingly.
#'
#'As an alternative, a semiparametric extension of the pure models
#'in the conditional variance can be implemented. If \code{use_nonpar = TRUE},
#'\code{meanspec} is omitted and before fitting a zero-mean model in the
#'conditional volatility following the remaining function arguments, a smooth scale function,
#'i.e. a function representing the unconditional standard deviation over time,
#'is being estimated following the specifications in \code{nonparspec} and
#'\code{control_nonpar}. This preliminary step stabilizes the input
#'series \code{rt}, as long-term changes in the unconditional variance
#'are being estimated and removed before the parametric step using
#'\code{\link[esemifar]{tsmoothlm}}. \code{control_nonpar} can be adjusted following
#'to make changes to the arguments of \code{\link[esemifar]{tsmoothlm}}
#'for long-memory specifications. These arguments specify settings
#'for the automated bandwidth selection algorithms implemented by this
#'function. By default, we use the settings
#'\code{pmin = 0}, \code{pmax = 1}, \code{qmin = 0},
#'\code{qmax = 1}, \code{InfR = "Nai"},
#'\code{bStart = 0.15}, \code{cb = 0.05}, and
#'\code{method = "lpr"} for \code{\link[esemifar]{tsmoothlm}}.
#'\code{\link{locpol_spec}} passed to \code{nonparspec} handles
#'more direct settings of the local polynomial smoother itself. See
#'the documentation for these functions to get a detailed overview
#'of these settings. Assume \eqn{\{r_t\}} to be the observed series, where
#'\eqn{t = 1, 2, \dots, n},
#'then \eqn{r_t^{*} = r_t - \bar{r}}, with \eqn{\bar{r}} being the arithmetic
#'mean over the observed \eqn{r_t}, is computed and subsequently
#'\eqn{y_t = \ln\left[\left(r_t^{*}\right)^2\right]}. The subtraction of
#'\eqn{\bar{r}} is necessary so that \eqn{r_t^{*}} are all different from zero
#'almost surely. Once \eqn{y_t} are available, its trend \eqn{m(x_t)},
#'with \eqn{x_t} as the rescaled time on the interval \eqn{[0, 1]}, is
#'being estimated using
#'\code{\link[esemifar]{tsmoothlm}} and denoted here by
#'\eqn{\hat{m}(x_t)}. Then from \eqn{\hat{\xi}_t = y_t - \hat{m}(x_t)}
#'obtain \eqn{\hat{C} = -\ln\left\{\sum_{t=1}^{n}\exp\left(\hat{\xi}_t\right)\right\}},
#'and obtain the estimated scale function as
#'\eqn{\hat{s}(x_t)=\exp\left[\left(\hat{\mu}(x_t) - \hat{C}\right) / 2\right]}.
#'The stabilized / standardized version of the series \eqn{\left\{r_t\right\}}
#'is then \eqn{\tilde{r}_t = r_t^{*} / \hat{s}(x_t)}, to which
#'a purely parametric volatility model following the remaining function arguments
#'is then
#'fitted. The estimated volatility at a given time point is then
#'the product of the estimate of the corresponding scale function value
#'and of the estimated conditional standard deviation (following the parametric
#'model part) for that same time point. See for example Feng et al. (2022)
#'or Letmathe et al. (2023) for more information on the semiparametric extension
#'of volatility models.
#'
#'@export
#'
#'@return
#'An object of S4-class \code{"fEGarch_fit_garch"}
#'is returned. It contains the following elements.
#'\describe{
#'\item{\code{pars}:}{a named numeric vector with the parameter estimates.}
#'\item{\code{se}:}{a named numeric vector with the obtained standard errors in accordance with the parameter estimates.}
#'\item{\code{vcov_mat}:}{the variance-covariance matrix of the parameter estimates with named columns and rows.}
#'\item{\code{rt}:}{the input object \code{rt} (or at least the training data, if \code{n_test} is greater than zero);
#'if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is kept.}
#'\item{\code{cmeans}:}{the estimated conditional means; if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{cmeans}.}
#'\item{\code{sigt}:}{the estimated conditional standard deviations (or for \code{use_nonpar = TRUE} the estimated total
#'volatilities, i.e. scale function value times conditional standard deviation); if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{sigt}.}
#'\item{\code{etat}:}{the obtained residuals; if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{etat}.}
#'\item{\code{orders}:}{a two-element numeric vector stating the considered model orders.}
#'\item{\code{cond_dist}:}{a character value stating the conditional distribution considered in the model fitting.}
#'\item{\code{long_memo}:}{a logical value stating whether or not long memory was considered in the model fitting.}
#'\item{\code{llhood}:}{the log-likelihood value obtained at the optimal parameter combination.}
#'\item{\code{inf_criteria}:}{a named two-element numeric vector with the corresponding AIC (first element) and BIC (second element) of the fitted parametric model part; for purely parametric models, these criteria are valid for the entire model; for semiparametric models, they are only valid for the parametric step and are not valid for the entire model.}
#'\item{\code{meanspec}:}{the settings for the model in the conditional mean; is an object
#'of class \code{"mean_spec"} that is identical to the object passed to the input argument
#'\code{meanspec}.}
#'\item{\code{test_obs}:}{the observations at the end up the input \code{rt} reserved for
#'testing following \code{n_test}.}
#'\item{\code{scale_fun}:}{the estimated scale function values, if \code{use_nonpar = TRUE}, otherwise
#'\code{NULL}; formatting of \code{rt} is reused.}
#'\item{\code{nonpar_model}:}{the estimation object returned by \code{\link[esemifar]{tsmoothlm}} for
#'\code{use_nonpar = TRUE}.}
#'\item{\code{trunc}:}{the input argument \code{trunc}.}
#'}
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- garch(rt)
#'model
#'
#'@references
#'\itemize{
#'\item{Bollerslev, T., & Mikkelsen, H. O. (1996). Modeling and pricing long
#'memory in stock market volatility. Journal of Econometrics, 73(1):
#'151-184. DOI: 10.1016/0304-4076(95)01736-4.}
#'}
#'
garch <- function(rt, orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                        meanspec = mean_spec(), Drange = c(0, 0.5 - 1e-6),
                        nonparspec = locpol_spec(), use_nonpar = FALSE,
                        n_test = 0,
                        start_pars = NULL, LB = NULL, UB = NULL, control = list(),
                        control_nonpar = list(),
                        parallel = TRUE, ncores = max(1, future::availableCores() - 1),
                        trunc = "none", Prange = c(1, 5)) {

  if (use_nonpar) {

    nonpar_out <- run_nonpar(nonparspec = nonparspec,
               control_nonpar = control_nonpar, rt = rt, n_test = n_test)
    rt <- nonpar_out$rt
    meanspec <- nonpar_out$meanspec
    nonpar_result <- nonpar_out$nonpar_result
    n_test <- nonpar_out$n_test

  }

  initial_split(rt, n_test)

  # if (!all(orders == c(1, 1))) {
  #   orders <- c(1, 1)
  #   message("Only orders p=1 with q=1 are currently supported for APARCH models.")
  # }

  n <- length(rt_core)
  if (Prange[[1]] < 1) {Prange[[1]] <- 1}
  trunc_i <- trunc
  if (is.character(trunc) && trunc == "none") {
    trunc <- n - 1
  } else if (is.numeric(trunc) && trunc >= n) {
    trunc <- n - 1
    message("trunc reduced to n - 1")
  }

  p <- orders[[1]]
  q <- orders[[2]]

  cond_d <- match.arg(cond_dist)
  lm <- FALSE

  # Settings for ARMA / FARIMA part in the mean; makes objects available in this
  # environment
  identify_arma_settings(meanspec = meanspec, Drange = Drange)

  dfun1 <- lookup_table$fun1[[cond_d]]
  dfun2 <- lookup_table$fun2[[cond_d]]

  skew_par <- lookup_table[["skew_par"]][[cond_d]]
  extra_par <- lookup_table[["extra_par"]][[cond_d]]

  # Adjust (if ARMA relevant)
  model_type <- "garch"
  sigt_fun_adj <- adj_sigt_fun(
    p_ar = p_ar, q_ma = q_ma, pg0 = pg0, qg0 = qg0,
    lm_arma = lm_arma, incl_mean = incl_mean,
    n_arma_pars = n_arma_pars, grab_fun = arma_grab_fun,
    arma_fun = cmean_fun,
    model_type = model_type
  )

  sig2_init <- 1   # Because of rescaling, var(rt_core) = 1

  goal_fun_s <- function(theta, rt, dfun1, dfun2) {

    goal_fun_creator_garch(
     theta = theta,
     x = rt,
     pdf_fun2 = dfun2,
     p = p,
     q = q,
     p_ar = p_ar,
     q_ma = q_ma,
     incl_mean = incl_mean,
     extra_par = extra_par,
     skew_par = skew_par,
     sig_fun = sigt_fun_adj,
     sig2_init = sig2_init,
     trunc = trunc
    )

  }

  # Get cond. distribution restrictions and starting values
  # as objects sval, lval and uval in the environment
  cond_d_vals <- cond_d_par_restr(cond_d)

  # Make further starting parameter values and restrictions available
  # using internally saved table function
  pars_and_restr <- start_pars_and_restr[[model_type]](drange, p, q)
  list2env(pars_and_restr, envir = environment())
  rm(pars_and_restr)

  if (is.null(start_pars) || is.null(LB) || is.null(UB)) {
    mean_rt <- list(NULL, mean(rt_core))[[incl_mean + 1]]
    min_rt <- list(NULL, min(rt_core))[[incl_mean + 1]]
    max_rt <- list(NULL, max(rt_core))[[incl_mean + 1]]
    skew_start <- list(NULL, 1)[[skew_par + 1]]
    skew_low <- list(NULL, 1e-15)[[skew_par + 1]]
    skew_up <- list(NULL, Inf)[[skew_par + 1]]

    set_start_LB_UB(
      mean_val = mean_rt,
      ar_start = ar_start,
      ma_start = ma_start,
      D_start = D_start,
      omega_val = omega_start,
      phi_start = phi_start,
      psi_start = psi_start,
      add_pars_start = add_pars_start,
      d_start = d_spar,
      sval = sval,
      skew_start = skew_start,
      mean_low = min_rt,
      ar_low = ar_lb,
      ma_low = ma_lb,
      D_low = D_low,
      omega_low = omega_low,
      phi_low = phi_low,
      psi_low = psi_low,
      add_low = add_low,
      d_low = d_low,
      lval = lval, skew_low = skew_low,
      mean_up = max_rt,
      ar_up = ar_ub, ma_up = ma_ub,
      D_up = D_up,
      omega_up = omega_up,
      phi_up = phi_up, psi_up = psi_up,
      add_up = add_up,
      d_up = d_up,
      uval = uval, skew_up = skew_up,
      start_pars = start_pars, LB = LB, UB = UB
    )


  }

  low_phi <- 2 + n_arma_pars
  up_phi <- low_phi - 1 + p
  low_beta <- up_phi + 1
  up_beta <- low_beta - 1 + q

  # For stationarity, restrict phi / AR values
  constr_fun <- function(theta, rt) {
    phi <- theta[low_phi:up_phi]
    beta <- theta[low_beta:up_beta]
    c1 <- sum(phi) + sum(beta)
    if (pg0) {
      ar <- c(1, -theta[low_ar:up_ar])
      c2 <- abs(polyroot(ar)) - 1
    } else {
      c2 <- numeric(0)
    }
    c(c1, c2)
  }

  if (is.null(control$trace)) {control$trace <- FALSE}
  #if (is.null(control$delta)) {control$delta <- 1e-10}

  p_all <- p_ar + 1
  ineqLB <- c(0, rep(1e-6, p_ar))
  ineqUB <- c(1 - 1e-6, rep(Inf, p_ar))

  if (cond_d %in% c("ald", "sald")) {

    # Makes "P_sel" and "result" available
    ald_fit(
      Prange = Prange,
      parallel = parallel,
      ncores = ncores,
      dfun1 = dfun1,
      dfun2 = dfun2,
      goal_fun_s = goal_fun_s,
      start_pars = start_pars,
      LB = LB, UB = UB,
      ineqfun = constr_fun,
      ineqLB = ineqLB,
      ineqUB = ineqUB,
      rt_core = rt_core,
      control = control
    )

  } else {

    P_sel <- NULL

    goal_fun <- function(theta, rt) {
      goal_fun_s(theta = theta, rt = rt, dfun1 = dfun1, dfun2 = dfun2)
    }

    result <- tryCatch(
      expr = {suppressWarnings(Rsolnp::solnp(
        pars = start_pars,
        fun = goal_fun,
        LB = LB,
        UB = UB,
        ineqfun = constr_fun,
        ineqLB = ineqLB,    # abs. values of roots must be outside of unit circle
        ineqUB = ineqUB,
        control = control,
        rt = rt_core
      ))},
      error = function(e1) {
        stop("Error during optimization. You may want to try different starting parameter and / or solver settings.", call. = FALSE)
      }
    )

  }

  if (result$convergence %in% c(1, 2)) {
    warning("Convergence failed. You may want to try different starting parameter and / or solver settings.", call. = FALSE)
  }

  pars <- result$pars

  # Compute log-likelihood which is shifted from the obtained negative
  # log-likelihood at the optimum; "-tail(result$values, 1)" is the
  # log-likelihood at the optimum for the scaled data, which is too large by
  # "n * log(scale_const)" in comparison to original data
  llhood <- -tail(result$values, 1) - n * log(scale_const)

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

  par_names <- c(
    list(NULL, "mu")[[incl_mean + 1]],
    names_ar,
    names_ma,
    D_par_name,
    "omega",
    paste0("phi", 1:p),
    beta_names,
    extra_par_name,
    list(NULL, "skew")[[skew_par + 1]]
  )

  if (cond_d %in% c("ald", "sald")) {

    # Fix dfun2
    dfun2_P <- function(x, mu, sigt, shape, skew) {
      dfun2(x = x, mu = mu, sigt = sigt, shape = P_sel, skew = skew)
    }

    # Fix the final goal function to optimize over
    # using dfun1_P and dfun2_P
    goal_fun <- function(theta, rt) {
      goal_fun_s(theta = theta, rt = rt, dfun2 = dfun2_P, dfun1 = NULL)
    }

  }

  # Compute the conditional standard deviations at the optimum
  final_series <- sigt_fun_adj(
    theta = unname(pars),
    x = rt_core,
    dfun = dfun1,
    p = p,
    q = q,
    p_ar = p_ar,
    q_ma = q_ma,
    incl_mean = incl_mean,
    extra_par = extra_par,
    skew_par = skew_par,
    trunc = trunc,
    sig2_init = sig2_init
  )

  delta <- 2                 # treat as APARCH with delta fixed at 2 for rescaling
  sc_delta <- scale_const^delta
  # Also updates "pars", if ald or sald
  compute_vcov(
    goal_fun = goal_fun,
    pars = pars,
    rt_core = rt_core,
    scale_const = scale_const,
    incl_mean = incl_mean,
    P_sel = P_sel,
    cond_d = cond_d,
    par_names = par_names,
    sc_delta = sc_delta,
    n_arma_pars = n_arma_pars,
    model_type = "aparch"     # treat as APARCH with delta fixed at 2 for rescaling
  )

  rescale_estim(
    pars = pars,
    rt_core_o = rt_core_o,
    sigt = final_series$sigt,
    cmeans = final_series$cmeans,
    scale_const = scale_const,
    n_arma_pars = n_arma_pars,
    incl_mean = incl_mean,
    model_type = "aparch",  # treat as APARCH with delta fixed at 2 for rescaling
    sc_delta = sc_delta
  )

  names(pars) <- par_names

  # Makes "aic" and "bic" available
  crit_calc(pars = pars, rt = rt_core_o, llhood = llhood)

  # Apply time series formatting (if relevant) to all
  # estimated series
  series_out <- format_applier_ts(
    rt = train_obs,
    list_of_ts = list(
      "cmeans" = cmeans,
      "sigt" = sigt,
      "etat" = etat
    )
  )
  cmeans <- series_out$cmeans
  sigt <- series_out$sigt
  etat <- series_out$etat

  list_out <- fEGarch_fit_garch(
    pars = pars,
    se = serrors,
    vcov_mat = vcov_mat,
    rt = train_obs,
    sigt = sigt,
    cmeans = cmeans,
    scale_fun = NULL,
    etat = etat,
    llhood = llhood,
    inf_criteria = c("aic" = aic, "bic" = bic),
    cond_dist = cond_d,
    orders = orders,
    long_memo = FALSE,
    meanspec = meanspec,
    test_obs = test_obs,
    nonpar_model = NULL,
    trunc = trunc_i
  )

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

#'@export
#'@rdname show-methods
#'@aliases show,fEgarch_fit_garch-method
setMethod(
  "show",
  "fEGarch_fit_garch",
  function(object) {
  x <- object
  par_names <- names(x@pars)
  se <- unname(x@se)
  pars <- unname(x@pars)

  tvals <- pars / se
  atvals <- abs(tvals)
  pvals <- 2 * (1 - pnorm(atvals))

  df <- data.frame(
    par = sprintf("%.4f", pars),
    se = sprintf("%.4f", se),
    tval = sprintf("%.4f", tvals),
    pval = sprintf("%.4f", pvals)
  )
  row.names(df) <- par_names

  check_nonpar <- !is.null(x@nonpar_model)
  b <- if (is.null(x@nonpar_model$b0)) {
    x@nonpar_model$b
  } else {
    x@nonpar_model$b0
  }
  bwidth <- if (check_nonpar) {
    paste0(
      " (poly_order = ", x@nonpar_model$p,"; bandwidth = ", sprintf("%.4f", b), ")"
    )
  } else {
    ""
  }

  part1 <- paste0(
    "************************************\n",
    "*        Fitted GARCH Model        *\n",
    "************************************\n",
    " \n",
    "Type: sGARCH\n",
    "Orders: (", x@orders[[1]], ", ", x@orders[[2]], ")\n",
    "Long memory: ", x@long_memo, "\n",
    "Cond. distribution: ", x@cond_dist, "\n",
    " \n",
    "ARMA orders (cond. mean): (", x@meanspec@orders[[1]], ", ", x@meanspec@orders[[2]], ")\n",
    "Long memory (cond. mean): ", x@meanspec@long_memo, "\n",
    " \n",
    "Scale estimation: ", paste0(check_nonpar, bwidth), "\n",
    " \n",
    "Fitted parameters:\n"
  )

  part2 <- paste0(
    " \n",
    "Information criteria (parametric part):\n",
    "AIC: ", sprintf("%.4f", unname(x@inf_criteria[[1]])),
    ", BIC: ", sprintf("%.4f", unname(x@inf_criteria[[2]]))
  )

  cat(part1, fill = TRUE)
  print(df)
  cat(part2, fill = TRUE)
  }
)
