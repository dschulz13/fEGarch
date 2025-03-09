aparch_grab_pars <- function(theta, p, q, incl_mean, extra_par, skew_par) {

  imean <- !incl_mean
  mu <- if (incl_mean) theta[[1]] else 0
  step0 <- 2 - imean
  omega <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q
  step3 <- step2 + p
  phi <- theta[(step0 + 1):step1]
  beta <- theta[(step1 + 1):step2]
  gamma <- theta[(step2 + 1):step3]
  delta <- theta[[step3 + 1]]
  shape <- if (extra_par) theta[[step3 + 2]] else 0
  skew <- if (skew_par) theta[[step3 + 2 + extra_par]] else 0

  list(
    mu = mu,
    omega = omega,
    phi = phi,
    beta = beta,
    gamma = gamma,
    delta = delta,
    shape = shape,
    skew = skew
  )

}

sigt_aparch_R <- function(theta, x, p, q, p_ar, q_ma,
              incl_mean, extra_par, skew_par, dfun, sigd_init, etransf_init,
              trunc) {

  all_pars <- aparch_grab_pars(
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
  gamma <- all_pars$gamma
  delta <- all_pars$delta
  shape <- all_pars$shape
  skew <- all_pars$skew

  # Initialization
  Esigd <- sd(x)^delta
  E_e <- 0

  sigt <- c(sigt_aparch_Cpp(
    x = x,
    mu = mu,
    phi = phi,
    beta = beta,
    gamma = gamma,
    delta = delta,
    omega = omega,
    sigd_init,
    etransf_init
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape, cmeans = rep(mu, length(sigt)))

}

goal_fun_creator_aparch <- function(
    theta,
    x,
    pdf_fun1,
    pdf_fun2,
    p,
    q,
    p_ar,
    q_ma,
    incl_mean,
    extra_par,
    skew_par,
    sig_fun,
    sigd_init,
    etransf_init,
    trunc
) {

  res_list <- sig_fun(theta, x, p, q, p_ar, q_ma,
              incl_mean, extra_par, skew_par, pdf_fun1,
              sigd_init, etransf_init, trunc)

  sigt <- res_list$sigt

  cmeans <- res_list$cmeans
  skew <- res_list$skew
  shape <- res_list$shape

  out <- -sum(log(pdf_fun2(x, cmeans, sigt, shape, skew)))

  out

}


#'APARCH Model Fitting
#'
#'Fit an APARCH model under the six
#'most common conditional distributions to observed data
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
#'
#'@details
#'Let \eqn{\left\{r_t\right\}}, with \eqn{t \in \mathbb{Z}} as the
#'time index, be a theoretical time series that follows
#'\deqn{r_t=\mu+\varepsilon_t \text{ with } \varepsilon_t=\sigma_t \eta_t \text{ and } \eta_t \sim \text{IID}(0,1), \text{ where}}
#'\deqn{\sigma_t^{\delta}=\omega+\left[1-\beta^{-1}(B)\phi(B)(1-B)^{d}\right]\left(\left|\varepsilon_t\right|-\gamma\varepsilon_t\right)^{\delta}.}
#'Here, \eqn{\eta_t\sim\text{IID}(0,1)} means that the innovations
#'\eqn{\eta_t} are independent but identically distributed (iid) with mean zero
#'and variance one, whereas \eqn{\sigma_t > 0} are the conditional standard
#'deviations in \eqn{r_t}.
#'Moreover, \eqn{B} is the backshift operator and
#'\eqn{\beta(B) = 1 - \sum_{j=1}^{q}\beta_j B^{j}}, where
#'\eqn{\beta_j}, \eqn{j=1,2,\dots, q}, are real-valued coefficients. Furthermore,
#'\eqn{\phi(B) = 1 - \sum_{i=1}^{p}\phi_i B^{i}}, where
#'\eqn{\phi_i}, \eqn{i=1,2,\dots, p}, are real-valued coefficients. \eqn{p}
#'and \eqn{q} are the model orders definable through the argument \code{orders},
#'where \eqn{p} is the first element and \eqn{q} is the second element in the
#'argument. In addition, we have \eqn{\mu = E\left(r_t\right)} as a
#'real-valued parameter, \eqn{\gamma \in (-1,1)} and \eqn{d \in [0,0.5)} as the
#'long-memory parameter. \eqn{\omega > 0} is the intercept. It is assumed that
#'all \eqn{\beta_j} and \eqn{\phi_i} are non-negative.
#'
#'The pre-sample values of \eqn{\left(\left|\varepsilon_t\right|-\gamma\varepsilon_t\right)^{\delta}}
#'are replaced by the in-sample arithmetic mean of
#'\eqn{\left(\left|r_t^{*} - \mu\right|-\gamma\left(r_t^{*}-\mu\right)\right)^{\delta}}, where
#'\eqn{r_t^{*}} are considered as the observations and \eqn{\mu}, \eqn{\gamma} and
#'\eqn{\delta} are not fixed values but also to be optimized over simultaneously
#'in the optimization
#'problem. The pre-sample is set to length 1000, and the infinite coefficient series
#'\eqn{\left[1-\beta^{-1}(B)\phi(B)(1-B)^{d}\right]} is, at each time point, always considered
#'as far back as to the first pre-sample value. Therefore, at observation time point
#'\eqn{t=1,\dots,n}, a truncation of 1000 + t - 1 coefficients is considered.
#'
#'Currently, only a model of orders \eqn{p=1} with \eqn{q=1} can be fitted;
#'to ensure the non-negativity of all of the infinite-order coefficient series
#'\eqn{\left[1-\beta^{-1}(B)\phi(B)(1-B)^{d}\right]}, which in combination
#'with \eqn{\omega>0} ensures that all the conditional volatilities
#'are greater than zero, we employ the sufficient conditions mentioned
#'in Footnote 6 in Tse (1998) in the optimization problem, which are
#'taken from Bollerslev and Mikkelsen (1996). Conrad and Haag (2006) provide
#'necessary and sufficient conditions for cases with simple model orders and
#'sufficient conditions for general model orders, however, due to
#'time efficiency, this function uses the set of sufficient conditions
#'by Bollerslev and Mikkelsen (1996) only.
#'
#'See also the reference section for sources on the APARCH (Ding et al., 1993)
#'and the FIAPARCH (Tse, 1998) models.
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
#'An ARMA-FIAPARCH or a FARIMA-FIAPARCH can be fitted by adjusting the
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
#'for the automated bandwidth selection algorithms implemented by thís
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
#'An object of S4-class \code{"fEGarch_fit_aparch"}
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
#'}
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- aparch(rt)
#'model
#'
#'@references
#'\itemize{
#'\item{Bollerslev, T., & Mikkelsen, H. O. (1996). Modeling and pricing long
#'memory in stock market volatility. Journal of Econometrics, 73(1):
#'151-184. DOI: 10.1016/0304-4076(95)01736-4.}
#'\item{Conrad, C., & Haag, B. R. (2006). Inequality constraints in the fractionally
#'integrated GARCH model. Journal of Financial Econometrics, 4(3):
#'413-449. DOI: 10.1093/jjfinec/nbj015.}
#'\item{Ding, Z., Granger, C. W. J., & Engle, R. F. (1993). A long memory property of stock market returns
#'and a new model. Journal of Empirical Finance, 1(1):
#'83-106. DOI: 10.1016/0927-5398(93)90006-D.}
#'\item{Feng, Y., Gries, T., Letmathe, S., & Schulz, D. (2022). The smoots Package in R for Semiparametric Modeling of
#'Trend Stationary Time Series. The R Journal,
#'14(1), 182-195. URL: https://journal.r-project.org/articles/RJ-2022-017/.}
#'\item{Letmathe, S., Beran, J., & Feng, Y. (2023). An extended exponential SEMIFAR model with application
#'in R. Comunications in Statistics - Theory and Methods,
#'53(22), 7914–7926. DOI: 10.1080/03610926.2023.2276049.}
#'\item{Tse, Y. K. (1998). The conditional heteroskedasticity of the
#'yen-dollar exchange rate. Journal of Applied Econometrics, 13(1):
#'49-55. DOI: 10.1002/(SICI)1099-1255(199801/02)13:1<49::AID-JAE459>3.0.CO;2-O.}
#'}
#'
aparch <- function(rt, orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                        meanspec = mean_spec(), Drange = c(0, 0.5 - 1e-6),
                        nonparspec = locpol_spec(), use_nonpar = FALSE,
                        n_test = 0,
                        start_pars = NULL, LB = NULL, UB = NULL, control = list(),
                        control_nonpar = list(),
                        parallel = TRUE, ncores = max(1, future::availableCores() - 1),
                        trunc = 1000) {

  if (use_nonpar) {

    nonpar_result <- nonpar_est(
      rt = rt,
      lm = FALSE,
      nonparspec = nonparspec,
      n_test = n_test,
      control_nonpar = control_nonpar
    )

    meanspec <- mean_spec(
      orders = c(0, 0),
      long_memo = FALSE,
      include_mean = FALSE
    )    # Force to not estimate the mean in parametric part
    n_test <- 0
    rt <- nonpar_result$et

  }

  split_rt <- split_ts(rt, n_test)
  rt_core_o <- zoo::coredata(split_rt$train)

  scale_const <- sd(rt_core_o)         # for rescaling later
  rt_core <- rt_core_o / scale_const   # now has sample variance 1

  # if (!all(orders == c(1, 1))) {
  #   orders <- c(1, 1)
  #   message("Only orders p=1 with q=1 are currently supported for APARCH models.")
  # }

  p <- orders[[1]]
  q <- orders[[2]]

  cond_d <- match.arg(cond_dist)

  # ARMA orders
  order_arma <- orders(meanspec)
  p_ar <- order_arma[[1]]
  q_ma <- order_arma[[2]]

  pdf_fun1 <- lookup_table$fun1[[cond_d]]
  pdf_fun2 <- lookup_table$fun2[[cond_d]]

  incl_mean <- meanspec@include_mean

  skew_par <- lookup_table[["skew_par"]][[cond_d]]
  extra_par <- lookup_table[["extra_par"]][[cond_d]]

  # Long-memory setting (for ARMA-part)
  lm_arma <- long_memo(meanspec)
  # Assign ARMA / FARIMA function
  arma_fun <- if (lm_arma) {
    grab_fun <- farima_grab_pars
    D_par_name <- "D"
    D_start <- 0.25 #mean(Drange)
    D_low <- Drange[[1]]
    D_up <- Drange[[2]]
    farima_fit
  } else {
    grab_fun <- arma_grab_pars
    D_par_name <- character(0)
    D_start <- D_low <- D_up <- numeric(0)
    arma_fit
  }


  n <- length(rt_core)
  sigt_fun <- sigt_aparch_R
  # Adjust (if ARMA relevant)
  pg0 <- p_ar > 0
  qg0 <- q_ma > 0
  check <- pg0 || qg0
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
        esemifar::farima_to_ar(ar = ar, ma = ma, d = arma_pars$D, max_i = trunc)[-1]
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
      out <- sigt_fun(theta2, res, p, q, p_ar, q_ma, FALSE, extra_par, skew_par, dfun, sigd_init, etransf_init, trunc)
      out$cmeans <- fitted_vs
      out
    }

  } else {
    sigt_fun             # ...else keep function for EGARCH-type
  }

  delta0 <- 2    # use delta0 = 2
  gamma0 <- 0  # use gamma0 = 0   as suitable initial guesses
  sigd_init <- sd(rt_core)^delta0
  rt_dm <- rt_core - mean(rt_core)
  etransf_init <- mean(abs(rt_dm) - gamma0 * rt_dm)^delta0

  goal_fun_s <- function(theta, rt, pdf_fun2) {

    goal_fun_creator_aparch(
     theta = theta,
     x = rt,
     pdf_fun1 = pdf_fun1,
     pdf_fun2 = pdf_fun2,
     p = p,
     q = q,
     p_ar = p_ar,
     q_ma = q_ma,
     incl_mean = incl_mean,
     extra_par = extra_par,
     skew_par = skew_par,
     sig_fun = sigt_fun_adj,
     sigd_init = sigd_init,
     etransf_init = etransf_init,
     trunc = trunc
    )

  }

  vals <- lookup_table[["vals"]][[cond_d]]

  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  if (is.null(start_pars)) {

    start_pars <- c(
      list(NULL, mean(rt_core))[[incl_mean + 1]],
      rep(0.05 / p_ar, p_ar),
      rep(0.05 / q_ma, q_ma),
      D_start,
      var(rt_core),
      rep(0.05 / p, p),
      rep(0.9 / q, q),
      rep(0.3, p), 2,
      sval,
      list(NULL, 1)[[skew_par + 1]]
    )

  }

  if (is.null(LB)) {

    LB <- c(
      list(NULL, min(rt_core))[[incl_mean + 1]],
      rep(-2, p_ar),
      rep(-2, q_ma),
      D_low,
      1e-10,
      rep(1e-6, p),
      rep(1e-6, q),
      rep(-1 + 1e-4, p), 0.3,
      lval,
      list(NULL, 1e-15)[[skew_par + 1]]
    )

  }

  if (is.null(UB)) {

    UB <- c(
      list(NULL, max(rt_core))[[incl_mean + 1]],
      rep(2, p_ar),
      rep(2, q_ma),
      D_up,
      10,
      rep(1 - 1e-6, p),
      rep(1 - 1e-6, q),
      rep(1 - 1e-4, p), 3,
      uval,
      list(NULL, Inf)[[skew_par + 1]]
    )

  }

  n_arma <- p_ar + q_ma
  n_arma_im <- n_arma + incl_mean + lm_arma
  low_ar <- 1 + incl_mean
  up_ar <- low_ar + p_ar - 1
  low_phi <- 2 + n_arma_im
  up_phi <- low_phi
  low_beta <- low_phi + 1
  up_beta <- low_beta
  #d_idx <- low_beta + 3

  # ineq_fun <- function(theta, rt) {
  #
  #   phi <- theta[low_phi:up_phi]
  #   beta <- theta[low_beta:up_beta]
  #   d <- theta[[d_idx]]
  #
  #   r1 <- phi - beta + d
  #   r2 <- (2 - d) / 3 - phi
  #   r3 <- beta * r1 - d * (phi - (1 - d) / 2)
  #
  #   if (pg0) {
  #     ar <- c(1, -theta[low_ar:up_ar])
  #     c2 <- abs(polyroot(ar)) - 1
  #   } else {
  #     c2 <- numeric(0)
  #   }
  #
  #   c(r1, r2, r3, c2)
  #
  # }

  if (is.null(control$trace)) {control$trace <- FALSE}
  #if (is.null(control$delta)) {control$delta <- 1e-10}

  if (cond_d %in% c("ald", "sald")) {

    P_range <- c(1, 8)
    oldplan <- future::plan()
    if (parallel) {
      future::plan(future::multisession, workers = ncores)
      on.exit(expr = {future::plan(oldplan)}, add = TRUE, after = TRUE)
    }

    P_seq <- P_range[[1]]:P_range[[2]]

    fitted_models <- furrr::future_map(
      .x = P_seq,
      .f = function(.x, start_pars, goal_fun_s, LB, UB, control, rt_core, dfun2, p_ar) {

        dfun2_P <- function(x, mu, sigt, shape, skew) {
          dfun2(x = x, mu = mu, sigt = sigt, shape = .x, skew = skew)
        }

        goal_fun <- function(theta, rt) {
          goal_fun_s(theta = theta, rt = rt, pdf_fun2 = dfun2_P)
        }

        est <- suppressWarnings(Rsolnp::solnp(
          pars = start_pars,
          fun = goal_fun,
          LB = LB,
          UB = UB,
          # ineqfun = ineq_fun,
          # ineqLB = c(0, 0, 0, rep(1e-6, p_ar)),
          # ineqUB = c(Inf, Inf, Inf, rep(Inf, p_ar)),
          control = control,
          rt = rt_core
        ))
        est
      }, start_pars = start_pars, goal_fun_s = goal_fun_s, LB = LB, UB = UB,
      control = control, rt_core = rt_core,#ineq_fun = ineq_fun,
      dfun2 = pdf_fun2, p_ar = p_ar,
      .progress = FALSE, .options = furrr::furrr_options(seed = NULL)
    )

    future::plan(oldplan)

    coll_nllh <- vapply(
      X = fitted_models,
      FUN = function(.x) {
        tail(.x$values, 1)
      },
      FUN.VALUE = numeric(1)
    )

    P_sel <- P_seq[[which.min(coll_nllh)]]

    result <- fitted_models[[which.min(coll_nllh)]]

  } else {

    goal_fun <- function(theta, rt) {
      goal_fun_s(theta = theta, rt = rt, pdf_fun2 = pdf_fun2)
    }

    # result <- tryCatch(
    #   expr = {
    #     suppressWarnings(Rsolnp::solnp(
    #       pars = start_pars,
    #       fun = goal_fun,
    #       LB = LB,
    #       UB = UB,
    #       control = control,
    #       # ineqfun = ineq_fun,
    #       # ineqLB = c(0, 0, 0, rep(1e-6, p_ar)),
    #       # ineqUB = c(Inf, Inf, Inf, rep(Inf, p_ar)),
    #       rt = rt_core
    #     ))
    #   },
    #   error = function(e1) {
    #     stop("Convergence failed.", call. = FALSE)
    #   }
    # )

    result <- suppressWarnings(Rsolnp::solnp(
          pars = start_pars,
          fun = goal_fun,
          LB = LB,
          UB = UB,
          control = control,
          # ineqfun = ineq_fun,
          # ineqLB = c(0, 0, 0, rep(1e-6, p_ar)),
          # ineqUB = c(Inf, Inf, Inf, rep(Inf, p_ar)),
          rt = rt_core
        ))

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

  names_ar <- list(character(0), paste0("ar", 1:p_ar))[[pg0 + 1]]
  names_ma <- list(character(0), paste0("ma", 1:q_ma))[[qg0 + 1]]

  par_names <- c(
    list(NULL, "mu")[[incl_mean + 1]],
    names_ar,
    names_ma,
    D_par_name,
    "omega",
    paste0("phi", 1:p),
    beta_names,
    paste0("gamma", 1:p),
    "delta",
    extra_par_name,
    list(NULL, "skew")[[skew_par + 1]]
  )

  nhessian <- tryCatch(
    expr = {hessCalc(func = goal_fun, pars = unname(pars), rt = rt_core,
                     control.args = list(
				eps = 1e-4,
				d = 1e-5,
				zero.tol = sqrt(.Machine$double.eps / 7e-7),
				r = 4,
				v = 2,
				show.details = FALSE))},
    error = function(e1) {
      warning("Unable to obtain standard errors.", call. = FALSE)
      matrix(NA_real_, ncol = length(pars), nrow = length(pars))
    }
  )

  I_mat <- diag(length(pars))
  delta <- pars[[up_beta + 2]]
  if (incl_mean) {
    I_mat[1, 1] <- scale_const
  }
  n_arma_im1 <- n_arma_im + 1
  sc_delta <- scale_const^delta
  I_mat[n_arma_im1, n_arma_im1] <- sc_delta

  vcov_mat <- tryCatch(
    expr = {I_mat %*% solve(nhessian) %*% I_mat},
    error = function(e1) {
      warning("Unable to obtain standard errors.", call. = FALSE)
      matrix(NA_real_, ncol = length(pars), nrow = length(pars))
    }
  )

  diag_vcov <- diag(vcov_mat)
  if (all(!is.na(diag_vcov)) && any(diag_vcov < 0)) {
    vcov_mat <- matrix(NA, nrow = length(pars), ncol = length(pars))
    diag_vcov <- rep(NA, length(pars))
    warning("Unable to compute Hessian matrix. No standard errors available as consequence.")
  }

  res_f <- sigt_fun_adj(unname(pars), rt_core, p, q, p_ar, q_ma,
              incl_mean, extra_par, skew_par, dfun = pdf_fun1,
              sigd_init, etransf_init, trunc)

  # Adjust parameters according to scaling where necessary
  if (incl_mean) {
    pars[[1]] <- pars[[1]] * scale_const    # mean adjustment
  }
  # omega
  pars[[n_arma_im1]] <- pars[[n_arma_im1]] * sc_delta

  sigt <- res_f$sigt * scale_const
  cmeans <- res_f$cmeans * scale_const

  #mu <- c(0, pars[[1]])[[incl_mean + 1]]
  etat <- (rt_core_o - cmeans) / sigt

  if (cond_d %in% c("ald", "sald")) {
    if (cond_d == "ald") {
      pars <- c(pars, P_sel)
      vcov_mat <- rbind(cbind(vcov_mat, NA), NA)
    } else {
      l_p <- length(pars)
      pars <- c(pars[1:(l_p - 1)], P_sel, pars[[l_p]])
      vcov_c <- cbind(vcov_mat[, 1:(l_p - 1)], NA, vcov_mat[, l_p])
      vcov_mat <- rbind(vcov_c[1:(l_p - 1), ], NA, vcov_c[l_p, ])
    }
  }

  serrors <- sqrt(diag_vcov)

  names(pars) <- par_names
  names(serrors) <- par_names
  colnames(vcov_mat) <- par_names
  rownames(vcov_mat) <- par_names

  k <- length(pars)
  n <- length(rt_core)

  m2llhood <- -2 * llhood
  bic <- (m2llhood + k * log(n)) / n
  aic <- (m2llhood + 2 * k) / n

  # Apply time series formatting (if relevant) to all
  # estimated series
  series_out <- format_applier_ts(
    rt = split_rt$train,
    list_of_ts = list(
      "cmeans" = cmeans,
      "sigt" = sigt,
      "etat" = etat
    )
  )
  cmeans <- series_out$cmeans
  sigt <- series_out$sigt
  etat <- series_out$etat

  list_out <- fEGarch_fit_aparch(
    pars = pars,
    se = serrors,
    vcov_mat = vcov_mat,
    rt = split_rt$train,
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
    test_obs = split_rt$test,
    nonpar_model = NULL
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
#'@aliases show,fEgarch_fit_aparch-method
setMethod(
  "show",
  "fEGarch_fit_aparch",
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
    "*************************************\n",
    "*        Fitted APARCH Model        *\n",
    "*************************************\n",
    " \n",
    "Type: apARCH\n",
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

aparch_pre <- function(rt, orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                        meanspec = mean_spec(), Drange = c(0, 0.5 - 1e-6),
                        nonparspec = locpol_spec(), use_nonpar = FALSE,
                        n_test = 0,
                        start_pars = NULL, LB = NULL, UB = NULL, control = list(),
                        control_nonpar = list(),
                        parallel = TRUE, ncores = max(1, future::availableCores() - 1),
                        trunc = 1000) {

  if (use_nonpar) {

    nonpar_result <- nonpar_est(
      rt = rt,
      lm = FALSE,
      nonparspec = nonparspec,
      n_test = n_test,
      control_nonpar = control_nonpar
    )

    meanspec <- mean_spec(
      orders = c(0, 0),
      long_memo = FALSE,
      include_mean = FALSE
    )    # Force to not estimate the mean in parametric part
    n_test <- 0
    rt <- nonpar_result$et

  }

  split_rt <- split_ts(rt, n_test)
  rt_core_o <- zoo::coredata(split_rt$train)

  scale_const <- sd(rt_core_o)         # for rescaling later
  rt_core <- rt_core_o / scale_const   # now has sample variance 1

  # if (!all(orders == c(1, 1))) {
  #   orders <- c(1, 1)
  #   message("Only orders p=1 with q=1 are currently supported for APARCH models.")
  # }

  p <- orders[[1]]
  q <- orders[[2]]

  cond_d <- match.arg(cond_dist)

  # ARMA orders
  order_arma <- orders(meanspec)
  p_ar <- order_arma[[1]]
  q_ma <- order_arma[[2]]

  pdf_fun1 <- lookup_table$fun1[[cond_d]]
  pdf_fun2 <- lookup_table$fun2[[cond_d]]

  incl_mean <- meanspec@include_mean

  skew_par <- lookup_table[["skew_par"]][[cond_d]]
  extra_par <- lookup_table[["extra_par"]][[cond_d]]

  # Long-memory setting (for ARMA-part)
  lm_arma <- long_memo(meanspec)
  # Assign ARMA / FARIMA function
  arma_fun <- if (lm_arma) {
    grab_fun <- farima_grab_pars
    D_par_name <- "D"
    D_start <- 0.25 #mean(Drange)
    D_low <- Drange[[1]]
    D_up <- Drange[[2]]
    farima_fit
  } else {
    grab_fun <- arma_grab_pars
    D_par_name <- character(0)
    D_start <- D_low <- D_up <- numeric(0)
    arma_fit
  }


  n <- length(rt_core)
  sigt_fun <- sigt_aparch_R
  # Adjust (if ARMA relevant)
  pg0 <- p_ar > 0
  qg0 <- q_ma > 0
  check <- pg0 || qg0
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
        esemifar::farima_to_ar(ar = ar, ma = ma, d = arma_pars$D, max_i = trunc)[-1]
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
      out <- sigt_fun(theta2, res, p, q, p_ar, q_ma, FALSE, extra_par, skew_par, dfun, sigd_init, etransf_init, trunc)
      out$cmeans <- fitted_vs
      out
    }

  } else {
    sigt_fun             # ...else keep function for EGARCH-type
  }

  delta0 <- 2    # use delta0 = 2
  gamma0 <- 0  # use gamma0 = 0   as suitable initial guesses
  sigd_init <- sd(rt_core)^delta0
  rt_dm <- rt_core - mean(rt_core)
  etransf_init <- mean(abs(rt_dm) - gamma0 * rt_dm)^delta0

  goal_fun_s <- function(theta, rt, pdf_fun2) {

    goal_fun_creator_aparch(
     theta = theta,
     x = rt,
     pdf_fun1 = pdf_fun1,
     pdf_fun2 = pdf_fun2,
     p = p,
     q = q,
     p_ar = p_ar,
     q_ma = q_ma,
     incl_mean = incl_mean,
     extra_par = extra_par,
     skew_par = skew_par,
     sig_fun = sigt_fun_adj,
     sigd_init = sigd_init,
     etransf_init = etransf_init,
     trunc = trunc
    )

  }

  vals <- lookup_table[["vals"]][[cond_d]]

  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  if (is.null(start_pars)) {

    start_pars <- c(
      list(NULL, mean(rt_core))[[incl_mean + 1]],
      rep(0.05 / p_ar, p_ar),
      rep(0.05 / q_ma, q_ma),
      D_start,
      var(rt_core),
      rep(0.05 / p, p),
      rep(0.9 / q, q),
      rep(0.3, p), 2,
      sval,
      list(NULL, 1)[[skew_par + 1]]
    )

  }

  if (is.null(LB)) {

    LB <- c(
      list(NULL, min(rt_core))[[incl_mean + 1]],
      rep(-2, p_ar),
      rep(-2, q_ma),
      D_low,
      1e-10,
      rep(1e-6, p),
      rep(1e-6, q),
      rep(-1 + 1e-4, p), 0.3,
      lval,
      list(NULL, 1e-15)[[skew_par + 1]]
    )

  }

  if (is.null(UB)) {

    UB <- c(
      list(NULL, max(rt_core))[[incl_mean + 1]],
      rep(2, p_ar),
      rep(2, q_ma),
      D_up,
      10,
      rep(1 - 1e-6, p),
      rep(1 - 1e-6, q),
      rep(1 - 1e-4, p), 3,
      uval,
      list(NULL, Inf)[[skew_par + 1]]
    )

  }

  n_arma <- p_ar + q_ma
  n_arma_im <- n_arma + incl_mean + lm_arma
  low_ar <- 1 + incl_mean
  up_ar <- low_ar + p_ar - 1
  low_phi <- 2 + n_arma_im
  up_phi <- low_phi
  low_beta <- low_phi + 1
  up_beta <- low_beta
  #d_idx <- low_beta + 3

  # ineq_fun <- function(theta, rt) {
  #
  #   phi <- theta[low_phi:up_phi]
  #   beta <- theta[low_beta:up_beta]
  #   d <- theta[[d_idx]]
  #
  #   r1 <- phi - beta + d
  #   r2 <- (2 - d) / 3 - phi
  #   r3 <- beta * r1 - d * (phi - (1 - d) / 2)
  #
  #   if (pg0) {
  #     ar <- c(1, -theta[low_ar:up_ar])
  #     c2 <- abs(polyroot(ar)) - 1
  #   } else {
  #     c2 <- numeric(0)
  #   }
  #
  #   c(r1, r2, r3, c2)
  #
  # }

  if (is.null(control$trace)) {control$trace <- FALSE}
  #if (is.null(control$delta)) {control$delta <- 1e-10}

  if (cond_d %in% c("ald", "sald")) {

    P_range <- c(1, 8)
    oldplan <- future::plan()
    if (parallel) {
      future::plan(future::multisession, workers = ncores)
      on.exit(expr = {future::plan(oldplan)}, add = TRUE, after = TRUE)
    }

    P_seq <- P_range[[1]]:P_range[[2]]

    fitted_models <- furrr::future_map(
      .x = P_seq,
      .f = function(.x, start_pars, goal_fun_s, LB, UB, control, rt_core, dfun2, p_ar) {

        dfun2_P <- function(x, mu, sigt, shape, skew) {
          dfun2(x = x, mu = mu, sigt = sigt, shape = .x, skew = skew)
        }

        goal_fun <- function(theta, rt) {
          goal_fun_s(theta = theta, rt = rt, pdf_fun2 = dfun2_P)
        }

        est <- suppressWarnings(Rsolnp::solnp(
          pars = start_pars,
          fun = goal_fun,
          LB = LB,
          UB = UB,
          # ineqfun = ineq_fun,
          # ineqLB = c(0, 0, 0, rep(1e-6, p_ar)),
          # ineqUB = c(Inf, Inf, Inf, rep(Inf, p_ar)),
          control = control,
          rt = rt_core
        ))
        est
      }, start_pars = start_pars, goal_fun_s = goal_fun_s, LB = LB, UB = UB,
      control = control, rt_core = rt_core,#ineq_fun = ineq_fun,
      dfun2 = pdf_fun2, p_ar = p_ar,
      .progress = FALSE, .options = furrr::furrr_options(seed = NULL)
    )

    future::plan(oldplan)

    coll_nllh <- vapply(
      X = fitted_models,
      FUN = function(.x) {
        tail(.x$values, 1)
      },
      FUN.VALUE = numeric(1)
    )

    P_sel <- P_seq[[which.min(coll_nllh)]]

    result <- fitted_models[[which.min(coll_nllh)]]

  } else {

    goal_fun <- function(theta, rt) {
      goal_fun_s(theta = theta, rt = rt, pdf_fun2 = pdf_fun2)
    }

    # result <- tryCatch(
    #   expr = {
    #     suppressWarnings(Rsolnp::solnp(
    #       pars = start_pars,
    #       fun = goal_fun,
    #       LB = LB,
    #       UB = UB,
    #       control = control,
    #       # ineqfun = ineq_fun,
    #       # ineqLB = c(0, 0, 0, rep(1e-6, p_ar)),
    #       # ineqUB = c(Inf, Inf, Inf, rep(Inf, p_ar)),
    #       rt = rt_core
    #     ))
    #   },
    #   error = function(e1) {
    #     stop("Convergence failed.", call. = FALSE)
    #   }
    # )

    result <- suppressWarnings(Rsolnp::solnp(
          pars = start_pars,
          fun = goal_fun,
          LB = LB,
          UB = UB,
          control = control,
          # ineqfun = ineq_fun,
          # ineqLB = c(0, 0, 0, rep(1e-6, p_ar)),
          # ineqUB = c(Inf, Inf, Inf, rep(Inf, p_ar)),
          rt = rt_core
        ))

  }

  pars <- result$pars

  c("gamma" = pars[[5]], "delta" = pars[[6]])

}

#'@export
#'@rdname show-methods
#'@aliases show,fEgarch_fit_aparch-method
setMethod(
  "show",
  "fEGarch_fit_aparch",
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
    "*************************************\n",
    "*        Fitted APARCH Model        *\n",
    "*************************************\n",
    " \n",
    "Type: apARCH\n",
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

