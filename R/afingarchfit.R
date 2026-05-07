afingarch_grab_pars <- function(theta, p, q, incl_mean, extra_par, skew_par) {

  # The same as for FINGARCH (so just call that function)
  fingarch_grab_pars(theta = theta, p = p, q = q, incl_mean = incl_mean,
                     extra_par = extra_par, skew_par = skew_par)

}

sigt_afingarch_R <- function(theta, x, p, q, p_ar, q_ma,
              incl_mean, extra_par, skew_par,
              trunc, presample, ...) {

  cur_env <- environment()
  list2env(list(...), envir = cur_env)
  stopifnot("No object 'dens_fun' passed to conditional standard deviation function." =
              (exists("dens_fun", envir = cur_env, inherits = FALSE)))

  all_pars <- afingarch_grab_pars(
    theta = theta,
    p = p,
    q = q,
    incl_mean = incl_mean,
    extra_par = extra_par,
    skew_par = skew_par
  )
  mu <- all_pars$mu
  omega_sig <- all_pars$omega_sig
  phi <- all_pars$phi
  beta <- all_pars$beta
  gamma <- rep(0, p)       # fixed for NGARCH-type
  delta <- all_pars$delta
  d <- all_pars$d
  shape <- all_pars$shape
  skew <- all_pars$skew

  n <- length(x) + presample

  coef_inf <- ar_infty(
    ar = phi,
    ma = beta,
    d = d,
    max_i = trunc
  )
  coef_inf[[1]] <- 0

  x_dm <- x - mu

  Czd <- comp_C_zd(
    shape = shape,
    skew = skew,
    gamma = gamma,
    delta = delta,
    dens_fun = dens_fun   # !!! Check later that this works out
  )

  if (is.na(Czd) || is.nan(Czd)) {
    message("Czd is NA or NaN on some iterations.")
  }

  z_sc <- afiaparch_adjust_returns(
    rt = x_dm,
    gamma = gamma,
    delta = delta,
    Czd = Czd
  )

  # if (!exists("presample_val", inherits = FALSE)) {
  #   presample_val <- mean(z_sc)
  # }

  lc <- length(coef_inf)
  np2 <- nextpow2(2 * n - 1)
  c_adj <- c(coef_inf, rep(0, np2 - lc))
  x_adj <- c(rep(presample_val, presample), z_sc, rep(0, np2 - n))

  conv_raw <- Re(
    stats::fft(
      stats::fft(c_adj) * stats::fft(x_adj),
      inverse = TRUE
    )
  )

  conv_part <- utils::tail(
    utils::head(conv_raw, n),
    length(x_dm)
  ) / np2

  # inside_sc <- omega_sig + conv_part

  # if (any(!is.finite(inside_sc)) || any(inside_sc <= 0)) {
  #   return(list(
  #     sigt = rep(NaN, length(x_dm)),
  #     mu = mu,
  #     skew = skew,
  #     shape = shape,
  #     cmeans = rep(mu, length(x_dm)),
  #     inside_sc = inside_sc,
  #     z_sc = z_sc,
  #     coef_inf = coef_inf,
  #     Czd = Czd,
  #     presample_val = presample_val
  #   ))
  # }

  list(
    sigt = (omega_sig + conv_part)^(1 / delta),
    mu = mu,
    skew = skew,
    shape = shape,
    cmeans = rep(mu, length(conv_part)),
    # inside_sc = inside_sc,
    # z_sc = z_sc,
    # coef_inf = coef_inf,
    Czd = Czd   #,        !!! Check later, if we need to return this
#    presample_val = presample_val
  )

}

goal_fun_creator_afingarch <- function(
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
    penalty,
    sig_fun,
    trunc,
    presample,
    presample_val,
    dens_fun
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
    presample = presample,
    presample_val = presample_val,
    dens_fun = dens_fun)

  nllhood_calc(res_list, x, pdf_fun2)

}


#'AFINGARCH Model Fitting
#'
#'Fit an adjusted fractionally integrated NGARCH model under the six
#'most common and further conditional distributions to observed data
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
#'an internally saved default set of values is used; see "Details" for the
#'order of elements; elements should be set with respect to a series rescaled
#'to have sample variance one.
#'@param drange a two-element numeric vector that gives the boundaries of the
#'search interval for the fractional differencing parameter \eqn{d} in the
#'conditional volatility model part; is
#'overwritten by the settings of the arguments \code{LB} and \code{UB}.
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
#'by \code{\link{locpol_spec}}; defines the settings of the nonparametric
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
#'an internally saved default set of values is used; see "Details" for the
#'order of elements; elements should be set with respect to a series rescaled
#'to have sample variance one.
#'@param UB the upper boundaries of the parameters in the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used; see "Details" for the
#'order of elements; elements should be set with respect to a series rescaled
#'to have sample variance one.
#'@param control a list that is passed to \code{control} of the
#'function \code{solnp} of the package \code{Rsolnp}.
#'@param control_nonpar a list containing changes to the arguments
#'for the hyperparameter estimation algorithm in the nonparametric
#'scale function estimation for
#'\code{use_nonpar = TRUE}; see "Details" for more information.
#'@param mean_after_nonpar only for \code{use_nonpar = TRUE}; considers the unconditional mean
#'of the parametric model part in the QMLE step in a semiparametric model; by default, a zero-mean
#'model is considered for the parametric part in a semiparametric model.
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
#'@param trunc a positive integer indicating the finite truncation length of the
#'infinite-order polynomials of the infinite-order representations of the
#'long-memory model parts; the character \code{"none"} is an optional input
#'that specifies that truncation should always be applied back to the first (presample) observation
#'time point, i.e. that maximum length filters should be applied at all times.
#'@param presample the presample length for initialization (for extended EGARCH- / Log-GARCH-type
#'models only relevant for the FARIMA-part, as series in log-transformed
#'conditional variance are initialized by zero).
#'@param Prange a two-element vector that indicates the search boundaries for
#'the parameter \eqn{P} in a (skewed) average Laplace distribution.
#'@param skip_vcov a logical indicating whether or not to skip the computation
#'of the variance-covariance matrix of the parameter estimators and therefore
#'also standard error computation.
#'
#'@details
#'Consider an adjusted FINGARCH(\eqn{p, d, q}), an AFINGARCH(\eqn{p, d, q}) for short,
#'with constant asymmetry term zero regardless
#'of the lag. The model was introduced by Feng et al. (2026) and is an adjusted version
#'of the fractionally integrated NGARCH (FINGARCH). Models from the original FIAPARCH class are often
#'not well defined overall, and they are usually not weakly stationary
#'(Feng et al., 2026). Hence, the adjusted FIAPARCH class is recommended.
#'
#'An AFINGARCH(\eqn{p, d, q}) is defined as follows. Let \eqn{\left\{r_t\right\}}, with \eqn{t \in \mathbb{Z}} as the
#'time index, be a theoretical time series that follows
#'\deqn{r_t=\mu+\varepsilon_t \text{ with } \varepsilon_t=\sigma_t \eta_t \text{ and } \eta_t \sim \text{IID}(0,1), \text{ where}}
#'\deqn{\sigma_t^{\delta}=\left[1-\beta^{-1}(B)\phi(B)(1-B)^{d}\right]C_z^{-1}\left| \varepsilon_t \right|^{\delta}.}
#'Here, \eqn{\eta_t\sim\text{IID}(0,1)} means that the innovations
#'\eqn{\eta_t} are independent and identically distributed (iid) innovations with mean zero
#'and variance one, whereas \eqn{\sigma_t > 0} are the conditional standard
#'deviations in \eqn{r_t}. \eqn{S_t=1} for \eqn{\varepsilon_t < 0} and \eqn{S_t=-1} for
#'\eqn{\varepsilon_t \geq 0}, and
#'\deqn{C_z=p_n \times  E\left(|\eta_t|^{\delta} \mid \eta_t < 0\right) + p_p \times  E\left(|\eta_t|^{\delta} \mid \eta_t \geq 0\right),}
#'where \eqn{p_n = \text{Pr}(S_t = 1)} and \eqn{p_p = 1 - p_n}. As we can see, the theoretical
#'process has no volatility intercept. However, since in practice a finite-length truncation is
#'considered for the infinite-length polynomial \eqn{\left[1-\beta^{-1}(B)\phi(B)(1-B)^{d}\right]},
#'the practical implementation does indeed consider a volatility intercept \eqn{\omega > 0}. While
#'\eqn{\omega} is in theory not a free parameter then, we still simply include it as a parameter
#'to optimize over in the estimation routine.
#'Moreover, \eqn{B} is the backshift operator and
#'\eqn{\beta(B) = 1 - \sum_{j=1}^{q}\beta_j B^{j}}, where
#'\eqn{\beta_j}, \eqn{j=1,2,\dots, q}, are real-valued coefficients. Furthermore,
#'\eqn{\phi(B) = 1 - \sum_{i=1}^{p}\phi_i B^{i}}, where
#'\eqn{\phi_i}, \eqn{i=1,2,\dots, p}, are real-valued coefficients. \eqn{p}
#'and \eqn{q} are the model orders definable through the argument \code{orders},
#'where \eqn{p} is the first element and \eqn{q} is the second element in the
#'argument. In addition, we have \eqn{\mu = E\left(r_t\right)} as a
#'real-valued parameter, \eqn{\gamma \in (-1,1)} and \eqn{d \in [0,1]} as the
#'parameter for the level of integration. With \eqn{d = 0} the model reduces
#'to a short-memory model, for \eqn{d=1} we have a full integration, and for
#'\eqn{d\in(0, 1)}, we have fractional integration, where \eqn{d\in(0, 0.5)} is
#'considered to describe a long-memory process. \eqn{\omega > 0} is the intercept.
#'It is assumed that
#'all elements of the infinite-length polynomial are non-negative. For further information
#'on this model and its theoretical advantages, we refer the reader to
#'Feng et al. (2026).
#'
#'The pre-sample values of \eqn{C_{z}^{-1}\left| \varepsilon_t \right|^{\delta}}
#'are replaced by the in-sample arithmetic mean of
#'\deqn{C_{z,0}^{-1}\left[\left| r_t^{*} -  \bar{r}_t^{*}\right|\right]^{\delta_0}}
#'where
#'\eqn{r_t^{*}} are considered as the observations and \eqn{\bar{r}_t^{*}} is the sample
#'mean of \eqn{r_t^{*}}. An initial value
#'for
#'\eqn{\delta_0} as an initial guess for \eqn{\delta} is obtained
#'from a previous fitting of a short-memory NGARCH. Analogously, an initial value
#'of \eqn{C_z} denoted by \eqn{C_{z,0}} is also obtained from that initial fit by
#'replacing \eqn{\delta} in the formula of \eqn{C_z} by
#'\eqn{\delta_0}, and by obtaining \eqn{p_n E(|\eta_t|^{\delta} \mid \eta_t < 0)}
#'and \eqn{p_p E(|\eta_t|^{\delta} \mid \eta_t \geq 0)} numerically via the
#'\code{integrate} function using the fitted
#'parameters from the initial fit.
#'
#'The truncated infinite order polynomial is computed following the idea
#'by Nielsen
#'and Noel (2021) as is the series of conditional variances for most
#'computational efficiency. To ensure stability of the first fitted in-sample
#'conditional standard deviations, we however use a small, but also adjustable (also
#'to length zero) presample, which may introduce biases into the parameter estimators.
#'
#'Currently, only a model of orders \eqn{p=1} with \eqn{q=1} can be fitted;
#'to ensure the non-negativity of all of the infinite-order coefficient series
#'\eqn{\left[1-\beta^{-1}(B)\phi(B)(1-B)^{d}\right]}, which in combination
#'with \eqn{\omega>0} ensures that all the conditional volatilities
#'are greater than zero, we employ inequality constraints ensuring that the first
#'50 coefficients of the infinite-order ARCH-representation are non-negative as an
#'approximation to ensuring that all of the coefficients are non-negative. To ensure
#'that they are non-negative, one may in theory consider the sufficient conditions
#'mentioned in Bollerslev and Mikkelsen (1996) or Tse (1998), which are however
#'sometimes restrictive, or the simultaneously necessary and sufficient conditions
#'by Conrad and Haag (2006), which are however complex to implement properly.
#'
#'See also the reference section for sources on the
#'NGARCH model (Higgins and Bera, 1992), the more general
#'APARCH (Ding et al., 1993) and the FIAPARCH (Tse, 1998) models.
#'
#'In the current package version, standard errors of parameter estimates are
#'computed from the Hessian at the optimum of the log-likelihood using
#'\code{\link[numDeriv]{hessian}}. To ensure numerical stability and
#'applicability to a huge variety of differently scaled data, parametric
#'models are first fitted to data that is scaled to have sample variance
#'1. Parameter estimates and other quantities are then either
#'retransformed or recalculated afterwards for the original data.
#'
#'For a conditional average Laplace distribution, an optimal model for each
#'distribution parameter \eqn{P} from 1 to 5 is estimated (assuming that
#'\eqn{P} is then fixed to the corresponding value). Afterwards, \eqn{P} is then
#'estimated by selecting the estimated model among the five fitted models that
#'has the largest log-likelihood. The five models are, by default, fitted
#'simultaneously using parallel programming techniques (see also the arguments
#'\code{parallel} and \code{ncores}, which are only relevant for a conditional
#'average Laplace distribution). After the optimal model (including
#'the estimate of \eqn{P} called \eqn{\hat{P}}) has been determined, \eqn{P=\hat{P}}
#'is seen as fixed to obtain the standard errors via the Hessian matrix for the
#'estimates of the continuous parameters. A standard error for \eqn{\hat{P}} is therefore
#'not obtained and the ones obtained for the remaining estimates do not account
#'for \eqn{\hat{P}}.
#'
#'An ARMA-AFINGARCH or a FARIMA-AFINGARCH can be fitted by adjusting the
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
#'The order for manual settings of \code{start_pars}, \code{LB} and \code{UB}
#'is crucial. The correct order is: \eqn{\mu}, \eqn{\text{ar}_1,\dots,\text{ar}_{p^{*}}},
#'\eqn{\text{ma}_1,\dots,\text{ma}_{q^{*}}},\eqn{D},\eqn{\omega},
#'\eqn{\phi}, \eqn{\beta}, \eqn{\delta}, \eqn{d},
#'shape parameter,
#'skewness parameter. Depending on the exact model specification,
#'parameters irrelevant for the specification at hand should be dropped
#'in \code{start_pars}, \code{LB} and \code{UB}.
#'
#'@export
#'
#'@return
#'An object of S4-class \code{"fEGarch_fit_afingarch"}
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
#'\item{\code{C_z}:}{the obtained constant \eqn{C_z} for the estimated model; helpful in
#'analyzing further properties of the fitted model.}
#'}
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- afingarch(rt)
#'model
#'
#'@references
#'\itemize{
#'\item{Ayensu, O. K., Feng, Y., & Schulz, D. (2026). Well-known and recent
#'long-memory GARCH models and their semiparametric extensions. Forthcoming
#'preprint, Paderborn University.}
#'\item{Bollerslev, T., & Mikkelsen, H. O. (1996). Modeling and pricing long
#'memory in stock market volatility. Journal of Econometrics, 73(1),
#'151-184. DOI: 10.1016/0304-4076(95)01736-4.}
#'\item{Conrad, C., & Haag, B. R. (2006). Inequality constraints in the fractionally
#'integrated GARCH model. Journal of Financial Econometrics, 4(3),
#'413-449. DOI: 10.1093/jjfinec/nbj015.}
#'\item{Ding, Z., Granger, C. W. J., & Engle, R. F. (1993). A long memory property of stock market returns
#'and a new model. Journal of Empirical Finance, 1(1),
#'83-106. DOI: 10.1016/0927-5398(93)90006-D.}
#'\item{Feng, Y., Gries, T., Letmathe, S., & Schulz, D. (2022). The smoots Package in R for Semiparametric Modeling of
#'Trend Stationary Time Series. The R Journal,
#'14(1), 182-195. URL: https://journal.r-project.org/articles/RJ-2022-017/.}
#'\item{Higgins, M. L., & Bera, A. K. (1992). A Class of Nonlinear ARCH Models.
#'International Economic Review, 33(1), 137–158. DOI: 10.2307/2526988.}
#'\item{Letmathe, S., Beran, J., & Feng, Y. (2023). An extended exponential SEMIFAR model with application
#'in R. Communications in Statistics - Theory and Methods,
#'53(22), 7914–7926. DOI: 10.1080/03610926.2023.2276049.}
#'\item{Nielsen, M. O., & Noel, A. L. (2021). To infinity and beyond: Efficient
#'computation of ARCH(\eqn{\infty}) models. Journal of Time Series Analysis,
#'42(3), 338-354. DOI: 10.1111/jtsa.12570.}
#'\item{Tse, Y. K. (1998). The conditional heteroskedasticity of the
#'yen-dollar exchange rate. Journal of Applied Econometrics, 13(1),
#'49-55. DOI: 10.1002/(SICI)1099-1255(199801/02)13:1<49::AID-JAE459>3.0.CO;2-O.}
#'}
#'
afingarch <- function(rt, orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                        drange = c(0, 1),
                        meanspec = mean_spec(), Drange = c(0, 1),
                        nonparspec = locpol_spec(), use_nonpar = FALSE,
                        n_test = 0,
                        start_pars = NULL, LB = NULL, UB = NULL, control = list(),
                        control_nonpar = list(), mean_after_nonpar = FALSE,
                        parallel = TRUE, ncores = max(1, future::availableCores() - 1),
                        trunc = "none", presample = 50, Prange = c(1, 5),
                        skip_vcov = FALSE) {

  inp_args <- mget(names(formals()), envir = environment())

  inp_args[["model_type"]] <- "afingarch"
  do.call(general_garch_fitting, args = inp_args)

}

#'@export
#'@rdname show-methods
#'@aliases show,fEgarch_fit_afingarch-method
setMethod(
  "show",
  "fEGarch_fit_afingarch",
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
    "*       Fitted AFINGARCH Model      *\n",
    "*************************************\n",
    " \n",
    "Type: afingarch\n",
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

  part1_2 <- paste0(
    " \n",
    "Scaling constant at optimal solution:\n",
    "Cz: ", sprintf("%.4f", unname(x@Cz))
  )

  part2 <- paste0(
    " \n",
    "Information criteria (parametric part):\n",
    "LLH: ", sprintf("%.4f", unname(x@llhood) / length(x@rt)),
    ", AIC: ", sprintf("%.4f", unname(x@inf_criteria[[1]])),
    ", BIC: ", sprintf("%.4f", unname(x@inf_criteria[[2]]))
  )

  cat(part1, fill = TRUE)
  print(df)
  cat(part1_2, fill = TRUE)
  cat(part2, fill = TRUE)
  }
)
