sigt_delta_fiaparch_R <- function(theta, x, p, q,
              incl_mean, extra_par, skew_par) {


  imean <- !incl_mean
  mu <- c(0, theta[[1]])[[incl_mean + 1]]
  step0 <- 2 - imean
  omega_sig <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q
  phi <- theta[(step0 + 1):step1]
  beta <- -theta[(step1 + 1):step2]
  gamma <- theta[[step2 + 1]]
  delta <- theta[[step2 + 2]]
  d <- theta[[step2 + 3]]
  if (extra_par) {
    shape <- theta[[step2 + 4]]
  } else {
    shape <- 0
  }
  if (skew_par) {
    skew <- theta[[step2 + 4 + extra_par]]
  } else {
    skew <- 0
  }

  n <- length(x)
  coef_inf <- -esemifar::farima_to_ar(ar = phi, ma = beta, d = d,
                                     max_i = n + 1000)[-1]

  sigt_delta <- c(sigt_delta_fiaparch_Cpp(
    x = x,
    mu = mu,
    coef_inf = coef_inf,
    gamma = gamma,
    delta = delta,
    omega_sig = omega_sig
  ))

  list(sigt_d = sigt_delta, mu = mu, skew = skew, shape = shape, d = d)

}

goal_fun_creator_fiaparch <- function(
    theta,
    x,
    pdf_fun2,
    p,
    q,
    incl_mean,
    extra_par,
    skew_par,
    penalty
) {

  res_list <- sigt_delta_fiaparch_R(theta, x, p, q,
              incl_mean, extra_par, skew_par)

  sigt_d <- res_list$sigt_d

  if (penalty && any(sigt_d <= 0)) {
    return(1e10)
  }

  d <- res_list$d
  sigt <- sigt_d^(1 / d)

  mu <- res_list$mu
  skew <- res_list$skew
  shape <- res_list$shape

  out <- -sum(log(pdf_fun2(x, mu, sigt, shape, skew)))

  out

}


#'FIAPARCH Model Fitting
#'
#'Fit a fractionally integrated APARCH model under the six
#'most common conditional distributions to observed data
#'using quasi maximum-likelihood estimation.
#'
#'@param rt the observed series ordered from past to present; can be
#'a numeric vector or a \code{"zoo"} class time series object.
#'@param orders a two-element numeric vector containing the two model
#'orders \eqn{p} and \eqn{q} (see Details for more information).
#'@param cond_dist the conditional distribution to consider as a
#'character object; the default is a conditional normal distribution
#'\code{"norm"}; available are also, however, a \eqn{t}-distribution
#'(\code{"std"}), a generalized error distribution (\code{"ged"}),
#'and their three skewed variants (\code{"snorm"}, \code{"sstd"},
#'\code{"sged"}).
#'@param include_mean a logical value indicating whether or not to include
#'the series' mean in the maximum-likelihood estimation process; for
#'\code{include_mean = FALSE}, the mean is fixed to zero.
#'@param start_pars the starting parameters for the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used.
#'@param LB the lower boundaries of the parameters in the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used.
#'@param UB the upper boundaries of the parameters in the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used.
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
#'See also the reference section for sources on the APARCH (Ding et al., 1993)
#'and the FIAPARCH (Tse, 1998) models.
#'
#'@export
#'
#'@return
#'An object of S4-class \code{"fEGarch_fit_fiaparch"}
#'is returned. It contains the following elements.
#'\describe{
#'\item{\code{pars}:}{a named numeric vector with the parameter estimates.}
#'\item{\code{se}:}{a named numeric vector with the obtained standard errors in accordance with the parameter estimates.}
#'\item{\code{rt}:}{the input object \code{rt}; if \code{rt} was a \code{"zoo"} object, the formatting is kept.}
#'\item{\code{sigt}:}{the estimated conditional standard deviations; if \code{rt} was a \code{"zoo"} object, the formatting is also applied to \code{sigt}.}
#'\item{\code{etat}:}{the obtained residuals; if \code{rt} was a \code{"zoo"} object, the formatting is also applied to \code{etat}.}
#'\item{\code{orders}:}{a two-element numeric vector stating the considered model orders.}
#'\item{\code{cond_dist}:}{a character value stating the conditional distribution considered in the model fitting.}
#'\item{\code{long_memo}:}{a logical value stating whether or not long memory was considered in the model fitting.}
#'\item{\code{llhood}:}{the log-likelihood value obtained at the optimal parameter combination.}
#'\item{\code{inf_criteria}:}{a named two-element numeric vector with the corresponding AIC (first element) and BIC (second element) of the fitted model.}
#'}
#'
#'@references
#'\itemize{
#'\item{Ding, Z., Granger, C. W. J., & Engle, R. F. (1993). A long memory property of stock market returns
#'and a new model. Journal of Empirical Finance, 1(1):
#'83-106. DOI: 10.1016/0927-5398(93)90006-D.}
#'\item{Tse, Y. K. (1998). The conditional heteroskedasticity of the
#'yen-dollar exchange rate. Journal of Applied Econometrics, 13(1):
#'49-55. DOI: 10.1002/(SICI)1099-1255(199801/02)13:1<49::AID-JAE459>3.0.CO;2-O.}
#'}
#'
fiaparch <- function(rt, orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
                        include_mean = TRUE,
                        start_pars = NULL, LB = NULL, UB = NULL) {

  rt_core <- zoo::coredata(rt)
  p <- orders[[1]]
  q <- orders[[2]]

  cond_d <- match.arg(cond_dist)

  pdf_fun2 <- lookup_table$fun2[[cond_d]]

  incl_mean = include_mean

  skew_par <- lookup_table[["skew_par"]][[cond_d]]
  extra_par <- lookup_table[["extra_par"]][[cond_d]]

  goal_fun <- function(theta, rt) {

    goal_fun_creator_fiaparch(
     theta = theta,
     x = rt,
     pdf_fun2 = pdf_fun2,
     p = p,
     q = q,
     incl_mean = incl_mean,
     extra_par = extra_par,
     skew_par = skew_par,
     penalty = TRUE
    )

  }

  # goal_fun2 <- function(theta, rt) {
  #
  #   goal_fun_creator_fiaparch(
  #    theta = theta,
  #    x = rt,
  #    pdf_fun2 = pdf_fun2,
  #    p = p,
  #    q = q,
  #    incl_mean = incl_mean,
  #    extra_par = extra_par,
  #    skew_par = skew_par,
  #    penalty = FALSE
  #   )
  #
  # }

  vals <- lookup_table[["vals"]][[cond_d]]

  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  d_spar <- 0.25
  d_low <- 0
  d_up <- 0.4999
  d_name <- "d"

  if (is.null(start_pars)) {

    start_pars <- c(
      list(NULL, 0)[[incl_mean + 1]],
      0.01,
      rep(0.2 / p, p),
      rep(0.4 / q, q),
      0.2, 2,
      d_spar,
      sval,
      list(NULL, 1)[[skew_par + 1]]
    )

  }

  if (is.null(LB)) {

    LB <- c(
      list(NULL, -2)[[incl_mean + 1]],
      1e-10,
      rep(1e-6, p),
      rep(1e-6, q),
      -0.9999, 1,
      d_low,
      lval,
      list(NULL, 1e-15)[[skew_par + 1]]
    )

  }

  if (is.null(UB)) {

    UB <- c(
      list(NULL, 2)[[incl_mean + 1]],
      2,
      rep(1 - 1e-6, p),
      rep(1 - 1e-6, q),
      0.9999, 5,
      d_up,
      uval,
      list(NULL, Inf)[[skew_par + 1]]
    )

  }


  result <- tryCatch(
    expr = {suppressWarnings(Rsolnp::solnp(
    pars = start_pars,
    fun = goal_fun,
    LB = LB,
    UB = UB,
    control = list(trace = TRUE, tol = 1e-6),
    rt = rt_core
  ))
  }, error = function(e1) {stop("Convergence failed.")})

  pars <- result$pars

  # ineq_fun <- function(theta, rt) {
  #   imean <- !incl_mean
  #   step0 <- 2 - imean
  #   step1 <- step0 + p
  #   step2 <- step1 + q
  #   phi <- theta[(step0 + 1):step1]
  #   beta <- -theta[(step1 + 1):step2]
  #
  #   if (p > q) {
  #     beta <- c(beta, rep(0, p - q))
  #   } else if (q > p) {
  #     phi <- c(phi, rep(0, q - p))
  #   }
  #
  #   beta - phi
  # }
  #
  # result <- suppressWarnings(Rsolnp::solnp(
  #   pars = pars,
  #   fun = goal_fun2,
  #   LB = LB,
  #   UB = UB,
  #   control = list(trace = FALSE, tol = 1e-6),
  #   ineqfun = ineq_fun,
  #   ineqLB = rep(0, max(p, q)),
  #   ineqUB = rep(0.9999, max(p, q)),
  #   rt = rt_core
  # ))
  #
  # pars <- result$pars

  llhood <- -tail(result$values, 1)

  res_f <- sigt_delta_fiaparch_R(unname(pars), rt_core, p, q,
              incl_mean, extra_par, skew_par)
  sigt_d <- res_f$sigt_d
  d_opt <- res_f$d
  sigt <- sigt_d^(1 / d_opt)

  mu <- c(0, pars[[1]])[[incl_mean + 1]]
  etat <- (rt - mu) / sigt

  beta_names <- paste0("beta", 1:q)

  par_names <- c(
    list(NULL, "mu")[[incl_mean + 1]],
    "omega",
    paste0("phi", 1:p),
    beta_names,
    "gamma",
    "delta",
    d_name,
    list(NULL, c("shape", "df")[[(cond_d %in% c("std", "sstd")) + 1]])[[extra_par + 1]],
    list(NULL, "skew")[[skew_par + 1]]
  )

  serrors <- tryCatch(
    expr = {sqrt(diag(solve(result$hessian)))},
    error = function(e1) {
      warning("Unable to obtain standard errors.")
      rep(NA, length(pars))
    }
  )
  names(pars) <- par_names
  names(serrors) <- par_names

  k <- length(pars)
  n <- length(rt)

  m2llhood <- -2 * llhood
  bic <- (m2llhood + k * log(n)) / n
  aic <- (m2llhood + 2 * k) / n

  if (inherits(rt, what = "zoo")) {
    tp <- as.Date(time(rt))
    sigt <- zoo::zoo(zoo::coredata(sigt), order.by = tp)
    etat <- zoo::zoo(zoo::coredata(etat), order.by = tp)
  }

  list_out <- fEGarch_fit_fiaparch(
    pars = pars,
    se = serrors,
    rt = rt,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = c("aic" = aic, "bic" = bic),
    cond_dist = cond_d,
    orders = orders,
    long_memo = TRUE
  )

  list_out

}

#'@export
setMethod(
  "show",
  "fEGarch_fit_fiaparch",
  function(object) {
  x <- object
  par_names <- names(x@pars)
  se <- unname(x@se)
  pars <- unname(x@pars)

  if (all(!is.na(se))) {
    tvals <- pars / se
    atvals <- abs(tvals)

    pvals <- 2 * (1 - pnorm(atvals))
  } else {

    tvals <- pvals <- rep(NA, length(se))

  }


  df <- data.frame(
    par = sprintf("%.4f", pars),
    se = sprintf("%.4f", se),
    tval = sprintf("%.4f", tvals),
    pval = sprintf("%.4f", pvals)
  )
  row.names(df) <- par_names

  part1 <- paste0(
    "*************************************\n",
    "*       Fitted FIAPARCH Model       *\n",
    "*************************************\n",
    " \n",
    "Type: fiapARCH\n",
    "Orders: (", x@orders[[1]], ", ", x@orders[[2]], ")\n",
    "Long memory: ", x@long_memo, "\n",
    "Cond. distribution: ", x@cond_dist, "\n",
    " \n",
    "Fitted parameters:\n"
  )

  part2 <- paste0(
    " \n",
    "Information criteria:\n",
    "AIC: ", sprintf("%.4f", unname(x@inf_criteria[[1]])),
    ", BIC: ", sprintf("%.4f", unname(x@inf_criteria[[2]]))
  )

  cat(part1, fill = TRUE)
  print(df)
  cat(part2, fill = TRUE)
  }
)
