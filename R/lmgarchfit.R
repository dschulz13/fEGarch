sigt_lmgarch_R <- function(theta, x, p, q,
              incl_mean, extra_par, skew_par) {


  imean <- !incl_mean
  mu <- c(0, theta[[1]])[[incl_mean + 1]]
  step0 <- 2 - imean
  step1 <- step0 + p
  step2 <- step1 + q
  omega <- theta[[step0]]
  phi <- theta[(step0 + 1):step1]
  beta <- -theta[(step1 + 1):step2]
  d <- theta[[step2 + 1]]
  if (extra_par) {
    shape <- theta[[step2 + 2]]
  } else {
    shape <- 0
  }
  if (skew_par) {
    skew <- theta[[step2 + 2 + extra_par]]
  } else {
    skew <- 0
  }

  n <- length(x)
  coef_inf <- esemifar::farima_to_ar(ar = phi, ma = beta, d = d,
                                     max_i = n + 1000)[-1]

  # For initialization
  E_e2 <- var(x)

  sigt <- c(sigt_lmgarch_Cpp(
    x = x,
    mu = mu,
    omega,
    coef_inf = coef_inf,
    E_e2
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape)

}

goal_fun_creator_lmgarch <- function(
    theta,
    x,
    pdf_fun2,
    p,
    q,
    incl_mean,
    extra_par,
    skew_par
) {

  res_list <- sigt_lmgarch_R(theta, x, p, q,
              incl_mean, extra_par, skew_par)

  sigt <- res_list$sigt
  mu <- res_list$mu
  skew <- res_list$skew
  shape <- res_list$shape

  out <- -sum(log(pdf_fun2(x, rep(mu, length(sigt)), sigt, shape, skew)))

  out

}


#'LMGARCH Model Fitting
#'
#'Fit a long-memory GARCH (LMGARCH) model under the six
#'most common conditional distributions to observed data
#'using quasi maximum-likelihood estimation.
#'
#'@param rt the observed series ordered from past to present; can be
#'a numeric vector or a \code{"zoo"} class time series object.
#'@param orders a two-element numeric vector containing the two model
#'orders \eqn{p} and \eqn{q} (see Details for more information); currently,
#'only the default \code{orders = c(1, 1)} is supported; other specifications
#'of a two-element numeric vector will lead to \code{orders = c(1, 1)} being
#'run and a warning message being returned.
#'@param cond_dist the conditional distribution to consider as a
#'character object; the default is a conditional normal distribution
#'\code{"norm"}; available are also, however, a \eqn{t}-distribution
#'(\code{"std"}), a generalized error distribution (\code{"ged"}),
#'an average Laplace distribution (\code{"ald"}),
#'and their four skewed variants (\code{"snorm"}, \code{"sstd"},
#'\code{"sged"}, \code{"sald"}).
#'@param include_mean a logical value indicating whether or not to include
#'the series' mean in the maximum-likelihood estimation process; for
#'\code{include_mean = FALSE}, the mean is fixed to zero.
#'@param drange a two-element numeric vector that gives the boundaries of the
#'search interval for the fractional differencing parameter \eqn{d}; is
#'overwritten by the settings of the arguments \code{LB} and \code{UB}.
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
#'@param control a list that is passed to \code{control} of the
#'function \code{solnp} of the package \code{Rsolnp}.
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
#'\deqn{\sigma_t^{2}=\omega+\left[1-\beta^{-1}(B)\phi(B)(1-B)^{d}\right]\varepsilon_t^2.}
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
#'real-valued parameter and \eqn{d \in [0,0.5)} as the
#'long-memory parameter. \eqn{\omega > 0} is the intercept. It is assumed that
#'all \eqn{\beta_j} and \eqn{\phi_i} are non-negative. Furthermore, we have
#'\eqn{E(\varepsilon_t^2)= \omega / \left(1-\psi(1)\right)}, where
#'\deqn{\psi(B) = \left[1-\beta^{-1}(B)\phi(B)(1-B)^{d}\right].}
#'
#'The pre-sample values of \eqn{\varepsilon_t^2} are taken as the sample
#'variance of the observed series. The pre-sample is set to length 1000, and the infinite coefficient series
#'\eqn{\psi(B)} is, at each time point, always considered
#'as far back as to the first pre-sample value. Therefore, at observation time point
#'\eqn{t=1,\dots,n}, a truncation of 1000 + t - 1 coefficients is considered.
#'
#'Currently, only a model of orders \eqn{p=1} with \eqn{q=1} can be fitted;
#'to ensure the non-negativity of all of the infinite-order coefficient series
#'\eqn{\psi(B)}, which in combination
#'with \eqn{\omega>0} ensures that all the conditional volatilities
#'are greater than zero, we employ the sufficient conditions mentioned
#'in Bollerslev and Mikkelsen (1996). Conrad and Haag (2006) provide
#'necessary and sufficient conditions for cases with simple model orders and
#'sufficient conditions for general model orders, however, due to
#'time efficiency, this function uses the set of sufficient conditions
#'by Bollerslev and Mikkelsen (1996) only.
#'
#'See also the reference section for sources on the LMGARCH model
#'(Karanasos et al., 2004; Conrad and Karanasos, 2006; Conrad and
#'Haag, 2006).
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
#'@export
#'
#'@return
#'An object of S4-class \code{"fEGarch_fit_lmgarch"}
#'is returned. It contains the following elements.
#'\describe{
#'\item{\code{pars}:}{a named numeric vector with the parameter estimates.}
#'\item{\code{se}:}{a named numeric vector with the obtained standard errors in accordance with the parameter estimates.}
#'\item{\code{vcov_mat}:}{the variance-covariance matrix of the parameter estimates with named columns and rows.}
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
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- lmgarch(rt)
#'model
#'
#'@importFrom stats var
#'
#'@references
#'\itemize{
#'\item{Bollerslev, T., & Mikkelsen, H. O. (1996). Modeling and pricing long
#'memory in stock market volatility. Journal of Econometrics, 73(1):
#'151-184. DOI: 10.1016/0304-4076(95)01736-4.}
#'\item{Conrad, C., & Haag, B. R. (2006). Inequality constraints in the fractionally
#'integrated GARCH model. Journal of Financial Econometrics, 4(3):
#'413-449. DOI: 10.1093/jjfinec/nbj015.}
#'\item{Conrad, C., & Karanasos, M. (2006). The impulse response function
#'of the long memory GARCH process. Economics Letters, 90(1):
#'34-41. DOI: 10.1016/j.econlet.2005.07.001.}
#'\item{Karanasos, M., Psaradakis, Z., & Sola, M. (2004). On the autocorrelation properties of
#'long-memory GARCH processes. Journal of Time Series Analysis, 25(2):
#'265-281. DOI: 10.1046/j.0143-9782.2003.00349.x.}
#'}
#'
lmgarch <- function(rt, orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                        include_mean = TRUE, drange = c(0, 0.5 - 1e-6),
                        start_pars = NULL, LB = NULL, UB = NULL, control = list(),
                        parallel = TRUE, ncores = max(1, future::availableCores() - 1)) {

  if (!all(orders == c(1, 1))) {
    orders <- c(1, 1)
    message("Only orders p=1 with q=1 are currently supported for FIAPARCH models.")
  }

  rt_core <- zoo::coredata(rt)
  p <- orders[[1]]
  q <- orders[[2]]

  cond_d <- match.arg(cond_dist)

  pdf_fun2 <- lookup_table$fun2[[cond_d]]

  incl_mean <- include_mean

  skew_par <- lookup_table[["skew_par"]][[cond_d]]
  extra_par <- lookup_table[["extra_par"]][[cond_d]]

  goal_fun_s <- function(theta, rt, pdf_fun2) {

    goal_fun_creator_lmgarch(
     theta = theta,
     x = rt,
     pdf_fun2 = pdf_fun2,
     p = p,
     q = q,
     incl_mean = incl_mean,
     extra_par,
     skew_par
    )

  }

  vals <- lookup_table[["vals"]][[cond_d]]

  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  d_low <- drange[[1]]
  d_up <- drange[[2]]
  d_spar <- d_low + abs(d_up - d_low) / 2
  d_name <- "d"

  if (is.null(start_pars)) {

    start_pars <- c(
      list(NULL, 0)[[incl_mean + 1]],
      0.001,
      rep(0.05 / p, p),
      rep(0.4 / q, q),
      d_spar,
      sval,
      list(NULL, 1)[[skew_par + 1]]
    )

  }

  if (is.null(LB)) {

    LB <- c(
      list(NULL, -2)[[incl_mean + 1]],
      1e-10,
      rep(0, p),
      rep(0, q),
      d_low,
      lval,
      list(NULL, 1e-15)[[skew_par + 1]]
    )

  }

  if (is.null(UB)) {

    UB <- c(
      list(NULL, 2)[[incl_mean + 1]],
      Inf,
      rep(2, p),
      rep(2, q),
      d_up,
      uval,
      list(NULL, Inf)[[skew_par + 1]]
    )

  }

  ineq_fun <- function(theta, rt) {

    imean <- !incl_mean
    step0 <- 2 - imean
    step1 <- step0 + p
    step2 <- step1 + q
    omega <- theta[[step0]]
    phi <- theta[(step0 + 1):step1]
    beta <- theta[(step1 + 1):step2]
    d <- theta[[step2 + 1]]

    r1 <- phi - beta + d
    r2 <- (2 - d) / 3 - phi
    r3 <- beta * r1 - d * (phi - (1 - d) / 2)

    c(r1, r2, r3)

  }

  if (is.null(control$trace)) {control$trace <- FALSE}

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
      .f = function(.x, start_pars, goal_fun_s, LB, UB, ineq_fun, control, rt_core, dfun2) {

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
          ineqfun = ineq_fun,
          ineqLB = c(0, 0, 0),
          ineqUB = c(Inf, Inf, Inf),
          control = control,
          rt = rt_core
        ))
        est
      }, start_pars = start_pars, goal_fun_s = goal_fun_s, LB = LB, UB = UB,
      ineq_fun = ineq_fun, control = control, rt_core = rt_core,
      dfun2 = pdf_fun2,
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

    result <- tryCatch(
      expr = {
        suppressWarnings(Rsolnp::solnp(
          pars = start_pars,
          fun = goal_fun,
          LB = LB,
          UB = UB,
          control = control,
          ineqfun = ineq_fun,
          ineqLB = c(0, 0, 0),
          ineqUB = c(Inf, Inf, Inf),
          rt = rt_core
        ))
      },
      error = function(e1) {
        stop("Convergence failed.", call. = FALSE)
      }
    )

  }

  pars <- result$pars
  llhood <- -tail(result$values, 1)

  sigt <- sigt_lmgarch_R(unname(pars), rt_core, p, q,
              incl_mean, extra_par, skew_par)$sigt

  mu <- c(0, pars[[1]])[[incl_mean + 1]]
  etat <- (rt - mu) / sigt

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
    "omega",
    paste0("phi", 1:p),
    paste0("beta", 1:q),
    d_name,
    extra_par_name,
    list(NULL, "skew")[[skew_par + 1]]
  )

  len_par <- length(pars)
  hessian <- result$hessian[1:len_par, 1:len_par]
  vcov_mat <- solve(hessian)

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

  serrors <- sqrt(diag(vcov_mat))
  names(pars) <- par_names
  names(serrors) <- par_names

  colnames(vcov_mat) <- par_names
  rownames(vcov_mat) <- par_names

  k <- length(pars)
  n <- length(rt)

  m2llhood <- -2 * llhood
  bic <- (m2llhood + k * log(n)) / n
  aic <- (m2llhood + 2 * k) / n

  if (inherits(rt, what = "zoo")) {
    tp <- as.Date(time(rt))
    sigt <- zoo::zoo(sigt, order.by = tp)
    etat <- zoo::zoo(etat, order.by = tp)
  }

  list_out <- fEGarch_fit_lmgarch(
    pars = pars,
    se = serrors,
    vcov_mat = vcov_mat,
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
#'@rdname show-methods
#'@aliases show,fEgarch_fit_lmgarch-method
setMethod(
  "show",
  "fEGarch_fit_lmgarch",
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

  part1 <- paste0(
    "*************************************\n",
    "*        Fitted LMGARCH Model       *\n",
    "*************************************\n",
    " \n",
    "Type: lmGARCH\n",
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

