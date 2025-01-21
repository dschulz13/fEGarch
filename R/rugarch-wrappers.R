rugarch_modelfit <- function(rt, model = "apARCH", orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
                        include_mean = TRUE,
                        start_pars = NULL) {

  rt_core <- zoo::coredata(rt)

  cond_dist <- match.arg(cond_dist)

  if (is.null(start_pars)) {
    start.pars <- list()
  } else {
    start.pars <- start_pars
  }

  if (model == "gjrGARCH") {
    fixed.pars <- list(delta = 2)
    model_s <- "apARCH"
  } else {
    fixed.pars <- list()
    model_s <- model
  }

  spec <- rugarch::ugarchspec(
    variance.model = list(model = model_s, garchOrder = orders),
    mean.model = list(armaOrder = c(0, 0), include.mean = include_mean),
    distribution.model = cond_dist, start.pars = start.pars,
    fixed.pars = fixed.pars
  )

  est <- rugarch::ugarchfit(spec = spec, data = rt_core)
  sig <- est@fit$sigma
  etat <- est@fit$residuals
  coef <- est@fit$coef
  if (model == "gjrGARCH") {
    coef <- coef[names(coef) != "delta"]
  }
  serrors <- est@fit$se.coef
  names(serrors) <- names(coef)
  vcov_mat <- est@fit$cvar
  llh <- est@fit$LLH
  inf_crit <- rugarch::infocriteria(est)[c(1, 2)]
  names(inf_crit) <- c("aic", "bic")

  if (inherits(rt, "zoo")) {
    tp <- as.Date(time(rt))
    sig <- zoo::zoo(sig, order.by = tp)
    etat <- zoo::zoo(etat, order.by = tp)
  } else if (inherits(rt, "ts")) {
    start_rt <- stats::start(rt)
    frequ_rt <- stats::frequency(rt)
    sig <- stats::ts(sig, start = start_rt, frequency = frequ_rt)
    etat <- stats::ts(etat, start = start_rt, frequency = frequ_rt)
  }

  out <- fEGarch_fit_rugarch_wrapper(
    rt = rt,
    sigt = sig,
    etat = etat,
    pars = coef,
    se = serrors,
    vcov_mat = vcov_mat,
    llhood = llh,
    orders = orders,
    cond_dist = cond_dist,
    inf_criteria = inf_crit,
    rugarch_model = est,
    long_memo = FALSE,
    model = model
  )
  out
}

#'Wrapper Functions for Selected \code{rugarch} GARCH Models
#'
#'Easy to use functions for fitting selected GARCH-type models from
#'the widely-known \code{rugarch} package by Ghalanos (2024). These
#'functions are meant as an easy to use way to compare the main results from
#'\code{rugarch} to the newly established models in the \code{fEGarch} package
#'and are by no means considered to be a replacement of \code{rugarch}.
#'
#'@param rt the input time series to fit the model to ordered from
#'past to present; can also be a \code{"zoo"} class object or a
#'\code{"ts"} class object.
#'@param orders the ARCH and GARCH orders of the model as a two-element
#'numeric vector.
#'@param cond_dist a single-element character vector with the conditional
#'distribution to consider.
#'@param include_mean a logical that indicates whether or not to include
#'the mean in the maximum-likelihood estimation. If \code{FALSE}, the mean
#'is fixed to zero.
#'@param start_pars a named list with starting parameters.
#'
#'@details
#'For details, please see the documentation of the \code{rugarch} package
#'(Ghalanos, 2024).
#'
#'@export
#'
#'@return
#'A list with the following named elements is returned.
#'\describe{
#'\item{\code{pars}:}{a named numeric vector with the parameter estimates.}
#'\item{\code{se}:}{a named numeric vector with the obtained standard errors in accordance with the parameter estimates.}
#'\item{\code{vcov_mat}:}{the variance-covariance matrix of the parameter estimates with named columns and rows.}
#'\item{\code{rt}:}{the input object \code{rt}; if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is kept.}
#'\item{\code{sigt}:}{the estimated conditional standard deviations; if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{sigt}.}
#'\item{\code{etat}:}{the obtained residuals; if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{etat}.}
#'\item{\code{orders}:}{a two-element numeric vector stating the considered model orders.}
#'\item{\code{cond_dist}:}{a character value stating the conditional distribution considered in the model fitting.}
#'\item{\code{llhood}:}{the log-likelihood value obtained at the optimal parameter combination.}
#'\item{\code{inf_criteria}:}{a named two-element numeric vector with the corresponding AIC (first element) and BIC (second element) of the fitted model.}
#'\item{\code{rugarch_model}:}{the estimation object returned by \code{ugarchfit()} of the \code{rugarch} package (Ghalanos, 2024).}
#'}
#'
#'@references
#'\itemize{
#'\item{Ghalanos, A. (2024). \code{rugarch}: Univariate GARCH models. R package
#'version 1.5-3. DOI: 10.32614/CRAN.package.rugarch.}
#'}
#'
#'@name rugarch_wrappers
#'
#'@examples
#'est <- gjrgarch_ru(SP500)
#'est$pars
#'est$se
#'plot(est$sigt)
#'
aparch_ru <- function(rt, orders = c(1, 1),
            cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
            include_mean = TRUE,
            start_pars = NULL) {

  rugarch_modelfit(
    rt = rt,
    orders = orders,
    cond_dist = cond_dist,
    model = "apARCH",
    include_mean = include_mean, start_pars = start_pars
  )

}

#'@export
#'@rdname rugarch_wrappers

gjrgarch_ru <- function(rt, orders = c(1, 1),
            cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
            include_mean = TRUE,
            start_pars = NULL) {

  rugarch_modelfit(
    rt = rt,
    orders = orders,
    cond_dist = cond_dist,
    model = "gjrGARCH",
    include_mean = include_mean, start_pars = start_pars
  )

}

#'@export
#'@rdname rugarch_wrappers

egarch_ru <- function(rt, orders = c(1, 1),
            cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
            include_mean = TRUE,
            start_pars = NULL) {

  rugarch_modelfit(
    rt = rt,
    orders = orders,
    cond_dist = cond_dist,
    model = "eGARCH",
    include_mean = include_mean, start_pars = start_pars
  )

}

#'@export
#'@rdname show-methods
#'@aliases show,fEgarch_fit_rugarch_wrapper-method
setMethod(
  "show",
  "fEGarch_fit_rugarch_wrapper",
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
    "* Fitted GARCH-Type Model (rugarch) *\n",
    "*************************************\n",
    " \n",
    paste0("GARCH-type model (rugarch): ", x@model,"\n"),
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

