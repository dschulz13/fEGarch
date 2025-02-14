setGeneric("measure_risk", function(object, measure = c("VaR", "ES"), level = c(0.975, 0.99), ...) {standardGeneric("measure_risk")})

#'VaR and ES Computation Following Fitted Models or Forecasts
#'
#'Provides easy access to value-at-risk (VaR) and expected shortfall (ES)
#'computation for available models in this package. VaR and ES can either
#'be computed based on (a) fitted conditional means and conditional
#'standard deviations for a training period, or following (b) point forecasts
#'(either multistep or rolling) of the conditional means and conditional
#'standard deviations.
#'
#'@param object either an object of class \code{"fEGarch_fit"} returned by the
#'fitting / estimation functions of this package like returned by for example
#'\code{\link{fEGarch}} among others, or an object of class \code{"fEGarch_forecast"}
#'as returned by \code{\link{predict,fEGarch_fit-method}} or
#'\code{\link{predict_roll,fEGarch_fit-method}}.
#'@param a character vector with element \code{"VaR"}, \code{"ES"} or both;
#'indicates, what risk measures should be computed; by default, both VaR and ES
#'are computed.
#'@param level a numeric vector of arbitrary length indicating the confidence
#'levels to compute the VaR and the ES at; by default, the levels
#'\code{0.975} and \code{0.99}, i.e. 97.5 percent and 99 percent, are considered.
#'@param ... currently without use.
#'
#'@details
#'Given a fitted model with fitted conditional means and conditional standard deviations
#'or given point forecasts of such series based on a fitted model, the risk measures
#'VaR and ES can be computed (at arbitrary confidence levels) following the conditional
#'loss distribution defined through the estimated / forecasted conditional mean value,
#'the estimated / forecasted conditional standard deviation value, and the assumed
#'conditional distribution (including potential estimates of distribution parameters).
#'
#'Let \eqn{\hat{\mu}_t} be the estimated / forecasted conditional mean and \eqn{\hat \sigma_t} be the
#'estimated / forecasted conditional standard deviation at some time point \eqn{t}. Furthermore,
#'define \eqn{\text{VaR}_{\eta,\alpha}} and \eqn{\text{ES}_{\eta,\alpha}} be
#'the time-invariant VaR and ES, respectively, of some identically but independently
#'distributed random variables \eqn{\eta_t} with mean zero and variance one. Given that
#'the relationship \eqn{r_t = \mu_t + \sigma_t\eta_t}, where \eqn{\mu_t} and \eqn{\sigma_t}
#'are the true conditional mean and conditional standard deviation at time \eqn{t}, is assumed
#'for some return series \eqn{\{r_t\}}, the estimated / forecasted conditional VaR and ES of \eqn{r_t} are simply
#'\deqn{\widehat{\text{VaR}}_{r,\alpha}(t) = \hat{\mu}_t + \hat{\sigma}_t \text{VaR}_{\eta,\alpha}  \hspace{3mm} \text{ and } \hspace{3mm} \widehat{\text{ES}}_{r,\alpha}(t) = \hat{\mu}_t + \hat{\sigma}_t \text{ES}_{\eta,\alpha}.}
#'This definition holds, when losses and therefore also
#'\eqn{\text{VaR}_{\eta,\alpha}(t)} and \eqn{\text{ES}_{\eta,\alpha}(t)} (for common \eqn{\alpha} such as
#'\eqn{\alpha = 0.975} or \eqn{\alpha = 0.99}) are considered
#'to be negative in sign.
#'
#'Define
#'\deqn{\text{VaR}_{\eta,\alpha} = f_{\eta}^{-1}(1-\alpha) \hspace{3mm} \text{ and } \hspace{3mm} \text{ES}_{\eta,\alpha} = (1-\alpha)^{-1}\int_{\alpha}^{\infty} \text{VaR}_{\eta, x} dx,}
#'which also need to be estimated for some distributions, if a distribution parameter needed to be estimated. \eqn{f} in the previous formula
#'is the cumulative distribution function of the random variables \eqn{\eta_t}. Therefore,
#'\eqn{f^{-1}_{\eta}(1-\alpha)} returns the quantile of the innovation distribution at level
#'\eqn{1-\alpha}.
#'
#'@rdname measure_risk
#'
#'@export
#'
#'@return
#'Returns a two-element list with elements \code{measures} and \code{observations}.
#'\code{observations} is the observed series at the time points, for which the
#'risk measures are calculated. \code{measures} is again a list with elements
#'\code{VaR} and \code{ES}, distinguishing between computed VaR and ES values.
#'These elements again are list with named elements representing the various
#'computed series.
#'
#'@examples
#'
#'# In-sample
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- fEGarch(egarch_spec(), rt)
#'risk <- measure_risk(model, measure = c("VaR", "ES"), level = c(0.95, 0.975, 0.99))
#'risk
#'
#'# Out-of-sample rolling point forecasts
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model2 <- fEGarch(egarch_spec(), rt, n_test = 250)
#'fcast <- predict_roll(model2)
#'risk2 <- measure_risk(fcast, measure = c("VaR", "ES"), level = c(0.95, 0.975, 0.99))
#'risk2
#'
#'
#'
#'
setMethod("measure_risk", "fEGarch_fit",
  function(object, measure = c("VaR", "ES"), level = c(0.975, 0.99), ...) {
    cond_d <- object@cond_dist

    par_name <- switch(
      cond_d,
      "norm" = "df",   # irrelevant
      "snorm" = "df",  # irrelevant
      "std" = "df",
      "sstd" = "df",
      "ged" = "shape",
      "sged" = "shape",
      "ald" = "P",
      "sald" = "P"
    )
    shape <- if (!(cond_d %in% c("norm", "snorm"))) {
      object@pars[[par_name]]
    } else {
      0
    }
    skew <- if (!(cond_d %in% c("norm", "std", "ged", "ald"))) {
      object@pars[["skew"]]
    } else {
      0
    }


    args <- list(
      level = level,
      dist = cond_d,
      skew = skew
    )
    args[[par_name]] <- shape

    if ("VaR" %in% measure) {
      q_VaR <- do.call(what = VaR_calc, args = args)
      names_q <- paste0("VaR", level)
      lVaR <- list()
      for (i in seq_along(q_VaR)) {
        lVaR[[names_q[[i]]]] <- zoo::zoo(zoo::coredata(object@cmeans + q_VaR[[i]] * object@sigt), order.by = time(object@sigt))
      }
    } else {
      lVaR <- list()
    }

    if ("ES" %in% measure) {
      q_ES <- do.call(what = ES_calc, args = args)
      names_q <- paste0("ES", level)
      lES <- list()
      for (i in seq_along(q_ES)) {
        lES[[names_q[[i]]]] <- zoo::zoo(zoo::coredata(object@cmeans + q_ES[[i]] * object@sigt), order.by = time(object@sigt))
      }
    } else {
      lES <- list()
    }

    fEGarch_risk(
      measures = list(VaR = lVaR, ES = lES),
      observations = object@rt
    )

  }
)

#'@export
#'@rdname measure_risk
setMethod("measure_risk", "fEGarch_forecast",
  function(object, measure = c("VaR", "ES"), level = c(0.975, 0.99), ...) {

    object2 <- object@model
    object2@sigt <- object@sigt
    object2@cmeans <- object@cmeans

    func <- selectMethod("measure_risk", "fEGarch_fit")
    out <- func(object = object2, measure = measure, level = level, ...)
    out@observations <- object2@test_obs
    out
  }
)


#'@export
varcast_fEGarch <- function(
    x,
    level_VaR = 0.99,
    level_ES = 0.975,
    par_model_spec = list(
      spec = egarch_spec(), include_mean = TRUE, drange = c(0, 0.5 - 1e-06), start_pars = NULL, LB = NULL, UB = NULL, control = list(), parallel = TRUE, ncores = max(1, future::availableCores() - 1)
    ),
    nonpar_lm_model_spec = list(p = 3, pmax = 1, qmax = 1, pmin = 0, qmin = 0, mu = 1, InfR = "Nai", bStart = 0.15, bb = 1, cb = 0.05),
    nonpar_sm_model_spec = list(p = 3, mu = 1, bStart = 0.15, alg = "B"),
    n.out = 250,
    smooth = FALSE
) {

  args_p <- list(
      spec = egarch_spec(), include_mean = TRUE, drange = c(0, 0.5 - 1e-06),
      start_pars = NULL, LB = NULL, UB = NULL, control = list(),
      parallel = TRUE, ncores = max(1, future::availableCores() - 1)
  )
  elems <- names(par_model_spec)
  for (nam in elems) {
    args_p[[nam]] <- par_model_spec[[nam]]
  }


  rt <- diff(log(x))

  a.v <- level_VaR
  a.e <- level_ES

  if (!smooth) {

    est <- fEGarch(spec = args_p$spec, rt = rt, drange = args_p$drange,
                   meanspec = mean_spec(include_mean = args_p$include_mean),
                   start_pars = ars_p$start_pars,
                   LB = args_p$LB,
                   UB = args_p$UB, n_test = n.out,
                   control = args_p$control,
                   parallel = args_p$parallel,
                   ncores = args_p$ncores)

    fcast <- predict_roll(est)
    risk <- measure_risk(fcast, level = c(a.v, a.e))
    name.v <- paste0("VaR", a.v)
    name.VaR.e <- paste0("VaR", a.e)
    name.ES.e <- paste0("ES", a.e)
    VaR.v <- -risk@measures$VaR[[name.v]]
    VaR.e <- -risk@measures$VaR[[name.VaR.e]]
    ES.e <- -risk@measures$ES[[name.ES.e]]
    mu <- if (args_p$include_mean) {
      est@pars[["mu"]]
    } else {
      0
    }
    est_np <- NULL
    sig.in <- est@sigt
    sig.fc <- fcast@sigt
    scale_fc <- NULL
    scale <- NULL
    est_np <- NULL
  } else {

    rt_tr <- utils::head(rt, -n.out)
    rt_te <- utils::tail(rt, n.out)

    mu <- mean(rt_tr)
    rt_dm <- rt_tr - mu
    yt <- zoo::coredata(log((rt_dm)^2))


    lm <- long_memo(args_p$spec)
    if (!lm) {
      func <- smoots::msmooth
      input_np <- nonpar_sm_model_spec
      args_np <- list(p = 3, mu = 1, bStart = 0.15, alg = "B")
    } else {
      func <- esemifar::tsmoothlm
      input_np <- nonpar_lm_model_spec
      args_np <- list(p = 3, pmax = 1, qmax = 1, pmin = 0, qmin = 0, mu = 1, InfR = "Nai", bStart = 0.15, bb = 1, cb = 0.05)
    }
    names_np <- names(input_np)
    for (j in names_np) {
      args_np[[j]] <- input_np[[j]]
    }

    args_np[["y"]] <- yt

    est_np <- do.call(what = func, args = args_np)
    res <- est_np$res
    C <- -log(mean(exp(res)))
    scale_f <- exp((est_np$ye - C) / 2)

    scale_fcast <- rep(utils::tail(scale_f, 1), n.out)

    eps <- c(rt_dm / scale_f, (rt_te - mu) / scale_fcast)

    est <- fEGarch(
      spec = args_p$spec,
      rt = eps,
      drange = args_p$drange,
      meanspec = mean_spec(include_mean = FALSE),
      n_test = n.out,
      start_pars = args_p$start_pars,
      LB = args_p$LB,
      UB = args_p$UB,
      control = args_p$control,
      parallel = args_p$parallel,
      ncores = args_p$ncores
    )

    fcast <- predict_roll(est)
    risk <- measure_risk(fcast, level = c(a.v, a.e))
    name.v <- paste0("VaR", a.v)
    name.VaR.e <- paste0("VaR", a.e)
    name.ES.e <- paste0("ES", a.e)
    VaR.v <- -risk@measures$VaR[[name.v]] * scale_fcast - mu
    VaR.e <- -risk@measures$VaR[[name.VaR.e]] * scale_fcast - mu
    ES.e <- -risk@measures$ES[[name.ES.e]] * scale_fcast - mu

    sig.in <- est@sigt * scale_f
    sig.fc <- fcast@sigt * scale_fcast
    scale <- scale_f
    scale_fc <- scale_fcast

  }

  list(
    mean = mu,
    ret_in = rt_tr,
    ret_out = rt_te,
    sig_in = sig.in,
    sig_fc = sig.fc,
    scale = scale,
    scale_fc = scale_fc,
    VaR_e = VaR.e,
    VaR_v = VaR.v,
    ES = ES.e,
    level_VaR = a.v,
    level_ES = a.e,
    model_par = est,
    model_nonpar = est_np
  )

}

#'@export
#'@importFrom ufRisk covtest
ufRisk::covtest

#'@export
#'@importFrom ufRisk trafftest
ufRisk::trafftest

#'@export
#'@importFrom ufRisk lossfunc
ufRisk::lossfunc

