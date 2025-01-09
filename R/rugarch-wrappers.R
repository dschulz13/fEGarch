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
  } else {
    fixed.pars <- list()
  }

  spec <- rugarch::ugarchspec(
    variance.model = list(model = "apARCH", garchOrder = orders),
    mean.model = list(armaOrder = c(0, 0), include.mean = include_mean),
    distribution.model = cond_dist, start.pars = start.pars,
    fixed.pars = fixed.pars
  )

  est <- rugarch::ugarchfit(spec = spec, data = rt_core)
  sig <- est@fit$sigma
  etat <- est@fit$residuals
  coef <- est@fit$coef
  if (model == "gjrGARCH") {
    coef <- head(coef, -1)
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
  }

  out <- list(
    obs = rt,
    sigt = sig,
    etat = etat,
    pars = coef,
    se = serrors,
    vcov_mat = vcov_mat,
    llhood = llh,
    orders = orders,
    cond_dist = cond_dist,
    inf_criteria = inf_crit
  )
  out
}

#'@export

aparch <- function(rt, orders = c(1, 1),
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

gjrgarch <- function(rt, orders = c(1, 1),
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
