#=============================================================#

# Setup for negative log-likelihood

nllhood_setup <- function(theta, y, pdf_fun, extra_par, skew_par, mean_par, sdev_par) {

  idx <- 1
  logic_extra <- is.logical(extra_par)
  if (logic_extra && extra_par) {
    shape <- theta[[idx]]
    idx <- idx + 1
  } else if (logic_extra && !extra_par) {
    shape <- 0
  } else if (!logic_extra && is.numeric(extra_par)) {
    shape <- extra_par
  }

  if (skew_par) {
    skew <- theta[[idx]]
    idx <- idx + 1
  } else {
    skew <- 1
  }

  if (is.null(mean_par)) {
    mean <- theta[[idx]]
    idx <- idx + 1
  } else {
    mean <- mean_par
  }

  if (is.null(sdev_par)) {
    sdev <- theta[[idx]]
  } else {
    sdev <- sdev_par
  }

  # Return negative log-likelihood
  -sum(log(pdf_fun(x = y, shape = shape, skew = skew, mean = mean, sdev = sdev)))

}

#==============================================#

# Estimation functions

### General estimation function

#'@export
distr_est <- function(x, dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                      fix_mean = NULL, fix_sdev = NULL, Prange = c(1, 5)) {

  dist <- match.arg(dist)

  dfun_s <- lookup_table[["fun1"]][[dist]]
  dfun <- function(x, shape, skew, mean, sdev) {
    dfun_s(x = (x - mean) / sdev, shape = shape, skew = skew) / sdev
  }

  extra_par <- switch(
    dist,
    "norm" = FALSE,
    "std" = TRUE,
    "ged" = TRUE,
    "ald" = TRUE,
    "snorm" = FALSE,
    "sstd" = TRUE,
    "sged" = TRUE,
    "sald" = TRUE
  )

  skew_par <- switch(
    dist,
    "norm" = FALSE,
    "std" = FALSE,
    "ged" = FALSE,
    "ald" = FALSE,
    "snorm" = TRUE,
    "sstd" = TRUE,
    "sged" = TRUE,
    "sald" = TRUE
  )

  # Fix the negative log-likelihood
  nllhood_s <- function(theta, y, extra_par) {
    nllhood_setup(
      theta = theta,
      y = y,
      pdf_fun = dfun,
      extra_par = extra_par,
      skew_par = skew_par,
      mean_par = fix_mean,
      sdev_par = fix_sdev
    )
  }

  if (is.logical(extra_par) && extra_par) {
    start_extra_par <- switch(
      dist,
      "std" = 10,
      "ged" = 2,
      "sstd" = 10,
      "sged" = 2
    )
    LB_extra_par <- switch(
      dist,
      "std" = 2,
      "ged" = 1,
      "sstd" = 2,
      "sged" = 1
    )
    UB_extra_par <- switch(
      dist,
      "std" = Inf,
      "ged" = 5,
      "sstd" = Inf,
      "sged" = 5
    )
    name_extra_par <- switch(
      dist,
      "norm" = character(0),
      "std" = "df",
      "ged" = "shape",
      "ald" = "P",
      "snorm" = character(0),
      "sstd" = "df",
      "sged" = "shape",
      "sald" = "P"
    )
  } else {
    start_extra_par <- LB_extra_par <- UB_extra_par <- numeric(0)
    name_extra_par <- character(0)
  }

  if (skew_par) {
    start_skew <- 1
    LB_skew <- 0.3
    UB_skew <- 3
    name_skew <- "skew"
  } else {
    start_skew <- LB_skew <- UB_skew <- numeric(0)
    name_skew <- character(0)
  }

  if (is.null(fix_mean)) {
    start_mean <- mean(x)
    LB_mean <- min(x)
    UB_mean <- max(x)
    name_mean <- "mean"
  } else {
    start_mean <- LB_mean <- UB_mean <- numeric(0)
    name_mean <- character(0)
  }

  if (is.null(fix_sdev)) {
    SD_x <- sd(x)
    start_sdev <- SD_x
    LB_sdev <- SD_x / 4
    UB_sdev <- SD_x * 4
    name_sdev <- "sdev"
  } else {
    start_sdev <- LB_sdev <- UB_sdev <- numeric(0)
    name_sdev <- character(0)
  }

  start <- c(
    start_extra_par,
    start_skew,
    start_mean,
    start_sdev
  )
  LB <- c(
    LB_extra_par,
    LB_skew,
    LB_mean,
    LB_sdev
  )
  UB <- c(
    UB_extra_par,
    UB_skew,
    UB_mean,
    UB_sdev
  )

  check_lstart <- length(start) > 0
  check_ald <- dist %in% c("ald", "sald")

  if (check_ald) {

    P_all <- seq(Prange[[1]], Prange[[2]], 1)
    i <- 1
    lP <- length(P_all)
    ests <- vector(mode = "list", length = lP)
    nllh <- rep(NA, lP)

    for (P0 in P_all) {

      nllhood <- function(theta, y) {
        nllhood_s(
          theta = theta,
          y = y,
          extra_par = P0
        )
      }

      ests[[i]] <- if (check_lstart) {
        suppressWarnings(Rsolnp::solnp(
          pars = start, fun = nllhood, LB = LB, UB = UB, y = x,
          control = list(trace = FALSE)
        ))
      } else {
        list(
          pars = numeric(0),
          values = nllhood(theta = numeric(0), y = x)
        )
      }

      nllh[[i]] <- utils::tail(ests[[i]]$values, 1)

      i <- i + 1
    }

    ind <- which.min(nllh)

    est <- ests[[ind]]
    P <- P_all[[ind]]

  } else {

    nllhood <- function(theta, y) {
      nllhood_s(
        theta = theta,
        y = y,
        extra_par = extra_par
      )
    }

    est <- if (check_lstart) {

      suppressWarnings(Rsolnp::solnp(
        pars = start, fun = nllhood, LB = LB, UB = UB, y = x,
        control = list(trace = FALSE)
      ))

    } else {

      list(
        pars = numeric(0),
        values = nllhood(theta = numeric(0), y = x)
      )

    }

  }

  par <- est$pars
  par_names <- c(
    name_extra_par,
    name_skew,
    name_mean,
    name_sdev
  )

  if (check_lstart) {
    if (check_ald) {
      nllhood <- function(theta, y) {
       nllhood_s(
         theta = theta,
         y = y,
         extra_par = P
       )
      }
    }
    hess <- numDeriv::hessian(nllhood, x = par, y = x)
    vcov <- solve(hess)
    serror <- sqrt(diag(vcov))
    if (check_ald) {
      par <- c(P, par)
      serror <- c(NA, serror)
      vcov <- rbind(NA, cbind(NA, vcov))
    }
    names(par) <- par_names
    names(serror) <- par_names
    rownames(vcov) <- par_names
    colnames(vcov) <- par_names
  } else {
    if (check_ald) {
      par <- P
      serror <- NA
      vcov <- matrix(NA, ncol = 1, nrow = 1)
      names(par) <- "P"
      names(serror) <- "P"
      rownames(vcov) <- "P"
      colnames(vcov) <- "P"
    } else {
      vcov <- serror <- numeric(0)
    }

  }

  k <- length(par)
  llhood <- -utils::tail(est$values, 1)
  n <- length(x)
  m2llhood <- -2 * llhood

  aic <- (m2llhood + 2 * k) / n
  bic <- (m2llhood + k * log(n)) / n

  if (!is.null(fix_mean) || !is.null(fix_sdev)) {
    fm <- if (!is.null(fix_mean)) {
      name_m <- "mean"
      fix_mean
    } else {
      name_m <- character(0)
      numeric(0)
    }

    sm <- if (!is.null(fix_sdev)) {
      name_s <- "sdev"
      fix_sdev
    } else {
      name_s <- character(0)
      numeric(0)
    }

    fixed <- c(fm, sm)
    names(fixed) <- c(name_m, name_s)

  } else {
    fixed <- numeric(0)
  }



  list(
    par = par, serror = serror, aic = aic, bic = bic, llhood = llhood, vcov = vcov,
    fixed = fixed
  )

}

# Special normal distribution estimation wrapper

#'@export
norm_est <- function(x, fix_mean = NULL, fix_sdev = NULL) {

  distr_est(x = x, dist = "norm", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = c(1, 5))

}

# Special t-distribution estimation wrapper

#'@export
std_est <- function(x, fix_mean = NULL, fix_sdev = NULL) {

  distr_est(x = x, dist = "std", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = c(1, 5))

}

# Special GED estimation wrapper

#'@export
ged_est <- function(x, fix_mean = NULL, fix_sdev = NULL) {

  distr_est(x = x, dist = "ged", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = c(1, 5))

}

# Special ALD estimation wrapper

#'@export
ald_est <- function(x, fix_mean = NULL, fix_sdev = NULL, Prange = c(1, 5)) {

  distr_est(x = x, dist = "ald", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = Prange)

}

# Special skewed normal distribution estimation wrapper

#'@export
snorm_est <- function(x, fix_mean = NULL, fix_sdev = NULL) {

  distr_est(x = x, dist = "snorm", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = c(1, 5))

}

# Special skewed t-distribution estimation wrapper

#'@export
sstd_est <- function(x, fix_mean = NULL, fix_sdev = NULL) {

  distr_est(x = x, dist = "sstd", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = c(1, 5))

}

# Special skewed GED estimation wrapper

#'@export
sged_est <- function(x, fix_mean = NULL, fix_sdev = NULL) {

  distr_est(x = x, dist = "sged", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = c(1, 5))

}

# Special skewed ALD estimation wrapper

#'@export
sald_est <- function(x, fix_mean = NULL, fix_sdev = NULL, Prange = c(1, 5)) {

  distr_est(x = x, dist = "sald", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = Prange)

}
