
#### MAFILog-GARCH (complete)

sigt_fiegarch <- function(theta, x, dfun1, dfun2, p, q, include_mean,
                             extra_par, skew_par) {

  imean <- !include_mean
  mu <- c(0, theta[[1]])[[include_mean + 1]]
  step0 <- 2 - imean
  Elnsig2 <- theta[[step0]]   # The expectation of ln(sig_t^2)
  step1 <- step0 + p
  step2 <- step1 + q - 1

  phi <- theta[(step0 + 1):step1]

  if (q > 1) {
    psi <- c(1, theta[(step1 + 1):step2])
    psi_2 <- psi[-1]
  } else {
    psi <- 1
    psi_2 <- numeric(0)
  }
  kappa <- theta[[(step2 + 1)]]
  gamma <- theta[[step2 + 2]]
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

  help_fun1 <- function(x, shape, skew, dfun) {
    abs(x) * dfun(x, shape, skew)
  }
  E1 <- integrate(help_fun1, -Inf, Inf, shape = shape, skew = skew,
                  dfun = dfun1,
                  abs.tol = 1e-15, subdivisions = 1000L)[[1]]

  n <- length(x)
  coef_inf <- esemifar::farima_to_ma(ar = phi, ma = psi_2, d = d,
                                     max_i = n - 1)

  alpha <- coef_inf * kappa
  beta <- coef_inf * gamma

  sigt <- c(sigt_fiegarchCpp(
    x = x,
    mu = mu,
    alpha = alpha,
    beta = beta,
    E1 = E1,
    Elnsig2 = Elnsig2
  ))

  list(sigt = sigt, mu = mu, shape = shape, skew = skew)
}

# p, q are now model orders
goal_fun_fiegarch <- function(theta, x, dfun1, dfun2, p, q, include_mean,
                             extra_par, skew_par) {

  res_list <- sigt_fiegarch(theta, x, dfun1, dfun2, p, q,
                include_mean, extra_par, skew_par)

  sigt <- res_list$sigt
  mu <- res_list$mu
  shape <- res_list$shape
  skew <- res_list$skew

  -sum(log(dfun2(x, mu, sigt, shape, skew)))

}

fiegarchFit <- function(rt, order = c(1L, 1L),
          cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
          include_mean = TRUE, start_pars = NULL, LB = NULL, UB = NULL) {

  stopifnot("Both orders must be at least 1." = any(order >= 1))
  cond_dist <- match.arg(cond_dist)

  p <- order[[1]]
  q <- order[[2]]
  force(lookup_table)
  force(include_mean)
  fun1 <- lookup_table[["fun1"]][[cond_dist]]
  fun2 <- lookup_table[["fun2"]][[cond_dist]]

  skew_par <- lookup_table[["skew_par"]][[cond_dist]]
  extra_par <- lookup_table[["extra_par"]][[cond_dist]]

  goal_fun <- function(theta, x) {
    goal_fun_fiegarch(
      theta = theta,
      x = x,
      dfun1 = fun1,
      dfun2 = fun2,
      p = p,
      q = q,
      include_mean = include_mean,
      extra_par = extra_par,
      skew_par = skew_par
    )
  }

  vals <- lookup_table[["vals"]][[cond_dist]]
  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  if (is.null(start_pars)) {

    start_pars <- c(
      list(NULL, 0)[[include_mean + 1]],
      -11,
      rep(0, p),
      rep(0, q - 1),
      rep(0, 2),
      0.1,
      sval,
      list(NULL, 1)[[skew_par + 1]]
    )

  }

  if (is.null(LB)) {

    LB <- c(
      list(NULL, -2)[[include_mean + 1]],
      -20,
      rep(-2, p),
      rep(-2, q - 1),
      rep(-2, 2),
      0,
      lval,
      list(NULL, 1e-15)[[skew_par + 1]]
    )

  }

  if (is.null(UB)) {

    UB <- c(
      list(NULL, 2)[[include_mean + 1]],
      20,
      rep(2, p),
      rep(2, q - 1),
      rep(2, 2),
      0.4999,
      uval,
      list(NULL, Inf)[[skew_par + 1]]
    )

  }

  result <- suppressWarnings(Rsolnp::solnp(
    pars = start_pars,
    fun = goal_fun,
    LB = LB,
    UB = UB,
    control = list(trace = FALSE),
    x = rt
  ))

  pars <- result$pars

  llhood <- -tail(result$values, 1)
  sigt <- sigt_fiegarch(
    theta = pars,
    x = rt,
    dfun1 = fun1,
    dfun2 = fun2,
    p = p,
    q = q,
    include_mean = include_mean,
    extra_par = extra_par,
    skew_par = skew_par
  )$sigt

  mu <- c(0, pars[[1]])[[include_mean + 1]]
  etat <- (rt - mu) / sigt

  if (q > 1) {
    psi_names <- paste0("psi", 2:q)
  } else {
    psi_names <- NULL
  }

  par_names <- c(
    list(NULL, "mu")[[include_mean + 1]],
    "omega_sig",
    paste0("phi", 1:p),
    psi_names,
    "kappa",
    "gamma",
    "d",
    list(NULL, c("shape", "df")[[(cond_dist %in% c("std", "sstd")) + 1]])[[extra_par + 1]],
    list(NULL, "skew")[[skew_par + 1]]
  )

  serrors <- sqrt(diag(solve(result$hessian)))
  names(pars) <- par_names
  names(serrors) <- par_names

  k <- length(pars)
  n <- length(rt)

  m2llhood <- -2 * llhood
  bic <- (m2llhood + k * log(n)) / n
  aic <- (m2llhood + 2 * k) / n

  list_out <- list(
    pars = pars,
    se = serrors,
    rt = rt,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = c("aic" = aic, "bic" = bic),
    cond_dist = cond_dist,
    orders = order
  )

  class(list_out) <- "fiegarch_fit"

  list_out

}

#'@export

print.fiegarch_fit <- function(x) {

  par_names <- names(x$pars)
  se <- unname(x$se)
  pars <- unname(x$pars)

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
    "*           Fitted FIEGARCH         *\n",
    "*************************************\n",
    " \n",
    "Orders: (", x$orders[[1]], ", ", x$orders[[2]], ")\n",
    "Cond. distribution: ", x$cond_dist, "\n",
    " \n",
    "Fitted parameters:\n"
  )

  part2 <- paste0(
    " \n",
    "Information criteria:\n",
    "AIC: ", sprintf("%.4f", unname(x$inf_criteria[[1]])),
    ", BIC: ", sprintf("%.4f", unname(x$inf_criteria[[2]]))
  )

  cat(part1, fill = TRUE)
  print(df)
  cat(part2, fill = TRUE)


}

##################################################

#### FILog-GARCH (complete)

sigt_filoggarch <- function(theta, x, dfun1, dfun2, p, q, include_mean,
                             extra_par, skew_par) {

  imean <- !include_mean
  mu <- c(0, theta[[1]])[[include_mean + 1]]
  step0 <- 2 - imean
  omega_sig <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q
  phi <- theta[(step0 + 1):step1]
  psi <- theta[(step1 + 1):step2]
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

  help_fun1 <- function(x, shape, skew, dfun) {
    log(x^2) * dfun(x, shape, skew)
  }
  E1 <- integrate(help_fun1, -Inf, Inf, shape = shape, skew = skew,
                  dfun = dfun1,
                  abs.tol = 1e-15, subdivisions = 1000L)[[1]]

  n <- length(x)
  coef_inf <- esemifar::farima_to_ma(ar = phi, ma = psi, d = d,
                                     max_i = n)[-1]

  Elnsig2 <- omega_sig


  sigt <- c(sigt_filoggarchCpp(
    x = x,
    mu = mu,
    coef_inf = coef_inf,
    E1 = E1,
    Elnsig2 = Elnsig2
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape)
}

# p, q are now model orders
goal_fun_filoggarch <- function(theta, x, dfun1, dfun2, p, q, include_mean,
                             extra_par, skew_par) {

  res_list <- sigt_filoggarch(theta, x, dfun1, dfun2, p, q,
              include_mean, extra_par, skew_par)
  sigt <- res_list$sigt
  mu <- res_list$mu
  skew <- res_list$skew
  shape <- res_list$shape

  out <- -sum(log(dfun2(x, mu, sigt, shape, skew)))

  out

}

filoggarchFit <- function(rt, order = c(1L, 1L),
          cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
          include_mean = TRUE, start_pars = NULL, LB = NULL, UB = NULL) {

  stopifnot("Both orders must be at least 1." = any(order >= 1))
  cond_dist <- match.arg(cond_dist)

  p <- order[[1]]
  q <- order[[2]]
  force(lookup_table)
  force(include_mean)
  fun1 <- lookup_table[["fun1"]][[cond_dist]]
  fun2 <- lookup_table[["fun2"]][[cond_dist]]

  skew_par <- lookup_table[["skew_par"]][[cond_dist]]
  extra_par <- lookup_table[["extra_par"]][[cond_dist]]

  goal_fun <- function(theta, x) {
    goal_fun_filoggarch(
      theta = theta,
      x = x,
      dfun1 = fun1,
      dfun2 = fun2,
      p = p,
      q = q,
      include_mean = include_mean,
      extra_par = extra_par,
      skew_par = skew_par
    )
  }

  vals <- lookup_table[["vals"]][[cond_dist]]
  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  if (is.null(start_pars)) {

    start_pars <- c(
      list(NULL, 0.0001)[[include_mean + 1]],
      -11,
      rep(0.2, p),
      rep(-0.1, q),
      0.1,
      sval,
      list(NULL, 1)[[skew_par + 1]]
    )

  }

  if (is.null(LB)) {

    LB <- c(
      list(NULL, -0.1)[[include_mean + 1]],
      -15,
      rep(-2, p),
      rep(-2, q),
      0,
      lval,
      list(NULL, 1e-15)[[skew_par + 1]]
    )

  }

  if (is.null(UB)) {

    UB <- c(
      list(NULL, 0.1)[[include_mean + 1]],
      -3,
      rep(2, p),
      rep(2, q),
      0.4999,
      uval,
      list(NULL, Inf)[[skew_par + 1]]
    )

  }

  result <- suppressWarnings(Rsolnp::solnp(
    pars = start_pars,
    fun = goal_fun,
    LB = LB,
    UB = UB,
    control = list(trace = FALSE),
    x = rt
  ))

  pars <- result$pars

  llhood <- -tail(result$values, 1)
  sigt <- sigt_filoggarch(
    theta = pars,
    x = rt,
    dfun1 = fun1,
    dfun2 = fun2,
    p = p,
    q = q,
    include_mean = include_mean,
    extra_par = extra_par,
    skew_par = skew_par
  )$sigt

  mu <- c(0, pars[[1]])[[include_mean + 1]]
  etat <- (rt - mu) / sigt

  par_names <- c(
    list(NULL, "mu")[[include_mean + 1]],
    "omega_sig",
    paste0("phi", 1:p),
    paste0("psi", 1:q),
    "d",
    list(NULL, c("shape", "df")[[(cond_dist %in% c("std", "sstd")) + 1]])[[extra_par + 1]],
    list(NULL, "skew")[[skew_par + 1]]
  )

  serrors <- sqrt(diag(solve(result$hessian)))
  names(pars) <- par_names
  names(serrors) <- par_names

  k <- length(pars)
  n <- length(rt)

  m2llhood <- -2 * llhood
  bic <- (m2llhood + k * log(n)) / n
  aic <- (m2llhood + 2 * k) / n

  list_out <- list(
    pars = pars,
    se = serrors,
    rt = rt,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = c("aic" = aic, "bic" = bic),
    cond_dist = cond_dist,
    orders = order
  )

  class(list_out) <- "filoggarch_fit"

  list_out

}

#'@export

print.filoggarch_fit <- function(x) {

  par_names <- names(x$pars)
  se <- unname(x$se)
  pars <- unname(x$pars)

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
    "*         Fitted FILog-GARCH        *\n",
    "*************************************\n",
    " \n",
    "Orders: (", x$orders[[1]], ", ", x$orders[[2]], ")\n",
    "Cond. distribution: ", x$cond_dist, "\n",
    " \n",
    "Fitted parameters:\n"
  )

  part2 <- paste0(
    " \n",
    "Information criteria:\n",
    "AIC: ", sprintf("%.4f", unname(x$inf_criteria[[1]])),
    ", BIC: ", sprintf("%.4f", unname(x$inf_criteria[[2]]))
  )

  cat(part1, fill = TRUE)
  print(df)
  cat(part2, fill = TRUE)


}

###########################################

#### MAFILog-GARCH (complete)

sigt_mafiloggarch <- function(theta, x, dfun1, dfun2, p, q, include_mean,
                             extra_par, skew_par) {

  imean <- !include_mean
  mu <- c(0, theta[[1]])[[include_mean + 1]]
  step0 <- 2 - imean
  Elnsig2 <- theta[[step0]]   # The expectation of ln(sig_t^2)
  step1 <- step0 + p
  step2 <- step1 + q - 1

  phi <- theta[(step0 + 1):step1]

  if (q > 1) {
    psi <- c(1, theta[(step1 + 1):step2])
    psi_2 <- psi[-1]
  } else {
    psi <- 1
    psi_2 <- numeric(0)
  }
  kappa <- theta[[(step2 + 1)]]
  gamma <- theta[[step2 + 2]]
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

  help_fun1 <- function(x, shape, skew, dfun) {
    (sign(x) * log(abs(x) + 1)) * dfun(x, shape, skew)
  }
  E1 <- integrate(help_fun1, -Inf, Inf, shape = shape, skew = skew,
                  dfun = dfun1,
                  abs.tol = 1e-15, subdivisions = 1000L)[[1]]

  help_fun2 <- function(x, E1 = E1, shape, skew, dfun) {
    abs(sign(x) * log(abs(x) + 1) - E1) * dfun(x, shape, skew)
  }
  E2 <- integrate(help_fun2, -Inf, Inf, E1 = E1, shape = shape,
                  skew = skew,
                  dfun = dfun1,
                  abs.tol = 1e-15, subdivisions = 1000L)[[1]]

  n <- length(x)
  coef_inf <- esemifar::farima_to_ma(ar = phi, ma = psi_2, d = d,
                                     max_i = n - 1)

  alpha <- coef_inf * kappa
  beta <- coef_inf * gamma

  sigt <- c(sigt_mafiloggarchCpp(
    x = x,
    mu = mu,
    alpha = alpha,
    beta = beta,
    E1 = E1,
    E2 = E2,
    Elnsig2 = Elnsig2
  ))

  list(sigt = sigt, mu = mu, shape = shape, skew = skew)
}

# p, q are now model orders
goal_fun <- function(theta, x, dfun1, dfun2, p, q, include_mean,
                             extra_par, skew_par, sigt_fun) {

  res_list <- sigt_fun(theta, x, dfun1, dfun2, p, q,
                include_mean, extra_par, skew_par)

  sigt <- res_list$sigt
  mu <- res_list$mu
  shape <- res_list$shape
  skew <- res_list$skew

  -sum(log(dfun2(x, mu, sigt, shape, skew)))

}

mafiloggarchFit <- function(rt, order = c(1L, 1L),
          cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
          include_mean = TRUE, start_pars = NULL, LB = NULL, UB = NULL) {

  stopifnot("Both orders must be at least 1." = any(order >= 1))
  cond_dist <- match.arg(cond_dist)

  p <- order[[1]]
  q <- order[[2]]
  force(lookup_table)
  force(include_mean)
  fun1 <- lookup_table[["fun1"]][[cond_dist]]
  fun2 <- lookup_table[["fun2"]][[cond_dist]]

  skew_par <- lookup_table[["skew_par"]][[cond_dist]]
  extra_par <- lookup_table[["extra_par"]][[cond_dist]]

  goal_fun <- function(theta, x) {
    goal_fun_mafiloggarch(
      theta = theta,
      x = x,
      dfun1 = fun1,
      dfun2 = fun2,
      p = p,
      q = q,
      include_mean = include_mean,
      extra_par = extra_par,
      skew_par = skew_par
    )
  }

  vals <- lookup_table[["vals"]][[cond_dist]]
  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  if (is.null(start_pars)) {

    start_pars <- c(
      list(NULL, 0)[[include_mean + 1]],
      -11,
      rep(0, p),
      rep(0, q - 1),
      rep(0, 2),
      0.1,
      sval,
      list(NULL, 1)[[skew_par + 1]]
    )

  }

  if (is.null(LB)) {

    LB <- c(
      list(NULL, -2)[[include_mean + 1]],
      -20,
      rep(-2, p),
      rep(-2, q - 1),
      rep(-2, 2),
      0,
      lval,
      list(NULL, 1e-15)[[skew_par + 1]]
    )

  }

  if (is.null(UB)) {

    UB <- c(
      list(NULL, 2)[[include_mean + 1]],
      20,
      rep(2, p),
      rep(2, q - 1),
      rep(2, 2),
      0.4999,
      uval,
      list(NULL, Inf)[[skew_par + 1]]
    )

  }

  result <- suppressWarnings(Rsolnp::solnp(
    pars = start_pars,
    fun = goal_fun,
    LB = LB,
    UB = UB,
    control = list(trace = FALSE),
    x = rt
  ))

  pars <- result$pars

  llhood <- -tail(result$values, 1)
  sigt <- sigt_mafiloggarch(
    theta = pars,
    x = rt,
    dfun1 = fun1,
    dfun2 = fun2,
    p = p,
    q = q,
    include_mean = include_mean,
    extra_par = extra_par,
    skew_par = skew_par
  )$sigt

  mu <- c(0, pars[[1]])[[include_mean + 1]]
  etat <- (rt - mu) / sigt

  if (q > 1) {
    psi_names <- paste0("psi", 2:q)
  } else {
    psi_names <- NULL
  }

  par_names <- c(
    list(NULL, "mu")[[include_mean + 1]],
    "omega_sig",
    paste0("phi", 1:p),
    psi_names,
    "kappa",
    "gamma",
    "d",
    list(NULL, c("shape", "df")[[(cond_dist %in% c("std", "sstd")) + 1]])[[extra_par + 1]],
    list(NULL, "skew")[[skew_par + 1]]
  )

  serrors <- sqrt(diag(solve(result$hessian)))
  names(pars) <- par_names
  names(serrors) <- par_names

  k <- length(pars)
  n <- length(rt)

  m2llhood <- -2 * llhood
  bic <- (m2llhood + k * log(n)) / n
  aic <- (m2llhood + 2 * k) / n

  list_out <- list(
    pars = pars,
    se = serrors,
    rt = rt,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = c("aic" = aic, "bic" = bic),
    cond_dist = cond_dist,
    orders = order
  )

  class(list_out) <- "mafiloggarch_fit"

  list_out

}

#'@export

print.mafiloggarch_fit <- function(x) {

  par_names <- names(x$pars)
  se <- unname(x$se)
  pars <- unname(x$pars)

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
    "*        Fitted MAFILog-GARCH       *\n",
    "*************************************\n",
    " \n",
    "Orders: (", x$orders[[1]], ", ", x$orders[[2]], ")\n",
    "Cond. distribution: ", x$cond_dist, "\n",
    " \n",
    "Fitted parameters:\n"
  )

  part2 <- paste0(
    " \n",
    "Information criteria:\n",
    "AIC: ", sprintf("%.4f", unname(x$inf_criteria[[1]])),
    ", BIC: ", sprintf("%.4f", unname(x$inf_criteria[[2]]))
  )

  cat(part1, fill = TRUE)
  print(df)
  cat(part2, fill = TRUE)


}

#################################

#### MEGARCH (complete)

sigt_megarch <- function(theta, x, dfun1, dfun2, p, q, include_mean,
                             extra_par, skew_par) {

  imean <- !include_mean
  mu <- c(0, theta[[1]])[[include_mean + 1]]
  step0 <- 2 - imean
  omega_sig <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q - 1

  phi <- theta[(step0 + 1):step1]

  if (q > 1) {
    psi <- c(1, theta[(step1 + 1):step2])
  } else {
    psi <- 1
  }

  kappa <- theta[[(step2 + 1)]]
  gamma <- theta[[step2 + 2]]
  if (extra_par) {
    shape <- theta[[step2 + 3]]
  } else {
    shape <- 0
  }
  if (skew_par) {
    skew <- theta[[step2 + 3 + extra_par]]
  } else {
    skew <- 0
  }

  help_fun1 <- function(x, shape, skew, dfun) {
    (sign(x) * log(abs(x) + 1)) * dfun(x, shape, skew)
  }
  E1 <- integrate(help_fun1, -Inf, Inf, shape = shape, skew = skew,
                  dfun = dfun1,
                  abs.tol = 1e-15, subdivisions = 1000L)[[1]]

  help_fun2 <- function(x, shape, skew, dfun) {
    abs(x) * dfun(x, shape, skew)
  }
  E2 <- integrate(help_fun2, -Inf, Inf, shape = shape, skew = skew,
                  dfun = dfun1,
                  abs.tol = 1e-15, subdivisions = 1000L)[[1]]

  alpha <- kappa * psi
  beta <- gamma * psi

  omega <- omega_sig * (1 - sum(phi))
  Elnsig2 <- omega_sig

  sigt <- c(sigt_megarchCpp(
    x = x,
    mu = mu,
    omega = omega,
    phi = phi,
    alpha = alpha,
    beta = beta,
    E1 = E1,
    E2 = E2,
    Elnsig2 = Elnsig2
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape)
}

# p, q are now model orders
goal_fun_megarch <- function(theta, x, dfun1, dfun2, p, q, include_mean,
                             extra_par, skew_par) {

  res_list <- sigt_megarch(theta, x, dfun1, dfun2, p, q,
               include_mean, extra_par, skew_par)

  sigt <- res_list$sigt
  mu <- res_list$mu
  shape <- res_list$shape
  skew <- res_list$skew

  -sum(log(dfun2(x, mu, sigt, shape, skew)))
}

megarchFit <- function(rt, order = c(1L, 1L),
                       cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
                       include_mean = TRUE, start_pars = NULL, LB = NULL, UB = NULL) {

  stopifnot("Both orders must be at least 1." = any(order >= 1))
  cond_dist <- match.arg(cond_dist)

  p <- order[[1]]
  q <- order[[2]]
  force(lookup_table)
  force(include_mean)
  fun1 <- lookup_table[["fun1"]][[cond_dist]]
  fun2 <- lookup_table[["fun2"]][[cond_dist]]

  skew_par <- lookup_table[["skew_par"]][[cond_dist]]
  extra_par <- lookup_table[["extra_par"]][[cond_dist]]

  goal_fun <- function(theta, x) {
    goal_fun_megarch(
      theta = theta,
      x = x,
      dfun1 = fun1,
      dfun2 = fun2,
      p = p,
      q = q,
      include_mean = include_mean,
      extra_par = extra_par,
      skew_par = skew_par
    )
  }

  vals <- lookup_table[["vals"]][[cond_dist]]
  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  if (is.null(start_pars)) {

    start_pars <- c(
      list(NULL, 0)[[include_mean + 1]],
      0,
      rep(0, p),
      rep(0, q - 1),
      rep(0, 2),
      sval,
      list(NULL, 1)[[skew_par + 1]]
    )

  }

  if (is.null(LB)) {

    LB <- c(
      list(NULL, -Inf)[[include_mean + 1]],
      -Inf,
      rep(-2, p),
      rep(-2, q - 1),
      rep(-2, 2),
      lval,
      list(NULL, 1e-15)[[skew_par + 1]]
    )

  }

  if (is.null(UB)) {

    UB <- c(
      list(NULL, Inf)[[include_mean + 1]],
      Inf,
      rep(2, p),
      rep(2, q - 1),
      rep(2, 2),
      uval,
      list(NULL, Inf)[[skew_par + 1]]
    )

  }

  result <- suppressWarnings(Rsolnp::solnp(
    pars = start_pars,
    fun = goal_fun,
    LB = LB,
    UB = UB,
    control = list(trace = FALSE),
    x = rt
  ))

  pars <- result$pars

  llhood <- -tail(result$values, 1)
  sigt <- sigt_megarch(
    theta = pars,
    x = rt,
    dfun1 = fun1,
    dfun2 = fun2,
    p = p,
    q = q,
    include_mean = include_mean,
    extra_par = extra_par,
    skew_par = skew_par
  )$sigt

  mu <- c(0, pars[[1]])[[include_mean + 1]]
  etat <- (rt - mu) / sigt

  if (q > 1) {
    psi_names <- paste0("psi", 2:q)
  } else {
    psi_names <- NULL
  }

  par_names <- c(
    list(NULL, "mu")[[include_mean + 1]],
    "omega_sig",
    paste0("phi", 1:p),
    psi_names,
    "kappa",
    "gamma",
    list(NULL, c("shape", "df")[[(cond_dist %in% c("std", "sstd")) + 1]])[[extra_par + 1]],
    list(NULL, "skew")[[skew_par + 1]]
  )

  serrors <- sqrt(diag(solve(result$hessian)))
  names(pars) <- par_names
  names(serrors) <- par_names

  k <- length(pars)
  n <- length(rt)

  m2llhood <- -2 * llhood
  bic <- (m2llhood + k * log(n)) / n
  aic <- (m2llhood + 2 * k) / n

  list_out <- list(
    pars = pars,
    se = serrors,
    rt = rt,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = c("aic" = aic, "bic" = bic),
    cond_dist = cond_dist,
    orders = order
  )

  class(list_out) <- "megarch_fit"

  list_out

}

#'@export

print.megarch_fit <- function(x) {

  par_names <- names(x$pars)
  se <- unname(x$se)
  pars <- unname(x$pars)

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
    "*           Fitted MEGARCH          *\n",
    "*************************************\n",
    " \n",
    "Orders: (", x$orders[[1]], ", ", x$orders[[2]], ")\n",
    "Cond. distribution: ", x$cond_dist, "\n",
    " \n",
    "Fitted parameters:\n"
  )

  part2 <- paste0(
    " \n",
    "Information criteria:\n",
    "AIC: ", sprintf("%.4f", unname(x$inf_criteria[[1]])),
    ", BIC: ", sprintf("%.4f", unname(x$inf_criteria[[2]]))
  )

  cat(part1, fill = TRUE)
  print(df)
  cat(part2, fill = TRUE)


}

##################################

#### MLog-GARCH (complete)

sigt_mloggarch <- function(theta, x, dfun1, dfun2, p, q, include_mean,
                             extra_par, skew_par) {

  imean <- !include_mean
  mu <- c(0, theta[[1]])[[include_mean + 1]]
  step0 <- 2 - imean
  omega_sig <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q - 1

  phi <- theta[(step0 + 1):step1]

  if (q > 1) {
    psi <- c(1, theta[(step1 + 1):step2])
  } else {
    psi <- 1
  }
  kappa <- theta[[(step2 + 1)]]
  gamma <- theta[[step2 + 2]]
  if (extra_par) {
    shape <- theta[[step2 + 3]]
  } else {
    shape <- 0
  }
  if (skew_par) {
    skew <- theta[[step2 + 3 + extra_par]]
  } else {
    skew <- 0
  }

  alpha <- psi * kappa
  beta <- psi * gamma

  help_fun1 <- function(x, shape, skew, dfun) {
    (sign(x) * log(abs(x) + 1)) * dfun(x, shape, skew)
  }
  E1 <- integrate(help_fun1, -Inf, Inf, shape = shape, skew = skew,
                  dfun = dfun1,
                  abs.tol = 1e-15, subdivisions = 1000L)[[1]]

  help_fun2 <- function(x, E1, shape, skew, dfun) {
    abs(sign(x) * log(abs(x) + 1) - E1) * dfun(x, shape, skew)
  }
  E2 <- integrate(help_fun2, -Inf, Inf, E1 = E1, shape = shape,
                  skew = skew, dfun = dfun1,
                  abs.tol = 1e-15, subdivisions = 1000L)[[1]]

  Elnsig2 <- omega_sig
  omega <- Elnsig2 * (1- sum(phi))

  sigt <- c(sigt_mloggarchCpp(
    x = x,
    mu = mu,
    omega = omega,
    phi = phi,
    alpha = alpha,
    beta = beta,
    E1 = E1,
    E2 = E2,
    Elnsig2
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape)
}

# p, q are now model orders
goal_fun_mloggarch <- function(theta, x, dfun1, dfun2, p, q, include_mean,
                             extra_par, skew_par) {

  res_list <- sigt_mloggarch(theta, x, dfun1, dfun2, p, q,
              include_mean, extra_par, skew_par)

  sigt <- res_list$sigt
  mu <- res_list$mu
  shape <- res_list$shape
  skew <- res_list$skew

  -sum(log(dfun2(x, mu, sigt, shape, skew)))
}

mloggarchFit <- function(rt, order = c(1L, 1L),
          cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
          include_mean = TRUE, start_pars = NULL, LB = NULL, UB = NULL) {

  stopifnot("Both orders must be at least 1." = any(order >= 1))
  cond_dist <- match.arg(cond_dist)

  p <- order[[1]]
  q <- order[[2]]
  force(lookup_table)
  force(include_mean)
  fun1 <- lookup_table[["fun1"]][[cond_dist]]
  fun2 <- lookup_table[["fun2"]][[cond_dist]]

  skew_par <- lookup_table[["skew_par"]][[cond_dist]]
  extra_par <- lookup_table[["extra_par"]][[cond_dist]]

  goal_fun <- function(theta, x) {
    goal_fun_mloggarch(
      theta = theta,
      x = x,
      dfun1 = fun1,
      dfun2 = fun2,
      p = p,
      q = q,
      include_mean = include_mean,
      extra_par = extra_par,
      skew_par = skew_par
    )
  }

  vals <- lookup_table[["vals"]][[cond_dist]]
  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  if (is.null(start_pars)) {

    start_pars <- c(
      list(NULL, 0)[[include_mean + 1]],
      0,
      rep(0, p),
      rep(0, q - 1),
      c(0, 0),
      sval,
      list(NULL, 1)[[skew_par + 1]]
    )

  }

  if (is.null(LB)) {

    LB <- c(
      list(NULL, -Inf)[[include_mean + 1]],
      -Inf,
      rep(-2, p),
      rep(-2, q - 1),
      rep(-2, 2),
      lval,
      list(NULL, 1e-15)[[skew_par + 1]]
    )

  }

  if (is.null(UB)) {

    UB <- c(
      list(NULL, Inf)[[include_mean + 1]],
      Inf,
      rep(2, p),
      rep(2, q - 1),
      rep(2, 2),
      uval,
      list(NULL, Inf)[[skew_par + 1]]
    )

  }

  result <- suppressWarnings(Rsolnp::solnp(
    pars = start_pars,
    fun = goal_fun,
    LB = LB,
    UB = UB,
    control = list(trace = FALSE),
    x = rt
  ))

  pars <- result$pars

  llhood <- -tail(result$values, 1)
  sigt <- sigt_mloggarch(
    theta = pars,
    x = rt,
    dfun1 = fun1,
    dfun2 = fun2,
    p = p,
    q = q,
    include_mean = include_mean,
    extra_par = extra_par,
    skew_par = skew_par
  )$sigt

  mu <- c(0, pars[[1]])[[include_mean + 1]]
  etat <- (rt - mu) / sigt

  if (q > 1) {
    psi_names <- paste0("psi", 2:q)
  } else {
    psi_names <- NULL
  }

  par_names <- c(
    list(NULL, "mu")[[include_mean + 1]],
    "omega_sig",
    paste0("phi", 1:p),
    psi_names,
    "kappa",
    "gamma",
    list(NULL, c("shape", "df")[[(cond_dist %in% c("std", "sstd")) + 1]])[[extra_par + 1]],
    list(NULL, "skew")[[skew_par + 1]]
  )

  serrors <- sqrt(diag(solve(result$hessian)))
  names(pars) <- par_names
  names(serrors) <- par_names

  k <- length(pars)
  n <- length(rt)

  m2llhood <- -2 * llhood
  bic <- (m2llhood + k * log(n)) / n
  aic <- (m2llhood + 2 * k) / n

  list_out <- list(
    pars = pars,
    se = serrors,
    rt = rt,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = c("aic" = aic, "bic" = bic),
    cond_dist = cond_dist,
    orders = order
  )

  class(list_out) <- "mloggarch_fit"

  list_out

}

#'@export

print.mloggarch_fit <- function(x) {

  par_names <- names(x$pars)
  se <- unname(x$se)
  pars <- unname(x$pars)

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
    "*         Fitted MLog-GARCH         *\n",
    "*************************************\n",
    " \n",
    "Orders: (", x$orders[[1]], ", ", x$orders[[2]], ")\n",
    "Cond. distribution: ", x$cond_dist, "\n",
    " \n",
    "Fitted parameters:\n"
  )

  part2 <- paste0(
    " \n",
    "Information criteria:\n",
    "AIC: ", sprintf("%.4f", unname(x$inf_criteria[[1]])),
    ", BIC: ", sprintf("%.4f", unname(x$inf_criteria[[2]]))
  )

  cat(part1, fill = TRUE)
  print(df)
  cat(part2, fill = TRUE)


}

