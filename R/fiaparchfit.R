sigt_fiaparch_R <- function(theta, x, p, q,
              incl_mean, extra_par, skew_par) {


  imean <- !incl_mean
  mu <- c(0, theta[[1]])[[incl_mean + 1]]
  step0 <- 2 - imean
  omega_sig <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q
  phi <- theta[(step0 + 1):step1]
  beta <- theta[(step1 + 1):step2]
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

  sigt <- c(sigt_fiaparch_Cpp(
    x = x,
    mu = mu,
    coef_inf = coef_inf,
    gamma = gamma,
    delta = delta,
    omega_sig = omega_sig
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape)

}

goal_fun_creator_fiaparch <- function(
    theta,
    x,
    pdf_fun2,
    p,
    q,
    incl_mean,
    extra_par,
    skew_par
) {

  res_list <- sigt_fiaparch_R(theta, x, p, q,
              incl_mean, extra_par, skew_par)

  sigt <- res_list$sigt
  mu <- res_list$mu
  skew <- res_list$skew
  shape <- res_list$shape

  out <- -sum(log(pdf_fun2(x, mu, sigt, shape, skew)))

  out

}


#'@export
fiaparchFit <- function(rt, orders = c(1, 1),
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
     extra_par,
     skew_par
    )

  }

  vals <- lookup_table[["vals"]][[cond_d]]

  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  d_spar <- 0.1
  d_low <- 0
  d_up <- 0.4999
  d_name <- "d"

  if (is.null(start_pars)) {

    start_pars <- c(
      list(NULL, 0)[[incl_mean + 1]],
      0.0001,
      rep(0.2, p),
      rep(0.1, q),
      0, 2,
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
      -0.9999, 1e-10,
      d_low,
      lval,
      list(NULL, 1e-15)[[skew_par + 1]]
    )

  }

  if (is.null(UB)) {

    UB <- c(
      list(NULL, 2)[[incl_mean + 1]],
      2,
      rep(2, p),
      rep(2, q),
      0.9999, 4,
      d_up,
      uval,
      list(NULL, Inf)[[skew_par + 1]]
    )

  }

  result <- Rsolnp::solnp(
    pars = start_pars,
    fun = goal_fun,
    LB = LB,
    UB = UB,
    control = list(trace = TRUE),
    rt = rt_core
  )

  pars <- result$pars
  llhood <- -tail(result$values, 1)

  sigt <- sigt_fiaparch_R(unname(pars), rt_core, p, q,
              incl_mean, extra_par, skew_par)$sigt

  mu <- c(0, pars[[1]])[[incl_mean + 1]]
  etat <- (rt - mu) / sigt

  beta_names <- paste0("beta", 1:q)

  par_names <- c(
    list(NULL, "mu")[[incl_mean + 1]],
    "omega_sig",
    paste0("phi", 1:p),
    beta_names,
    "gamma",
    "delta",
    d_name,
    list(NULL, c("shape", "df")[[(cond_d %in% c("std", "sstd")) + 1]])[[extra_par + 1]],
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

  if (inherits(rt, what = "zoo")) {
    tp <- as.Date(time(rt))
    sigt <- zoo::zoo(sigt, order.by = tp)
    etat <- zoo::zoo(etat, order.by = tp)
  }

  list_out <- list(
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
