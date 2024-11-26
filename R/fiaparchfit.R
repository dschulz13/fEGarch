sigt_fiaparch_R <- function(theta, x, dfun1, p, q,
              incl_mean, extra_par, skew_par,
              int_fun) {


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

  E_g_delta <- integrate(f = int_fun, lower = -Inf, upper = Inf, shape = shape, skew = skew,
                  pdf_fun = dfun1, delta = delta, gamma = gamma,
                  subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                  stop.on.error = FALSE)[[1]]

  n <- length(x)
  coef_inf <- -esemifar::farima_to_ar(ar = phi, ma = beta, d = d,
                                     max_i = n + 1000)[-1]

  sigt <- c(sigt_fiaparch_Cpp(
    x = x,
    mu = mu,
    coef_inf = coef_inf,
    gamma = gamma,
    delta = delta,
    E_g_delta = E_g_delta,
    omega_sig = omega_sig
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape)

}

goal_fun_creator_fiaparch <- function(
    theta,
    x,
    pdf_fun1,
    pdf_fun2,
    p,
    q,
    incl_mean,
    extra_par,
    skew_par,
    int_fun
) {

  res_list <- sigt_fiaparch_R(theta, x, dfun1, p, q,
              incl_mean, extra_par, skew_par,
              int_fun)

  sigt <- res_list$sigt
  mu <- res_list$mu
  skew <- res_list$skew
  shape <- res_list$shape

  out <- -sum(log(dfun2(x, mu, sigt, shape, skew)))

  out

}

fiaparchFit <- function(rt, orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
                        include_mean = TRUE) {

  cond_dist <- match.arg(cond_dist)

  pdf_fun1 <- lookup_table$fun1[[cond_dist]]
  pdf_fun2 <- lookup_table$fun2[[cond_dist]]
  int_fun <- function(x, delta, gamma, shape, skew, pdf_fun) {
    (abs(x) - gamma * x)^delta * pdf_fun(x, skew = skew, shape = shape)
  }

  incl_mean = include_mean

  #### continue here with skewpar and extra_par definitions

  goal_fun <- function(theta, x) {

    goal_fun_creator_fiaparch(
     theta = theta,
     x = x,
     pdf_fun1 = pdf_fun1,
     pdf_fun2 = pdf_fun2,
     p = p,
     q = q,
     incl_mean = incl_mean,
     extra_par,
     skew_par,
     int_fun
    )

  }




}
