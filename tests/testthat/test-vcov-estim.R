vcov_args <- function(...) {
  utils::modifyList(
    list(
      goal_fun = function(theta, rt) sum(theta^2),
      pars = c(mu = 0.1, omega = 0.2),
      rt_core = c(0.1, -0.2, 0.3),
      scale_const = 1,
      incl_mean = TRUE,
      P_sel = 0.8,
      cond_d = "norm",
      par_names = c("mu", "omega"),
      model_type = "egarch",
      skip_vcov = TRUE
    ),
    list(...)
  )
}

run_vcov <- function(...) {
  do.call(compute_vcov, vcov_args(...))
}

test_that("skip_vcov returns an NA variance-covariance matrix", {

  pars <- c(mu = 0.1, omega = 0.2)

  result <- run_vcov(
    pars = pars,
    skip_vcov = TRUE
  )

  expect_equal(dim(result$vcov_mat), c(2L, 2L))
  expect_true(all(is.na(result$vcov_mat)))

  expect_equal(length(result$serrors), 2L)
  expect_true(all(is.na(result$serrors)))
  expect_named(result$serrors, c("mu", "omega"))

  expect_equal(unname(result$pars), unname(pars))
  expect_equal(rownames(result$vcov_mat), c("mu", "omega"))
  expect_equal(colnames(result$vcov_mat), c("mu", "omega"))
})

test_that("compute_vcov calculates and rescales the inverse Hessian", {

  pars <- c(mu = 0.1, omega = 0.2)
  rt <- c(0.1, -0.2, 0.3)

  H <- matrix(
    c(
      4, 1,
      1, 3
    ),
    nrow = 2
  )

  seen <- new.env()

  # Replace hessCalc() in the environment where compute_vcov() is defined.
  local_mocked_bindings(
    hessCalc = function(func, pars, rt) {
      seen$pars <- pars
      seen$rt <- rt
      H
    },
    .env = environment(compute_vcov)
  )

  result <- run_vcov(
    pars = pars,
    rt_core = rt,
    scale_const = 10,
    incl_mean = TRUE,
    skip_vcov = FALSE
  )

  # The mean parameter is rescaled by scale_const.
  I_mat <- diag(c(10, 1))
  expected_vcov <- I_mat %*% solve(H) %*% I_mat
  rownames(expected_vcov) <- c("mu", "omega")
  colnames(expected_vcov) <- c("mu", "omega")

  expect_equal(result$vcov_mat, expected_vcov)
  expect_equal(
    result$serrors,
    setNames(sqrt(diag(expected_vcov)), c("mu", "omega"))
  )

  # hessCalc() receives unnamed parameters and the original residuals.
  expect_identical(seen$pars, unname(pars))
  expect_identical(seen$rt, rt)

  expect_equal(unname(result$pars), unname(pars))
  expect_equal(rownames(result$vcov_mat), c("mu", "omega"))
  expect_equal(colnames(result$vcov_mat), c("mu", "omega"))
})


test_that("compute_vcov returns NA when the covariance diagonal is negative", {

  # The inverse of this Hessian has a negative diagonal element.
  H <- diag(c(-1, 1))

  local_mocked_bindings(
    hessCalc = function(func, pars, rt) H,
    .env = environment(compute_vcov)
  )

  expect_warning(
    result <- run_vcov(
      skip_vcov = FALSE
    ),
    regexp = "Unable to compute Hessian matrix"
  )

  expect_true(all(is.na(result$vcov_mat)))
  expect_true(all(is.na(result$serrors)))
})



test_that("APARCH and FIAPARCH correctly transform parameters and objective", {

  pars <- c(mu = 2, omega = 0.5, delta = 3)
  rt <- c(0.1, -0.2, 0.3)
  par_names <- c("mu", "omega", "delta")
  H <- diag(c(2, 4, 7))

  hess_seen <- list()
  adj_seen <- list()
  goal_seen <- list()

  local_mocked_bindings(
    hessCalc = function(func, pars, rt) {

      # Record the parameters passed to hessCalc().
      # Also evaluate the transformed objective function.
      hess_seen[[length(hess_seen) + 1L]] <<- list(
        pars = pars,
        rt = rt,
        objective = func(theta = pars, rt = rt)
      )

      H
    },
    .env = environment(compute_vcov)
  )

  goal_fun <- function(theta, rt) {

    goal_seen[[length(goal_seen) + 1L]] <<- list(
      theta = theta,
      rt = rt
    )

    sum(theta^2)
  }

  adj_fun_delta <- function(theta, scale_const, delta_pos, delta_inp) {

    adj_seen[[length(adj_seen) + 1L]] <<- list(
      theta = theta,
      scale_const = scale_const,
      delta_pos = delta_pos,
      delta_inp = delta_inp
    )

    3
  }

  for (model in c("aparch", "fiaparch")) {

    result <- run_vcov(
      goal_fun = goal_fun,
      pars = pars,
      rt_core = rt,
      scale_const = 10,
      incl_mean = TRUE,
      model_type = model,
      skip_vcov = FALSE,
      n_arma_pars = 1,
      sc_delta = 2,
      delta_pos = 3,
      adj_fun_delta = adj_fun_delta,
      delta_inp = "delta-input",
      par_names = c("mu", "omega", "delta")
    )

    sH <- solve(H)
    rownames(sH) <- c("mu", "omega", "delta")
    colnames(sH) <- c("mu", "omega", "delta")
    expect_equal(result$vcov_mat, sH)
    expect_equal(unname(result$pars), unname(pars))
  }

  # Both model types should have passed transformed parameters to hessCalc().
  expect_length(hess_seen, 2L)
  expect_equal(hess_seen[[1]]$pars, c(20, 1, 3))
  expect_equal(hess_seen[[2]]$pars, c(20, 1, 3))

  expect_identical(hess_seen[[1]]$rt, rt)
  expect_identical(hess_seen[[2]]$rt, rt)

  # The transformed objective converts the mean back to its original scale.
  expect_equal(adj_seen[[1]]$theta, c(2, 1, 3))
  expect_equal(adj_seen[[1]]$scale_const, 10)
  expect_equal(adj_seen[[1]]$delta_pos, 3)
  expect_equal(adj_seen[[1]]$delta_inp, "delta-input")

  expect_equal(goal_seen[[1]]$theta, c(2, 1/3, 3))
  expect_identical(goal_seen[[1]]$rt, rt)

  # The objective also includes n * log(scale_const).
  expected_objective <- sum(c(2, 1/3, 3)^2) + length(rt) * log(10)

  expect_equal(hess_seen[[1]]$objective, expected_objective)
  expect_equal(hess_seen[[2]]$objective, expected_objective)
})


test_that("ALD adds P_sel as the last parameter", {

  pars <- c(mu = 1, omega = 2)
  P_sel <- 0.75

  H <- diag(c(0.25, 0.5))
  expected_vcov <- solve(H)

  local_mocked_bindings(
    hessCalc = function(func, pars, rt) H,
    .env = environment(compute_vcov)
  )

  result <- run_vcov(
    pars = pars,
    P_sel = P_sel,
    cond_d = "ald",
    par_names = c("mu", "omega", "P_sel"),
    skip_vcov = FALSE
  )

  expected <- matrix(NA_real_, nrow = 3, ncol = 3)
  expected[1:2, 1:2] <- expected_vcov
  colnames(expected) <- c("mu", "omega", "P_sel")
  rownames(expected) <- c("mu", "omega", "P_sel")

  expect_equal(result$vcov_mat, expected)
  expect_equal(unname(result$pars), c(1, 2, P_sel))

  expect_equal(
    result$serrors,
    setNames(c(sqrt(4), sqrt(2), NA_real_), c("mu", "omega", "P_sel"))
  )

  expect_equal(rownames(result$vcov_mat), c("mu", "omega", "P_sel"))
  expect_equal(colnames(result$vcov_mat), c("mu", "omega", "P_sel"))
})


test_that("SALD inserts P_sel before the last parameter", {

  pars <- c(mu = 1, omega = 2, delta = 3)
  P_sel <- 0.75

  # Use a covariance matrix with distinct entries to identify the
  # parameter positions after insertion.
  expected_vcov <- matrix(
    c(
      4,   0.1, 0.2,
      0.1, 9,   0.3,
      0.2, 0.3, 16
    ),
    nrow = 3,
    byrow = TRUE
  )

  H <- solve(expected_vcov)

  local_mocked_bindings(
    hessCalc = function(func, pars, rt) H,
    .env = environment(compute_vcov)
  )

  result <- run_vcov(
    pars = pars,
    P_sel = P_sel,
    cond_d = "sald",
    par_names = c("mu", "omega", "P_sel", "delta"),
    skip_vcov = FALSE
  )

  # P_sel is inserted at position 3. The original third parameter
  # therefore moves to position 4 in the covariance matrix.
  expected <- matrix(NA_real_, nrow = 4, ncol = 4)
  original_positions <- c(1, 2, 4)

  expected[
    original_positions,
    original_positions
  ] <- expected_vcov

  colnames(expected) <- c("mu", "omega", "P_sel", "delta")
  rownames(expected) <- c("mu", "omega", "P_sel", "delta")

  expect_equal(result$vcov_mat, expected)

  expect_equal(
    unname(result$pars),
    c(1, 2, P_sel, 3)
  )

  expect_equal(
    result$serrors,
    setNames(c(2, 3, NA_real_, 4), c("mu", "omega", "P_sel", "delta"))
  )

  expect_equal(
    rownames(result$vcov_mat),
    c("mu", "omega", "P_sel", "delta")
  )

  expect_equal(
    colnames(result$vcov_mat),
    c("mu", "omega", "P_sel", "delta")
  )
})

