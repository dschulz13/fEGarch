# Generics
setGeneric("ljung_box_test", function(object, m_max = 20, adj = TRUE, weights = function(lag) {m <- max(lag); (m + 1 - lag) / m}, silent = FALSE, ...) {standardGeneric("ljung_box_test")})


#'Weighted Ljung-Box Test for Autocorrelation
#'
#'Apply a (weighted) Ljung-Box test to check
#'the standardized residuals of a fitted model
#'from this package for remaining autocorrelation.
#'Two different options allow to check either the
#'simple residuals or the squared residuals.
#'
#'@param object an object \code{"fEGarch_fit"} as returned
#'by the fitting functions of this package, for example
#'by \code{\link{fEGarch}}.
#'@param m_max the maximum lag; tests will be conducted for 1 up to
#'\code{m_max}.
#'@param adj a logical value indicating whether or not to adjust
#'the degrees of freedom in the chi-square distribution according
#'to the number of estimated parameters.
#'@param weights a function with argument \code{lag} stating how
#'weights should be calculated.
#'@param silent a logical value reflecting whether or not test results
#'should be printed in a well-formatted manner to the console.
#'@param ... currently without use.
#'
#'@export
#'
#'@return
#'Returns a numeric matrix invisibly.
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'spec <- fEGarch_spec()
#'model <- fEGarch(spec, rt)
#'ljung_box_test(model)
#'
setMethod("ljung_box_test", "fEGarch_fit", function(object, m_max = 20, adj = TRUE, weights = function(lag) {m <- max(lag); (m + 1 - lag) / m}, silent = FALSE, type = c("simple", "squared"), ...) {

  type <- match.arg(type)
  et <- object@etat

  fun <- list(
    "simple" = function(et) {et},
    "squared" = function(et) {et^2}
  )[[type]]

  et_transf <- fun(et)

  cors <- c(autocorCpp(et_transf, lag_max = m_max))
  n <- length(et)

  stats <- rep(NA, m_max)

  for (i in 1:m_max) {
    w <- weights(lag = 1:i)
    stats[[i]] <- (n * (n + 2)) * sum((cors[1:i]^2 * w) / (n - (1:i)))
  }

  n_par <- if (adj) {
    length(object@pars)
  } else {
    0
  }

  dfr <- (1:m_max) - n_par
  n_df_pos <- sum(dfr > 0)
  p_vals <- c(rep(NA_real_, m_max - n_df_pos), 1 - pchisq(tail(stats, n_df_pos), df = tail(dfr, n_df_pos)))


  m <- 1:m_max

  df <- matrix(NA_real_, nrow = m_max, ncol = 3)
  df[, 1] <- m
  df[, 2] <- dfr
  df[, 3] <- suppressWarnings(p_vals)
  rownames(df) <- rep("", m_max)
  colnames(df) <- c("m", "df", "pval")

  if (!silent) {
    out1 <- paste0(
      "\n************************************************************\n",
      "* Weighted Ljung-Box test for autocorrelation in residuals *\n",
      "************************************************************\n\n",
      "H0: autocorrelations for lags 1 to m simultaneously equal to zero\n",
      "Mode: ", type, "\n"
    )
    cat(out1, fill = TRUE)
    df2 <- df
    df2[, 3] <- suppressWarnings(as.numeric(sprintf("%.4f", p_vals)))
    print(df2)
  }
  invisible(df)

})


# Generics
setGeneric("sign_bias_test", function(object, silent = FALSE, ...) {standardGeneric("sign_bias_test")})

#'Sign Bias Test
#'
#'Apply a sign bias test to check
#'the standardized residuals of a fitted model
#'from this package for remaining significant
#'sign effects.
#'
#'@param object an object \code{"fEGarch_fit"} as returned
#'by the fitting functions of this package, for example
#'by \code{\link{fEGarch}}.
#'@param silent a logical value reflecting whether or not test results
#'should be printed in a well-formatted manner to the console.
#'@param ... currently without use.
#'
#'@export
#'
#'@return
#'Returns a numeric matrix invisibly.
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'spec <- fEGarch_spec()
#'model <- fEGarch(spec, rt)
#'sign_bias_test(model)
#'
setMethod("sign_bias_test", "fEGarch_fit", function(object, silent = FALSE, ...) {

  et <- zoo::coredata(object@etat)
  et_h <- head(et, -1)

  It <- as.numeric(et_h < 0)
  It2 <- 1 - It

  V1 <- It
  V2 <- It * et_h
  V3 <- It2 * et_h
  Y <- et[-1]^2

  n <- length(Y)

  reg <- stats::lm(Y ~ V1 + V2 + V3)

  sum_stats <- summary(reg)
  R2 <- sum_stats$r.squared

  pvals_reg <- unname(sum_stats$coefficients[, 4, drop = TRUE][-1])

  p_joint <- 1 - pchisq(n * R2, df = 3)
  p_all <- c(pvals_reg, p_joint)

  mat <- matrix(NA, ncol = 1, nrow = 4)
  rownames(mat) <- c(
    "Sign bias",
    "Negative sign bias",
    "Positive sign bias",
    "Joint"
  )
  mat[, 1] <- p_all
  colnames(mat) <- "pval"

  if (!silent) {
    out1 <- paste0(
      "\n******************\n",
      "* Sign bias test *\n",
      "******************\n\n",
      "H0 (sign bias):      negative shocks have no significant effect\n",
      "                     on volatility beyond symmetric model\n",
      "H0 (neg. sign bias): magnitude of negative shocks has no\n",
      "                     additional impact on volatility\n",
      "H0 (pos. sign bias): magnitude of positive shocks has no\n",
      "                     additional impact on volatility\n",
      "H0 (joint):          all three of the above simultaneously\n"
    )
    cat(out1, fill = TRUE)
    mat2 <- mat
    mat2[, 1] <- as.numeric(sprintf("%.4f", p_all))
    print(mat2)
  }
  invisible(mat)

})

# Generics
setGeneric("goodn_of_fit_test", function(object, n_bins = c(20, 30, 40, 50), adj = TRUE, silent = FALSE, ...) {standardGeneric("goodn_of_fit_test")})


#'Adjusted Pearson Goodness-of-Fit Test for Standardized Model Residuals
#'
#'Consider a probability integral transform on the magnitudes of the
#'standardized residuals of a fitted model and apply a
#'Pearson goodness-of-fit test to the binned data using a
#'selection of predefined number of bins.
#'
#'@param object an object of class \code{"fEGarch_fit"} as returned
#'by the fitting functions of this package like for example
#'\code{\link{fEGarch}}.
#'@param n_bins a numeric vector giving the number of bins to use.
#'@param adj a logical value indicating whether or not to adjust
#'the number of degrees of freedom using the number of estimated
#'model parameters.
#'@param silent a logical indicating whether or not to print
#'the test results in a well-formatted manner to the console.
#'@param ... currently without purpose.
#'
#'@details
#'Use a probability integral transform on the magnitudes of the
#'standardized residuals of a fitted model. This is then the
#'basis to conduct a Pearson goodness-of-fit chi-square test.
#'
#'@return
#'Returns a numeric matrix invisibly.
#'
#'@export
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'spec <- fEGarch_spec()
#'model <- fEGarch(spec, rt)
#'goodn_of_fit_test(model)
#'
setMethod("goodn_of_fit_test", "fEGarch_fit", function(object, n_bins = c(10, 20, 30, 40, 50, 60), adj = TRUE, silent = FALSE, ...) {

  model <- object
  pars <- model@pars
  cond_d <- model@cond_dist
  n_par <- length(pars)

  int_fun_s <- lookup_table[["fun1"]][[cond_d]]

  int_fun <- function(x, int_fun_s, shape, skew) {
    xa <- abs(x)
    int_fun_s(x = xa, shape = shape, skew = skew) + int_fun_s(x = -xa, shape = shape, skew = skew)
  }

  args <- list(
    "dfun" = int_fun,
    "int_fun_s" = int_fun_s
  )

  args[["shape"]] <- if (!(cond_d %in% c("norm", "snorm"))) {
    par_name <- switch(
      cond_d,
      "std" = "df",
      "sstd" = "df",
      "ged" = "shape",
      "sged" = "shape",
      "ald" = "P",
      "sald" = "P"
    )
    pars[[par_name]]
  } else {
    0
  }

  args[["skew"]] <- if (cond_d %in% c("snorm", "sstd", "sged", "sald")) {
    pars[["skew"]]
  } else {
    0
  }

  lb <- length(n_bins)
  out_mat <- matrix(NA, nrow = lb, ncol = 3)
  rownames(out_mat) <- rep("", lb)
  colnames(out_mat) <- c("n_bins", "df", "pval")
  et <- zoo::coredata(model@etat)
  et_a <- abs(et)
  n <- length(et)
  adjust <- n_par * adj

  for (i in 1:lb) {

    nb <- n_bins[[i]]
    one_div_nb <- 1 / nb
    p <- seq(from = one_div_nb, to = 1 - one_div_nb, by = one_div_nb)

    args[["p"]] <- p

    quantiles <- -rev(do.call(inv_cdf_fun, args = args))   # Get quantiles corresponding to class bounds

    np1 <- length(p) + 1
    frequ_et <- rep(NA, np1)
    q_adj <- c(-Inf, quantiles, Inf)

    for (j in seq_along(frequ_et)) {
      frequ_et[[j]] <- sum(et_a >= q_adj[[j]] & et_a < q_adj[[j + 1]])
    }

    expec <- n * one_div_nb   # expected frequencies in each class

    stat <- sum((frequ_et - expec)^2 / expec)   # test statistic

    df <- nb - 1 - adjust

    pval <- if (df >= 1) {
      1 - pchisq(stat, df = df)
    } else {
      NA
    }

    out_mat[i, ] <- c(
      nb,
      df,
      pval
    )

  }



  if (!silent) {
    out1 <- paste0(
    "\n************************************\n",
      "* Adjusted Pearson goodness-of-fit *\n",
      "*         test on residuals        *\n",
      "************************************\n\n",
      "H0: standardized model residuals align with\n",
      "    assumed innovation distribution\n"
    )
    cat(out1, fill = TRUE)
    mat2 <- out_mat
    mat2[, 3] <- as.numeric(sprintf("%.4f", out_mat[, 3]))
    print(mat2)
  }

  invisible(out_mat)
})
