# Generics
setGeneric("ljung_box_test", function(object, m_max = 20, adj = FALSE, weights = function(lag) {m <- max(lag); (m + 1 - lag) / m}, silent = FALSE, ...) {standardGeneric("ljung_box_test")})



#'@export
setMethod("ljung_box_test", "fEGarch_fit", function(object, m_max = 20, adj = FALSE, weights = function(lag) {m <- max(lag); (m + 1 - lag) / m}, silent = FALSE, type = c("single", "squared"), ...) {

  type <- match.arg(type)
  et <- object@etat

  fun <- list(
    "single" = function(et) {et},
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
  df[, 3] <- suppressWarnings(as.numeric(sprintf("%.4f", p_vals)))
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
    print(df)
  }
  invisible(df)

})


# Generics
setGeneric("sign_bias_test", function(object, silent = FALSE, ...) {standardGeneric("sign_bias_test")})



#'@export
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
  mat[, 1] <- as.numeric(sprintf("%.4f", p_all))
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
    print(mat)
  }
  invisible(mat)

})
