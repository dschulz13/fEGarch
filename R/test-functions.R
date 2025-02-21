# Generics
setGeneric("ljung_box", function(object, m_max = 20, adj = FALSE, weights = function(lag) {m <- max(lag); (m + 1 - lag) / m}, ...) {standardGeneric("ljung_box")})



#'@export
setMethod("ljung_box", "fEGarch_fit", function(object, m_max = 20, adj = FALSE, type = c("single", "squared"), ...) {

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

  out1 <- paste0(
    "\n*************************************************\n",
    "*Ljung-Box test for autocorrelation in residuals*\n",
    "*************************************************\n\n",
    "H0: autocorrelations for lags 1 to m simultaneously equal to zero\n",
    "Mode: ", type, "\n"
  )
  cat(out1, fill = TRUE)
  print(df)

  invisible(df)

})


