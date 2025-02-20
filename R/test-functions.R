# Generics
setGeneric("ljung_box", function(object, m_max = 10, ...) {standardGeneric("ljung_box")})



#'@export
setMethod("ljung_box", "fEGarch_fit", function(object, m_max = 10, type = c("single", "squared"), ...) {

  type <- match.arg(type)
  et <- object@etat

  fun <- list(
    "single" = function(et) {et},
    "squared" = function(et) {et^2}
  )[[type]]

  et_transf <- fun(et)

  cors <- c(autocorCpp(et_transf, lag_max = m_max))
  n <- length(et)

  stats <- (n * (n + 2)) * cumsum(cors^2 / (n - (1:m_max)))

  p_vals <- 1 - pchisq(stats, df = 1:m_max)


  m <- 1:m_max

  df <- matrix(NA, nrow = m_max, ncol = 2)
  df[, 1] <- m
  df[, 2] <- as.numeric(sprintf("%.4f", p_vals))
  rownames(df) <- rep("", m_max)
  colnames(df) <- c("m", "pval")

  out1 <- paste0(
    "H0: autocorrelations for lags 1 to m simultaneously equal to zero\n",
    "Mode: ", type, "\n"
  )
  cat(out1, fill = TRUE)
  print(df)

  invisible(NULL)

})


