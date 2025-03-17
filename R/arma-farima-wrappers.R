arma_fit <- function(x, mu, ar, ma, coef_inf) {
  c(fitted_arma_Cpp(x = x, mu = mu, ma = ma, ar = ar))
}

# Using convolution techniques via FFT
farima_fit <- function(x, mu, ar, ma, coef_inf) {
  x_dm <- x - mu
  n <- length(x)
  np2 <- nextpow2(2 * n - 1)
  cl <- length(coef_inf)
  c_adj <- c(coef_inf, rep(0, np2 - cl))
  x_adj <- c(x_dm, rep(0, np2 - n))
  mu + Re(utils::head(stats::fft(stats::fft(c_adj) * stats::fft(x_adj), inverse = TRUE), n)) / np2
}

