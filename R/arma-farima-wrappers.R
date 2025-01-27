arma_fit <- function(x, mu, ar, ma, coef_inf) {
  c(fitted_arma_Cpp(x = x, mu = mu, ma = ma, ar = ar))
}

farima_fit <- function(x, mu, ar, ma, coef_inf) {
  c(fitted_farima_Cpp(x = x, mu = mu, coef_inf = coef_inf))
}
