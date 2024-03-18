skew_density <- function(x, pdf, skew_par) {

  2 / (skew_par + 1 / skew_par) * (pdf(skew_par * x) * (-x >= 0) + pdf((1 / skew_par) * x) * (x >= 0))
}
