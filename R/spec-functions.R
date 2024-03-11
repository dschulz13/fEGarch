e_spec <- function(
  variance.model = list(
    model = "eGARCH",
    garchOrder = c(1, 1),
    powerTransf = c(TRUE, TRUE),
    longMemo = FALSE
  ),
  mean.model = list(
    include.mean = TRUE
  ),
  distribution.model = "norm",
  start.pars = list(),
  fixed.pars = list(),
  ...
) {
  vm <- variance.model
  mm <- mean.model
  egarch_spec(
    model = vm$model,
    garch_order = vm$garchOrder,
    power_transf = vm$powerTransf,
    long_memo = vm$longMemo,
    include_mean = mm$include.mean,
    distribution = distribution.model
  )
}
