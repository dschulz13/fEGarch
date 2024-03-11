# For class "egarch_spec"

model_checker <- function(model) {
  !(model %in% c("eGARCH", "fieGARCH") & length(model) == 1)
}

order_checker <- function(garch_order) {
  !((length(garch_order) == 2) & (garch_order[[1]] > 0 || garch_order[[2]] > 0))
}

power_checker <- function(power_transf) {
  !(length(power_transf) == 2)
}

lmemo_checker <- function(long_memo) {
  !(length(long_memo) == 1)
}

modulus_checker <- function(modulus) {
  !(length(modulus) == 1)
}

distr_checker <- function(distribution) {
  !((distribution %in% c("norm", "snorm", "ged", "sged", "std", "sstd", "snig")) & (length(distribution) == 1))
}

mean_checker <- function(include_mean) {
  !(length(include_mean) == 1)
}
