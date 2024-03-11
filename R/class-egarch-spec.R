setClass("egarch_spec",
  slots = c(
    model = "character",
    garch_order = "numeric",
    power_transf = "logical",
    long_memo = "logical",
    modulus = "logical",
    distribution = "character",
    include_mean = "logical"
  ),
  prototype = list(
    model = "eGARCH",
    garch_order = c(1, 1),
    power_transf = c(NA, NA),
    long_memo = NA,
    modulus = NA,
    distribution = "norm",
    include_mean = TRUE
  )
)

egarch_spec <- function(
  model = "eGARCH",
  garch_order = c(1, 1),
  power_transf = c(NA, NA),
  long_memo = NA,
  modulus = NA,
  distribution = "norm",
  include_mean = TRUE
) {
  if (is.null(model)) {
    model <- "eGARCH"
  }
  if (is.null(garch_order)) {
    garch_order <- c(1, 1)
  }
  if (is.null(power_transf) || all(is.na(power_transf))) {
    power_transf <- c(model_defaults[model, "powerTransf1"], model_defaults[model, "powerTransf2"])
  }
  if (is.null(long_memo) || is.na(long_memo)) {
    long_memo <- model_defaults[model, "longMemo"]
  }
  if (is.null(modulus) || is.na(modulus)) {
    modulus <- model_defaults[model, "modulus"]
  }
  if (is.null(include_mean)) {
    include_mean <- TRUE
  }

  new("egarch_spec",
    model = model,
    garch_order = garch_order,
    power_transf = power_transf,
    long_memo = long_memo,
    modulus = modulus,
    distribution = distribution,
    include_mean = include_mean
  )
}

setValidity("egarch_spec", function(object){
  if (model_checker(object@model)) {
    "@model is not from the allowed models or longer than one"
  } else if (order_checker(object@garch_order)) {
    "@garch_order must be of length two and at least one of its elements must be greater than zero"
  } else if (power_checker(object@power_transf)) {
    "@power_transf should have length two"
  } else if (lmemo_checker(object@long_memo)) {
    "@long_memo must be of length one"
  } else if (modulus_checker(object@modulus)) {
    "@modulus must be of length one"
  } else if (distr_checker(object@distribution)) {
    "@distribution must be of length one and from the available model options"
  } else if (mean_checker(object@include_mean)) {
    "@include_mean must be of length one"
  } else {
    TRUE
  }
})
