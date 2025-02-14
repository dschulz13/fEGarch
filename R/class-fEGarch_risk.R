
setClass("fEGarch_risk",
  slots = c(
    measures = "list",
    observations = "ANY"
  ),
  prototype = list(
    measures = list(1),
    observations = numeric(1)
  )
)

fEGarch_risk <- function(
    measures,
    observations
) {

  new("fEGarch_risk",
    measures = measures,
    observations = observations
  )

}

PoT_calc <- function(object, level = NULL) {
  VaR <- object@measures$VaR
  VaR_names <- if (is.null(level)) {
    names(VaR)
  } else {
    paste0("VaR", level)
  }
  out <- vector(mode = "list", length = length(VaR_names))
  names_adj <- paste0("PoT_", VaR_names)
  names(out) <- names_adj
  for (i in seq_along(VaR_names)) {
    selector <- object@observations < VaR[[VaR_names[[i]]]]
    out[[names_adj[[i]]]] <- object@observations[selector]
  }
  out
}

setMethod("plot", "fEGarch_risk", function(x, which = NULL, ...) {

  which <- check_which(which)

  if (which == "0") {
    return(invisible(NULL))
  }

  plot_fun <- switch(
    which,
    "returns" = plot_returns,
    "means" = plot_cmeans,
    "standard_deviations" = plot_sigt,
    "1" = plot_returns,
    "2" = plot_cmeans,
    "3" = plot_sigt
  )

  plot_fun(x, ...)

})
