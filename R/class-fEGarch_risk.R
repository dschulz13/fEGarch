
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

check_which_risk <- function(which, object) {

    if (is.null(which)) {

    VaRs <- if (!is.null(object@measures$VaR)) {
      paste0(substring(names(object@measures$VaR), 4), collapse = ", ")
    } else {
      ""
    }

    cat("\nThe following VaR levels were found (exit with 0):\n")
    cat(paste0("\n", VaRs, "\n"))

    cat("\n")

    which <- readline("Enter level to display: ")

    cat("\n")
    }

    which

}

setMethod("plot", "fEGarch_risk", function(x, which = NULL, ...) {

  which <- check_which_risk(which, x)

  VaR_name <- paste0("VaR", which)
  ES_name <- paste0("ES", which)

  VaR <- x@measures$VaR[[VaR_name]]
  ES <- x@measures$ES[[ES_name]]
  obs <- x@observations

  ES_check <- is.null(ES)

  if (which == "0") {
    return(invisible(NULL))
  }

  all_series <- cbind(
    obs, VaR, ES
  )

  tp <- c(time(obs))
  dim_s <- dim(all_series)
  all_series <- matrix(
    as.numeric(all_series),
    nrow = dim_s[[1]],
    ncol = dim_s[[2]]
  )

  input <- list(...)
  if (is.null(input[["type"]])) {
    input[["type"]] <- list("hll", "hl")[[ES_check + 1]]
  }
  if (is.null(input[["lty"]])) {
    input[["lty"]] <- list(c(1, 1, 1), c(1, 1))[[ES_check + 1]]
  }
  if (is.null(input[["col"]])) {
    input[["col"]] <- list(c(1, "gold2", 4), c(1, "gold2"))[[ES_check + 1]]
  }
  if (is.null(input[["xlab"]])) {
    input[["xlab"]] <- "Time"
  }
  if (is.null(input[["ylab"]])) {
    input[["ylab"]] <- "Value"
  }
  if (is.null(input[["main"]])) {
    input[["main"]] <- "The test series together with risk measures"
  }

  input[["x"]] <- tp
  input[["y"]] <- all_series

  do.call(matplot, args = input)

  abline(h = 0, lty = 3, col = "lightgray")

  idx <- obs < VaR

  obs <- as.numeric(obs)[idx]
  tp2 <- c(time(VaR))[idx]
  if (sum(idx) >= 1) {
    lines(tp2, obs, type = "h", col = "red", lwd = 1.5)
  }


})
