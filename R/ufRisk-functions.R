setGeneric("measure_risk", function(object, measure = c("VaR", "ES"), level = 0.99, ...) {standardGeneric("measure_risk")})

#'@export
setMethod("measure_risk", "fEGarch_fit",
  function(object, measure = c("VaR", "ES"), level = 0.99, ...) {
    cond_d <- object@cond_dist

    par_name <- switch(
      cond_d,
      "norm" = "df",   # irrelevant
      "snorm" = "df",  # irrelevant
      "std" = "df",
      "sstd" = "df",
      "ged" = "shape",
      "sged" = "shape",
      "ald" = "P",
      "sald" = "P"
    )
    shape <- if (!(cond_d %in% c("norm", "snorm"))) {
      object@pars[[par_name]]
    } else {
      0
    }
    skew <- if (!(cond_d %in% c("norm", "std", "ged", "ald"))) {
      object@pars[["skew"]]
    } else {
      0
    }


    args <- list(
      level = level,
      dist = cond_d,
      skew = skew
    )
    args[[par_name]] <- shape

    if ("VaR" %in% measure) {
      q_VaR <- do.call(what = VaR_calc, args = args)
      names_q <- paste0("VaR", level)
      lVaR <- list()
      for (i in seq_along(q_VaR)) {
        lVaR[[names_q[[i]]]] <- zoo::zoo(zoo::coredata(object@cmeans - q_VaR[[i]] * object@sigt), order.by = time(object@sigt))
      }
    } else {
      lVaR <- list()
    }

    if ("ES" %in% measure) {
      q_ES <- do.call(what = ES_calc, args = args)
      names_q <- paste0("ES", level)
      lES <- list()
      for (i in seq_along(q_ES)) {
        lES[[names_q[[i]]]] <- zoo::zoo(zoo::coredata(object@cmeans - q_ES[[i]] * object@sigt), order.by = time(object@sigt))
      }
    } else {
      lES <- list()
    }

    list(
      measures = list(VaR = lVaR, ES = lES),
      observations = object@rt
    )

  }
)

#'@export
setMethod("measure_risk", "fEGarch_forecast",
  function(object, measure = c("VaR", "ES"), level = 0.99, ...) {

    object2 <- object@model
    object2@sigt <- object@sigt
    object2@cmeans <- object@cmeans

    func <- selectMethod("measure_risk", "fEGarch_fit")
    out <- func(object2)
    out[["observations"]] <- object2@test_obs
  }
)



