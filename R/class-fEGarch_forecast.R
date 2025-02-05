setClass("fEGarch_forecast",
  slots = c(
    cmeans = "ANY",
    sigt = "ANY",
    model = "fEGarch_fit"
  ),
  prototype = list(
    cmeans = numeric(1),
    sigt = numeric(1),
    model = new("fEGarch_fit")
  )
)


fEGarch_forecast <- function(
    cmeans,
    sigt,
    model
) {

  new("fEGarch_forecast",
    cmeans = cmeans,
    sigt = sigt,
    model = model
  )

}
