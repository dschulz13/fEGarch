#'Fitting Function for Models of the Broader EGARCH Family
#'
#'Use maximum-likelihood estimation to fit a model from the
#'broader EGARCH family to some observed time series.
#'
#'@param spec an S4 object of class \code{"egarch_type_spec"}
#'or \code{"loggarch_type_spec"} as returned by the various spec
#'functions of this package, for example \link[fEGarch]{fEGarch_spec}
#'or its wrappers like \link[fEGarch]{megarch_spec} or
#'\link[fEGarch]{mloggarch_spec}, among others.
#'@param rt the input time series to fit the model to ordered from
#'past to present; can also be a \code{"zoo"} class object.
#'@param start_pars a vector with starting parameters for the
#'optimization; must be of the same length as the output vector of
#'parameters; the default NULL uses internally saved default sets
#'of starting values.
#'@param LB a vector with lower boundaries for parameters; must be
#'of the same length as the output vector of
#'parameters; the default NULL uses internally saved default sets
#'of lower boundary values.
#'@param UB a vector with upper boundaries for parameters; must be
#'of the same length as the output vector of
#'parameters; the default NULL uses internally saved default sets
#'of upper boundary values.
#'
#'@export
#'
#'@importFrom Rcpp evalCpp
#'@useDynLib fEGarch
#'
#'@details
#'For details on the models, see \link[fEGarch]{fEGarch_spec}. The
#'specified models are then fitted using quasi maximum likelihood estimation,
#'where initial values are obtained from their theoretical expectations.
#'


fEGarch <- function(spec = egarch_spec(),
    rt,
    start_pars = NULL, LB = NULL, UB = NULL) {

  fEGarch_fit(spec, rt, start_pars, LB, UB)

}

