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
#'@return
#'An object of S4-class \code{"fEGarch_fit_egarch"} \code{"fEGarch_fit_loggarch"}
#'is returned depending on the selected model type in the model
#'specification. It contains the following elements.
#'\describe{
#'\item{\code{pars}:}{a named numeric vector with the parameter estimates.}
#'\item{\code{se}:}{a named numeric vector with the obtained standard errors in accordance with the parameter estimates.}
#'\item{\code{vcov_mat}:}{the variance-covariance matrix of the parameter estimates with named columns and rows.}
#'\item{\code{rt}:}{the input object \code{rt}; if \code{rt} was a \code{"zoo"} object, the formatting is kept.}
#'\item{\code{sigt}:}{the estimated conditional standard deviations; if \code{rt} was a \code{"zoo"} object, the formatting is also applied to \code{sigt}.}
#'\item{\code{etat}:}{the obtained residuals; if \code{rt} was a \code{"zoo"} object, the formatting is also applied to \code{etat}.}
#'\item{\code{orders}:}{a two-element numeric vector stating the considered model orders.}
#'\item{\code{cond_dist}:}{a character value stating the conditional distribution considered in the model fitting.}
#'\item{\code{long_memo}:}{a logical value stating whether or not long memory was considered in the model fitting.}
#'\item{\code{llhood}:}{the log-likelihood value obtained at the optimal parameter combination.}
#'\item{\code{inf_criteria}:}{a named two-element numeric vector with the corresponding AIC (first element) and BIC (second element) of the fitted model.}
#'\item{\code{powers}:}{a two-element numeric vector stating the powers considered for the asymmetry term (first element) and the magnitude term (second element); only exists, if a type I model was fitted.}
#'\item{\code{modulus}:}{a two-element logical vector stating the modulus transformations (\code{TRUE}, otherwise \code{FALSE}) considered for the asymmetry term (first element) and the magnitude term (second element); only exists, if a type I model was fitted.}
#'}
#'
#'@details
#'For details on the models, see \link[fEGarch]{fEGarch_spec}. The
#'specified models are then fitted using quasi maximum likelihood estimation,
#'where initial values are obtained from their theoretical expectations.
#'
#'The initialization of pre-sample values is always done according to the
#'theoretical expectations of the corresponding quantities; in the pre-sample
#'values of \eqn{g(\cdot)} are set to its expectation zero. All other possibly
#'required pre-sample values are set to their expectation using their formal
#'definition including unknown model parameters. These definitions are included
#'in the log-likelihood functions and are thus optimized over simultaneously.
#'
#'See the references section for sources on the EGARCH (Nelson, 1991),
#'FIEGARCH (Bollerslev and Mikkelsen, 1996),
#'Log-GARCH (Geweke, 1986; Pantual, 1986; Milhoj, 1987)
#'and FILog-GARCH (Feng et al., 2020) models.
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'spec <- fEGarch_spec()
#'model <- fEGarch(spec, rt)
#'model
#'
#'@references
#'\itemize{
#'\item{Bollerslev, T., & Mikkelsen, H. O. (1996). Modeling and pricing long memory in stock market volatility. Journal of Econometrics,
#'73(1), 151–184. DOI: 10.1016/0304-4076(95)01749-6.}
#'\item{Feng, Y., Beran, J., Ghosh, S., & Letmathe, S. (2020). Fractionally integrated Log-GARCH with application to value at risk and expected shortfall.
#'Working Papers CIE No. 137, Paderborn University, Center for International Economics.
#'URL: http://groups.uni-paderborn.de/wp-wiwi/RePEc/pdf/ciepap/WP137.pdf.}
#'\item{Geweke, J. (1986). Modeling the persistence of conditional variances: A comment. Econometric Reviews, 5(1),
#'57-61. DOI: 10.1080/07474938608800088.}
#'\item{Milhoj, A. (1987). A Multiplicative Parameterization of ARCH Models. University of Copenhagen, Denmark.}
#'\item{Nelson, D. B. (1991). Conditional Heteroskedasticity in Asset Returns: A New Approach. Econometrica,
#'59(2), 347–370. DOI: 10.2307/2938260.}
#'\item{Pantula, S. G. (1986). Modeling the persistence of conditional variances: A comment. Econometric Reviews, 5(1),
#'71-74. DOI: 10.1080/07474938608800089.}
#'}

fEGarch <- function(spec = egarch_spec(),
    rt,
    start_pars = NULL, LB = NULL, UB = NULL) {

  fEGarch_fit(spec, rt, start_pars, LB, UB)

}

