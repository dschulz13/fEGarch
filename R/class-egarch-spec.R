setClass("egarch_spec",
  slots = c(
    model = "character",
    orders = "numeric",
    powers = "numeric",
    long_memo = "logical",
    modulus = "logical",
    model_type = "numeric",
    cond_dist = "character",
    include_mean = "logical"
  ),
  prototype = list(
    model = "eGARCH",
    orders = c(1, 1),
    powers = c(NA_real_, NA_real_),
    long_memo = NA,
    modulus = c(NA, NA),
    model_type = NA_real_,
    cond_dist = "norm",
    include_mean = TRUE
  )
)

#'EGARCH Family Model Specification
#'
#'Create an object with specifications for a model from
#'the broader EGARCH family.
#'
#'@param model_type a character value (either \code{"feGARCH"} or
#'\code{"logGARCH"}) that indicates the type of model to
#'implement (see Details for more information).
#'@param orders a two-element numeric vector with the model
#'orders; the first element is the order for the term based on
#'\eqn{\ln\left(\sigma_t^2\right)}, i.e. the log-transformed
#'conditional variance, while the second element is the order for
#'the innovation-based term (see Details below for more information).
#'@param long_memo a logical value that indicates whether the
#'long-memory version of the model should be considered or not.
#'@param cond_dist a character value stating the underlying
#'conditional distribution to consider; available are a normal
#'distribution (\code{"norm"}), a \eqn{t}-distribution
#'(\code{"std"}), a generalized error distribution
#'(\code{"ged"}) and the skewed versions of them
#'(\code{"snorm"}, \code{"sstd"}, \code{"sged"}).
#'@param include_mean a logical value indicating whether or
#'not to include the mean in the estimation procedure; for
#'\code{include_mean = FALSE}, the mean of the series
#'\code{rt} is fixed to zero and not being estimated.
#'@param powers a two-element numeric vector that states the
#'exponents in the power-transformations of the asymmetry and the
#'magnitude terms in that order (see Details for more information).
#'@param modulus a two-element logical vector indicating if the
#'innovations in the asymmetry and the magnitude terms (in that
#'order) should use a log-modulus transformation (see Details for
#'more information).
#'
#'@export
#'
#'@details
#'Let \eqn{\left\{r_t\right\}}, with \eqn{t \in \mathbb{Z}} as the
#'time index, be a theoretical time series that follows
#'\deqn{r_t=\mu+\sigma_t \eta_t \text{ with } \eta_t \sim \text{IID}(0,1), \text{ where}}
#'\deqn{\ln\left(\sigma_t^2\right)=\omega_{\sigma}+\phi^{-1}(B)\left(1-B\right)^{-d}\psi(B)g(\eta_{t-1}).}
#'Here, \eqn{\eta_t\sim\text{IID}(0,1)} means that the innovations
#'\eqn{\eta_t} are independent but identically distributed (iid) with mean zero
#'and variance one, whereas \eqn{\sigma_t > 0} are the conditional standard
#'deviations in \eqn{r_t}. Note that \eqn{\ln\left(\cdot\right)} denotes the natural logarithm.
#'Moreover,
#'\deqn{\phi(B) = 1-\sum_{i=1}^{p}\phi_{i}B^{i} \text{ and}}
#'\deqn{\psi(B) = 1+\sum_{j=1}^{q-1}\psi_{j}B^{j},}
#'where \eqn{B} is the backshift operator, are characteristic polynomials
#'with real coefficients \eqn{\phi_{0},\dots,\phi_{p},\psi_{0},\dots,\psi_{q-1}},
#'by fixing \eqn{\phi_{0}=\psi_{0}=1}, and without common roots. Furthermore,
#'we have \eqn{\mu = E\left(r_t\right)} as a real-valued parameter and the
#'long-memory parameter \eqn{d \in [0,0.5)}. If we consider
#'\eqn{\left\{g\left(\eta_t\right)\right\}} to be a series of iid variables
#'with mean zero and finite variance, then the real-valued parameter
#'\eqn{\omega_{\sigma}} is in fact
#'\eqn{\omega_{\sigma}=E\left[\ln\left(\sigma_t^2\right)\right]}.
#'
#'\eqn{g\left(\cdot\right)} can be defined in different ways. Commonly,
#'following a (FI)EGARCH specification, we have
#'\deqn{g\left(\eta_t\right)=\kappa \left\{g_a\left(\eta_t\right) - E\left[g_a\left(\eta_t\right)\right] \right\} + \gamma\left\{g_m\left(\eta_t\right)-E\left[ g_m\left(\eta_t\right)\right]\right\}}
#'with \eqn{g_a\left(\eta_t\right)=g_m\left(\eta_t\right)=\eta_t} and consequently
#'\eqn{E\left[g_a\left(\eta_t\right)\right] = E\left(\eta_t\right) = 0} and
#'\eqn{E\left[\left|g_m\left( \eta_t \right)\right|\right] = E\left(\left| \eta_t\right|\right)},
#'where \eqn{\kappa} and \eqn{\gamma} are two additional real-valued
#'parameters. Nonetheless, the functions \eqn{g_a\left(\cdot\right)} and
#'\eqn{g_m\left(\cdot\right)} can be any suitable transformation of
#'\eqn{\eta_t}.
#'
#'Generally, we consider two cases:
#'\deqn{g_{1,a}\left(\eta_t\right)=\text{sgn}\left(\eta_t\right)\left|\eta_t\right|^{p_a}/p_a \text{ and}}
#'\deqn{g_{2,a}\left(\eta_t\right)=\text{sgn}\left(\eta_t\right)\left(\left|\eta_t\right| + 1\right)^{p_a} / p_{a}}
#'whereas
#'\deqn{g_{1,m}\left(\eta_t\right)=\left|\eta_t\right|^{p_m}/p_m \text{ and}}
#'\deqn{g_{2,m}\left(\eta_t\right)=\left(\left|\eta_t\right| + 1\right)^{p_m} / p_{m}.}
#'Note that \eqn{\text{sgn}\left(\eta_t\right)} denotes the sign of
#'\eqn{\eta_t}. \eqn{g_{1,\cdot}} incorporates a power transformation
#'and \eqn{g_{2,\cdot}} a modulus transformation together with a power
#'transformation. The choices
#'\eqn{g_{1,a}} and \eqn{g_{2,a}} correspond to
#'setting the first element in \code{modulus} to \code{FALSE} or
#'\code{TRUE}, respectively, under
#'\code{model_type = "fegarch"}, where \eqn{p_a} can be selected
#'via the first element in \code{powers}. As a special case, for
#'\eqn{p_a = 0}, a log-transformation is employed and the division
#'through \eqn{p_a} is dropped, i.e.
#'\eqn{g_{1,a}\left(\eta_t\right)=\text{sgn}\left(\eta_t\right)\ln\left(\left|\eta_t\right|\right)}
#'and \eqn{g_{2,a}\left(\eta_t\right)=\text{sgn}\left(\eta_t\right)\ln\left(\left|\eta_t\right|+1\right)}
#'are employed for \eqn{p_a=0}. Completely analogous
#'thoughts hold for \eqn{g_{1,m}} and \eqn{g_{2,m}} and the second
#'elements in the arguments \code{modulus} and \code{powers}. The aforementioned
#'model family is a type I model selectable through
#'\code{model_type = "feGARCH"}.
#'
#'As additional specifications of a Log-GARCH and a FILog-GARCH, we redefine
#'\deqn{\psi(B) = 1+\sum_{j=1}^{q}\psi_{j}B^{j}}
#'now as a polynomial of order \eqn{q} and
#'\deqn{\ln\left(\sigma_t^2\right)=\omega_{\sigma}+\left[\phi^{-1}(B)(1-B)^{-d}\psi(B) - 1\right]g\left(\eta_t\right),}
#'where \eqn{g(\cdot)} is now defined through
#'\deqn{g\left(\eta_t\right)=\ln\left(\eta_t^2\right)-E\left[\ln\left(\eta_t^2\right)\right].}
#'Everything else is defined as before. This describes the
#'established Log-GARCH models as close relatives to the EGARCH family
#'(type II models; \code{model_type = "logGARCH"}).
#'
#'While the arguments \code{powers} and \code{modulus} are only relevant
#'under a type I model, i.e. for \code{model_type = "feGARCH"}, the arguments
#'\code{orders}, \code{long_memo}, \code{cond_dist} and
#'\code{include_mean} are meaningful for both \code{model_type = "feGARCH"}
#'and \code{model_type = "logGARCH"}, i.e. both under type I and II models.
#'The first element of the two-element vector orders is the order \eqn{p},
#'while the second element is the order \eqn{q}. Furthermore, for
#'\code{long_memo = TRUE}, the metnioned models are kept as they are, while
#'for \code{long_memo = FALSE} the parameter \eqn{d} is set to zero.
#'\code{cond_dist} controls the conditional distribution and for
#'\code{include_mean = FALSE}, \eqn{\mu} is not being estimated and fixed
#'to zero.
#'
#'@examples
#'spec <- fegarch_spec()
#'
#'



fegarch_spec <- function(
  model_type = c("feGARCH", "logGARCH"),
  orders = c(1, 1),
  long_memo = FALSE,
  cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
  include_mean = TRUE,
  powers = c(1, 1),
  modulus = c(FALSE, FALSE)
) {
  model <- match.arg(model)
  if (is.null(model) || is.na(model)) {
    model <- "eGARCH"
  }
  if (is.null(orders) || all(is.na(orders))) {
    orders <- c(1, 1)
  }
  if ((is.null(powers) || all(is.na(powers))) && !(model == "feGARCH")) {
    powers <- c(model_defaults[model, "powerTransf1"], model_defaults[model, "powerTransf2"])
  }
  if (is.null(long_memo) || is.na(long_memo)) {
    long_memo <- model_defaults[model, "longMemo"]
  }
  if (is.null(modulus) || all(is.na(modulus))) {
    modulus <- c(model_defaults[model, "modulus1"], model_defaults[model, "modulus2"])
  }
  if (is.null(model_type) || is.na(model_type)) {
    model_type <- model_defaults[model, "model_type"]
  }
  cond_dist <- match.arg(cond_dist)
  if (is.null(include_mean) || is.na(include_mean)) {
    include_mean <- TRUE
  }

  new("egarch_spec",
    model = model,
    orders = orders,
    powers = powers,
    long_memo = long_memo,
    modulus = modulus,
    model_type = model_type,
    cond_dist = cond_dist,
    include_mean = include_mean
  )
}

setValidity("egarch_spec", function(object){
  if (model_checker(object@model)) {
    "@model is not from the allowed models or longer than one"
  } else if (order_checker(object@orders)) {
    "@orders must be of length two and both elements must be at least one"
  } else if (power_checker(object@powers)) {
    "@powers should have length two"
  } else if (lmemo_checker(object@long_memo)) {
    "@long_memo must be of length one"
  } else if (modulus_checker(object@modulus)) {
    "@modulus must be of length two"
  } else if (model_type_checker(object@model_type)) {
    "@model_type must be either 1 or 2"
  } else if (distr_checker(object@cond_dist)) {
    "@cond_dist must be of length one and from the available model options"
  } else if (mean_checker(object@include_mean)) {
    "@include_mean must be of length one"
  } else {
    TRUE
  }
})

#'Accessors for Class \code{"egarch_spec"}
#'
#'Access and change elements in objects of class \code{"egarch_spec"}.
#'The method names represent the name of the element to access /
#'manipulate.
#'
#'@param x the input object or object to modify.
#'@param value the value to modify the object \code{x} with.
#'
#'@details
#'These methods are intended to be used for accessing or manipulating
#'individual elements of objects of class \code{"egarch_spec"}.
#'
#'@name spec_methods
#'
#'@return
#'These methods either return an object of class \code{"egarch_spec"}
#'or the corresponding element of an object of class
#'\code{"egarch_spec"}.
#'
#'@export
#'
#'@examples
#'test_obj <- egarch_spec()
#'orders(test_obj)
#'orders(test_obj) <- c(2, 1)
#'orders(test_obj)
#'
setMethod("model", "egarch_spec", function(x) {x@model})
#'@export
#'@rdname spec_methods
setMethod("orders", "egarch_spec", function(x) {x@orders})
#'@export
#'@rdname spec_methods
setMethod("powers", "egarch_spec", function(x) {x@powers})
#'@export
#'@rdname spec_methods
setMethod("long_memo", "egarch_spec", function(x) {x@long_memo})
#'@export
#'@rdname spec_methods
setMethod("modulus", "egarch_spec", function(x) {x@modulus})
#'@export
#'@rdname spec_methods
setMethod("model_type", "egarch_spec", function(x) {x@model_type})
#'@export
#'@rdname spec_methods
setMethod("cond_dist", "egarch_spec", function(x) {x@cond_dist})
#'@export
#'@rdname spec_methods
setMethod("include_mean", "egarch_spec", function(x) {x@include_mean})


#'@export
#'@rdname spec_methods
setMethod("model<-", "egarch_spec", function(x, value) {x@model <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("orders<-", "egarch_spec", function(x, value) {x@orders <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("powers<-", "egarch_spec", function(x, value) {x@powers <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("long_memo<-", "egarch_spec", function(x, value) {x@long_memo <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("modulus<-", "egarch_spec", function(x, value) {x@modulus <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("model_type<-", "egarch_spec", function(x, value) {x@model_type <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("cond_dist<-", "egarch_spec", function(x, value) {x@cond_dist <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("include_mean<-", "egarch_spec", function(x, value) {x@include_mean <- value; validObject(x); x})
