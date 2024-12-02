setClass("base_garch_spec",
  slots = c(
    orders = "numeric",
    long_memo = "logical",
    cond_dist = "character",
    include_mean = "logical"
  ),
  prototype = list(
    orders = c(1, 1),
    long_memo = FALSE,
    cond_dist = "norm",
    include_mean = TRUE
  )
)

setValidity("base_garch_spec", function(object){
  if (order_checker(object@orders)) {
    "@orders must be of length two and both elements must be at least one"
  } else if (lmemo_checker(object@long_memo)) {
    "@long_memo must be of length one"
  } else if (distr_checker(object@cond_dist)) {
    "@cond_dist must be of length one and from the available model options"
  } else if (mean_checker(object@include_mean)) {
    "@include_mean must be of length one"
  } else {
    TRUE
  }
})

base_garch_spec <- function(
    orders = c(1, 1),
    long_memo = TRUE,
    cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
    include_mean = TRUE
) {
  cond_dist <- match.arg(cond_dist)

  new("base_garch_spec",
    orders = orders,
    long_memo = long_memo,
    cond_dist = cond_dist,
    include_mean = include_mean
  )

}

setClass("loggarch_type_spec",
  contains = "base_garch_spec"
)

setClass("egarch_type_spec",
  contains = "base_garch_spec",
  slots = c(
    powers = "numeric",
    modulus = "logical"
  ),
  prototype = list(
    powers = c(NA_real_, NA_real_),
    modulus = c(NA, NA)
  )
)

setValidity("egarch_type_spec", function(object){
  if (power_checker(object@powers)) {
    "@powers should have length two"
  } else if (modulus_checker(object@modulus)){
    "@modulus must be a logical two-element vector"
  } else {
    TRUE
  }
})

#'@export
egarch_type_spec <- function(
  orders = c(1, 1),
  long_memo = TRUE,
  cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
  include_mean = TRUE,
  powers = c(1, 1),
  modulus = c(FALSE, FALSE)
) {

  cond_dist <- match.arg(cond_dist)

  new("egarch_type_spec",
    orders = orders,
    long_memo = long_memo,
    cond_dist = cond_dist,
    include_mean = include_mean,
    powers = powers,
    modulus = modulus
  )

}

#'@export
loggarch_type_spec <- function(
  orders = c(1, 1),
  long_memo = TRUE,
  cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
  include_mean = TRUE
) {

  cond_dist <- match.arg(cond_dist)

  new("loggarch_type_spec",
    orders = orders,
    long_memo = long_memo,
    cond_dist = cond_dist,
    include_mean = include_mean
  )

}

#'General EGARCH Family Model Specification
#'
#'Create an object with specifications for a model from
#'the broader EGARCH family.
#'
#'@param model_type a character value (either \code{"eGARCH"} or
#'\code{"logGARCH"}) that indicates the type of model to
#'implement (see Details for more information).
#'@param orders a two-element numeric vector with the model
#'orders; the first element is the order \eqn{p} for the term based on
#'\eqn{\ln\left(\sigma_t^2\right)}, i.e. the log-transformed
#'conditional variance, while the second element is the order \eqn{q} for
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
#'order) should use a modulus transformation (see Details for
#'more information).
#'
#'@export
#'
#'@details
#'Let \eqn{\left\{r_t\right\}}, with \eqn{t \in \mathbb{Z}} as the
#'time index, be a theoretical time series that follows
#'\deqn{r_t=\mu+\sigma_t \eta_t \text{ with } \eta_t \sim \text{IID}(0,1), \text{ where}}
#'\deqn{\ln\left(\sigma_t^2\right)=\omega_{\sigma}+\theta(B)g(\eta_{t-1}).}
#'Here, \eqn{\eta_t\sim\text{IID}(0,1)} means that the innovations
#'\eqn{\eta_t} are independent but identically distributed (iid) with mean zero
#'and variance one, whereas \eqn{\sigma_t > 0} are the conditional standard
#'deviations in \eqn{r_t}. Note that \eqn{\ln\left(\cdot\right)} denotes the natural logarithm.
#'Moreover, \eqn{B} is the backshift operator and
#'\eqn{\theta(B) = 1 +\sum_{i=1}^{\infty}\theta_i B^{i}}, where
#'\eqn{\theta_i}, \eqn{i=1,2,\dots}, are real-valued coefficients.
#'\eqn{g\left(\eta_{t-1}\right)} is a suitable function in \eqn{\eta_{t-1}}.
#'Generally, \eqn{\left\{g\left(\eta_{t}\right)\right\}} should be an iid
#'zero-mean sequence with finite variance.
#'We have \eqn{\mu = E\left(r_t\right)} as a real-valued parameter.
#'The real-valued parameter
#'\eqn{\omega_{\sigma}} is in fact
#'\eqn{\omega_{\sigma}=E\left[\ln\left(\sigma_t^2\right)\right]}.
#'This previous set of
#'equations defines the broader family of EGARCH models, from which
#'subtypes are described in the following that depend on the choice of
#'\eqn{g}.
#'
#'\eqn{\textbf{Type I}:}
#'
#'We have \eqn{\theta(B) = \phi^{-1}(B)(1-B)^{-d}\psi(B)}, where
#'\deqn{\phi(B) = 1-\sum_{i=1}^{p}\phi_{i}B^{i} \text{ and}}
#'\deqn{\psi(B) = 1+\sum_{j=1}^{q-1}\psi_{j}B^{j},}
#'are characteristic polynomials
#'with real coefficients \eqn{\phi_{0},\dots,\phi_{p},\psi_{0},\dots,\psi_{q-1}},
#'by fixing \eqn{\phi_{0}=\psi_{0}=1}, and without common roots. Furthermore,
#'the
#'long-memory parameter is \eqn{d \in [0,0.5)}.
#'
#'\eqn{g\left(\cdot\right)} can be defined in different ways.
#'Following a type I specification (\code{model_type = "eGARCH"}), we have
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
#'\code{model_type = "eGARCH"}, where \eqn{p_a} can be selected
#'via the first element in \code{powers}. As a special case, for
#'\eqn{p_a = 0}, a log-transformation is employed and the division
#'through \eqn{p_a} is dropped, i.e.
#'\eqn{g_{1,a}\left(\eta_t\right)=\text{sgn}\left(\eta_t\right)\ln\left(\left|\eta_t\right|\right)}
#'and \eqn{g_{2,a}\left(\eta_t\right)=\text{sgn}\left(\eta_t\right)\ln\left(\left|\eta_t\right|+1\right)}
#'are employed for \eqn{p_a=0}. Completely analogous
#'thoughts hold for \eqn{g_{1,m}} and \eqn{g_{2,m}} and the second
#'elements in the arguments \code{modulus} and \code{powers}. The aforementioned
#'model family is a type I model selectable through
#'\code{model_type = "eGARCH"}.
#'
#'\eqn{\textbf{Type II}:}
#'
#'As additional specifications of a Log-GARCH and a FILog-GARCH, which belong
#'to the broader EGARCH family, we redefine
#'\deqn{\psi(B) = 1+\sum_{j=1}^{q}\psi_{j}B^{j}}
#'now as a polynomial of order \eqn{q} and
#'\deqn{\ln\left(\sigma_t^2\right)=\omega_{\sigma}+\left[\phi^{-1}(B)(1-B)^{-d}\psi(B) - 1\right]\xi_t,}
#'where
#'\deqn{\xi_t=\ln\left(\eta_t^2\right)-E\left[\ln\left(\eta_t^2\right)\right].}
#'Everything else is defined as before. Since
#'\deqn{\phi^{-1}(B)(1-B)^{-d}\psi(B) - 1=\sum_{i=1}^{\infty}\gamma_i B^{i}=\gamma_1 B\left[\sum_{i=1}^{\infty}(\gamma_i/\gamma_1)B^{i-1}\right]=\gamma_1 B\left[1+\sum_{i=1}^{\infty}\theta_i B^{i}\right]=\theta(B)\gamma_1 B,}
#'where \eqn{\theta_i = \gamma_{i+1}/\gamma_{1}}, \eqn{i=0,1,\dots}, and by defining
#'\deqn{g\left(\eta_{t-1}\right)=\gamma_1\left\{\ln\left(\eta_{t-1}^2\right)-E\left[\ln\left(\eta_{t-1}^2\right)\right]\right\},}
#'the equation of \eqn{\ln\left(\sigma_t^2\right)} can be stated to be
#'\deqn{\ln\left(\sigma_t^2\right)=\omega_{\sigma}+\theta(B)g(\eta_{t-1})}
#'as in the broad EGARCH family at the very beginning.
#'
#'This describes the
#'established Log-GARCH models as part of the broad EGARCH family
#'(type II models; \code{model_type = "logGARCH"}).
#'
#'\eqn{\textbf{General information}:}
#'
#'While the arguments \code{powers} and \code{modulus} are only relevant
#'under a type I model, i.e. for \code{model_type = "eGARCH"}, the arguments
#'\code{orders}, \code{long_memo}, \code{cond_dist} and
#'\code{include_mean} are meaningful for both \code{model_type = "eGARCH"}
#'and \code{model_type = "logGARCH"}, i.e. both under type I and II models.
#'The first element of the two-element vector orders is the order \eqn{p},
#'while the second element is the order \eqn{q}. Furthermore, for
#'\code{long_memo = TRUE}, the metnioned models are kept as they are, while
#'for \code{long_memo = FALSE} the parameter \eqn{d} is set to zero.
#'\code{cond_dist} controls the conditional distribution and for
#'\code{include_mean = FALSE}, \eqn{\mu} is not being estimated and fixed
#'to zero.
#'
#'See also the closely related spec-functions that immediately create
#'specifications of specific submodels of the broad EGARCH family. These
#'functions are \code{egarch_spec()}, \code{fiegarch_spec()},
#'\code{loggarch_spec()},\code{filoggarch_spec()}, \code{megarch_spec()},
#'\code{mloggarch_spec()} and \code{mafiloggarch_spec()}, which are all
#'wrappers for \code{fEGarch_spec()}.
#'
#'See the references section for sources on the EGARCH (Nelson, 1991),
#'FIEGARCH (Bollerslev and Mikkelsen, 1996),
#'Log-GARCH (Geweke, 1986; Pantual, 1986; Milhoj, 1987)
#'and FILog-GARCH (Feng et al., 2020) models.
#'
#'@return
#'An object of either class \code{"egarch_type_spec"} or
#'\code{"loggarch_type_spec"} is returned, depending on the
#'choice for the input argument \code{model_type}.
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
#'
#'@examples
#'# EGARCH(1, 1) with cond. normal distribution
#'spec1 <- fEGarch_spec()
#'# EGARCH(2, 1) with cond. t-distribution
#'spec2 <- fEGarch_spec(orders = c(2, 1), cond_dist = "std")
#'# FIEGARCH(1, 1) with cond. normal distribution
#'spec3 <- fEGarch_spec(long_memo = TRUE)
#'# MEGARCH(1, 1) with cond. generalized error distribution
#'spec4 <- fEGarch_spec(modulus = c(TRUE, FALSE), powers = c(0, 1))
#'# Some unnamed specification
#'spec5 <- fEGarch_spec(
#'  model_type = "eGARCH",
#'  orders = c(1, 1),
#'  long_memo = TRUE,
#'  cond_dist = "std",
#'  include_mean = TRUE,
#'  powers = c(0.25, 0.75),
#'  modulus = c(TRUE, FALSE)
#')
#'

fEGarch_spec <- function(
  model_type = c("eGARCH", "logGARCH"),
  orders = c(1, 1),
  long_memo = FALSE,
  cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
  include_mean = TRUE,
  powers = c(1, 1),
  modulus = c(FALSE, FALSE)
) {
  model_type <- match.arg(model_type)
  if (is.null(model_type) || is.na(model_type)) {
    model_type <- "eGARCH"
  }
  if (is.null(orders) || all(is.na(orders))) {
    orders <- c(1, 1)
  }
  if (is.null(long_memo) || is.na(long_memo)) {
    long_memo <- FALSE
  }
  if ((is.null(modulus) || all(is.na(modulus))) && (model_type == "eGARCH")) {
    modulus <- c(FALSE, FALSE)
  }
  cond_dist <- match.arg(cond_dist)
  if (is.null(include_mean) || is.na(include_mean)) {
    include_mean <- TRUE
  }
  if ((is.null(powers) || all(is.na(powers))) && (model_type == "eGARCH")) {
    powers <- c(1, 1)
  }

  if (model_type == "eGARCH") {
    out <- egarch_type_spec(
      orders = orders,
      long_memo = long_memo,
      cond_dist = cond_dist,
      include_mean = include_mean,
      powers = powers,
      modulus = modulus
    )
  } else if (model_type == "logGARCH") {
    out <- loggarch_type_spec(
      orders = orders,
      long_memo = long_memo,
      cond_dist = cond_dist,
      include_mean = include_mean
    )
  }

  out

}



#'Accessors for Classes \code{"base_garch_spec"} and \code{"egarch_spec"}
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
#'These methods either return an object of class \code{"egarch_type_spec"}
#'or \code{"loggarch_type_spec"} or the corresponding element of
#'object of such class objects.
#'
#'@export
#'
#'@examples
#'test_obj <- egarch_spec()
#'orders(test_obj)
#'orders(test_obj) <- c(2, 1)
#'orders(test_obj)
#'
setMethod("orders", "base_garch_spec", function(x) {x@orders})
#'@export
#'@rdname spec_methods
setMethod("powers", "egarch_type_spec", function(x) {x@powers})
#'@export
#'@rdname spec_methods
setMethod("long_memo", "base_garch_spec", function(x) {x@long_memo})
#'@export
#'@rdname spec_methods
setMethod("modulus", "egarch_type_spec", function(x) {x@modulus})
#'@export
#'@rdname spec_methods
setMethod("cond_dist", "base_garch_spec", function(x) {x@cond_dist})
#'@export
#'@rdname spec_methods
setMethod("include_mean", "base_garch_spec", function(x) {x@include_mean})

#'@export
#'@rdname spec_methods
setMethod("orders<-", "base_garch_spec", function(x, value) {x@orders <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("powers<-", "egarch_type_spec", function(x, value) {x@powers <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("long_memo<-", "base_garch_spec", function(x, value) {x@long_memo <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("modulus<-", "egarch_type_spec", function(x, value) {x@modulus <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("cond_dist<-", "base_garch_spec", function(x, value) {x@cond_dist <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("include_mean<-", "base_garch_spec", function(x, value) {x@include_mean <- value; validObject(x); x})

#'EGARCH Family Submodel Specification
#'
#'Wrappers of \code{fEGarch_spec()} that create
#'specifications of specific submodels of the broad
#'EGARCH family.
#'
#'@param orders a two-element numeric vector with the model
#'orders.
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
#'
#'@export
#'
#'@details
#'Available are shortcut specification functions for
#'EGARCH \code{egarch_spec()}, FIEGARCH \code{fiegarch_spec()},
#'MEGARCH \code{megarch_spec()}, Log-GARCH \code{loggarch_spec()},
#'FILog-GARCH \code{filoggarch_spec()}, MLog-GARCH \code{mloggarch_spec()} and
#'MAFILog-GARCH \code{mafiloggarch_spec()}.
#'
#'The following descriptions are following the descriptions in the
#'documentation of the more general \code{fEGarch_spec()}. Please go there
#'first to understand the following descriptions on the arguments of
#'\code{fEGarch_spec()} to obtain these wrappers.
#'
#'\eqn{\textbf{EGARCH:}}
#'
#'\code{model_type = "eGARCH"}, \code{long_memo = FALSE},
#'\code{powers = c(1, 1)}, \code{modulus = c(FALSE, FALSE)}
#'
#'\eqn{\textbf{FIEGARCH:}}
#'
#'\code{model_type = "eGARCH"}, \code{long_memo = TRUE},
#'\code{powers = c(1, 1)}, \code{modulus = c(FALSE, FALSE)}
#'
#'\eqn{\textbf{MEGARCH:}}
#'
#'\code{model_type = "eGARCH"}, \code{long_memo = FALSE},
#'\code{powers = c(0, 1)}, \code{modulus = c(TRUE, FALSE)}
#'
#'\eqn{\textbf{Log-GARCH:}}
#'
#'\code{model_type = "logGARCH"}, \code{long_memo = FALSE}
#'
#'\eqn{\textbf{FILog-GARCH:}}
#'
#'\code{model_type = "logGARCH"}, \code{long_memo = TRUE}
#'
#'\eqn{\textbf{MLog-GARCH:}}
#'
#'\code{model_type = "eGARCH"}, \code{long_memo = FALSE},
#'\code{powers = c(0, 0)}, \code{modulus = c(TRUE, TRUE)}
#'
#'\eqn{\textbf{MAFILog-GARCH:}}
#'
#'\code{model_type = "eGARCH"}, \code{long_memo = TRUE},
#'\code{powers = c(0, 0)}, \code{modulus = c(TRUE, TRUE)}
#'
#'@name submodel-specs
#'
#'@return
#'Depending on the spec-fun function, either an object of class
#'\code{"egarch-type-spec"} or \code{"loggarch-type-spec"} is returned.
#'
#'@examples
#'spec <- megarch_spec(cond_dist = "std")
#'

megarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
    include_mean = TRUE
) {
  fEGarch_spec(
    model_type = "eGARCH",
    orders = orders,
    long_memo = FALSE,
    cond_dist = cond_dist,
    include_mean = include_mean,
    powers = c(0, 1),
    modulus = c(TRUE, FALSE)
  )
}

#'@export
#'@rdname submodel-specs
egarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
    include_mean = TRUE
) {
  fEGarch_spec(
    model_type = "eGARCH",
    orders = orders,
    long_memo = FALSE,
    cond_dist = cond_dist,
    include_mean = include_mean,
    powers = c(1, 1),
    modulus = c(FALSE, FALSE)
  )
}

#'@export
#'@rdname submodel-specs
fiegarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
    include_mean = TRUE
) {
  fEGarch_spec(
    model_type = "eGARCH",
    orders = orders,
    long_memo = TRUE,
    cond_dist = cond_dist,
    include_mean = include_mean,
    powers = c(1, 1),
    modulus = c(FALSE, FALSE)
  )
}

#'@export
#'@rdname submodel-specs
mloggarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
    include_mean = TRUE
) {
  fEGarch_spec(
    model_type = "eGARCH",
    orders = orders,
    long_memo = FALSE,
    cond_dist = cond_dist,
    include_mean = include_mean,
    powers = c(0, 0),
    modulus = c(TRUE, TRUE)
  )
}

#'@export
#'@rdname submodel-specs
mafiloggarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
    include_mean = TRUE
) {
  fEGarch_spec(
    model_type = "eGARCH",
    orders = orders,
    long_memo = TRUE,
    cond_dist = cond_dist,
    include_mean = include_mean,
    powers = c(0, 0),
    modulus = c(TRUE, TRUE)
  )
}

#'@export
#'@rdname submodel-specs
loggarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
    include_mean = TRUE
) {
  fEGarch_spec(
    model_type = "logGARCH",
    orders = orders,
    long_memo = FALSE,
    cond_dist = cond_dist,
    include_mean = include_mean
  )
}

#'@export
#'@rdname submodel-specs
filoggarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
    include_mean = TRUE
) {
  fEGarch_spec(
    model_type = "logGARCH",
    orders = orders,
    long_memo = TRUE,
    cond_dist = cond_dist,
    include_mean = include_mean
  )
}

sigt_egarch_short_R <- function(theta, x, dfun1, p, q,
              incl_mean, extra_par, skew_par,
              int_fun_asy, int_fun_mag, pows, mods,
              mode) {

  imean <- !incl_mean

  mu <- c(0, theta[[1]])[[incl_mean + 1]]
  step0 <- 2 - imean
  omega_sig <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q - 1

  phi <- theta[(step0 + 1):step1]

  if (q > 1) {
    psi <- c(1, theta[(step1 + 1):step2])
  } else {
    psi <- 1
  }
  kappa <- theta[[(step2 + 1)]]
  gamma <- theta[[step2 + 2]]
  if (extra_par) {
    shape <- theta[[step2 + 3]]
  } else {
    shape <- 0
  }
  if (skew_par) {
    skew <- theta[[step2 + 3 + extra_par]]
  } else {
    skew <- 0
  }

  alpha <- psi * kappa
  beta <- psi * gamma

  E_asy <- integrate(f = int_fun_asy, lower = -Inf, upper = Inf, shape = shape, skew = skew,
                  dfun = dfun1, pow = pows[[1]],
                  subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                  stop.on.error = FALSE)[[1]]
  E_mag <- integrate(f = int_fun_mag, lower = -Inf, upper = Inf, shape = shape, skew = skew,
                  dfun = dfun1, pow = pows[[2]],
                  subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                  stop.on.error = FALSE)[[1]]

  Elnsig2 <- omega_sig
  omega <- Elnsig2 * (1 - sum(phi))

  sigt <- c(sigt_egarch_shortCpp(
    x = x,
    mu = mu,
    omega = omega,
    phi = phi,
    alpha = alpha,
    beta = beta,
    E_asy = E_asy,
    E_mag = E_mag,
    Elnsig2,
    pows, mods, mode
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape)

}

sigt_egarch_long_R <- function(theta, x, dfun1, p, q,
              incl_mean, extra_par, skew_par,
              int_fun_asy, int_fun_mag, pows, mods,
              mode) {


  imean <- !incl_mean

  mu <- c(0, theta[[1]])[[incl_mean + 1]]
  step0 <- 2 - imean
  Elnsig2 <- theta[[step0]]   # The expectation of ln(sig_t^2)
  step1 <- step0 + p
  step2 <- step1 + q - 1

  phi <- theta[(step0 + 1):step1]

  if (q > 1) {
    psi <- c(1, theta[(step1 + 1):step2])
    psi_2 <- psi[-1]
  } else {
    psi <- 1
    psi_2 <- numeric(0)
  }
  kappa <- theta[[(step2 + 1)]]
  gamma <- theta[[step2 + 2]]
  d <- theta[[step2 + 3]]
  if (extra_par) {
    shape <- theta[[step2 + 4]]
  } else {
    shape <- 0
  }
  if (skew_par) {
    skew <- theta[[step2 + 4 + extra_par]]
  } else {
    skew <- 0
  }

  n <- length(x)
  coef_inf <- esemifar::farima_to_ma(ar = phi, ma = psi_2, d = d,
                                     max_i = n - 1)

  E_asy <- integrate(f = int_fun_asy, lower = -Inf, upper = Inf, shape = shape, skew = skew,
                  dfun = dfun1, pow = pows[[1]],
                  subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                  stop.on.error = FALSE)[[1]]
  E_mag <- integrate(f = int_fun_mag, lower = -Inf, upper = Inf, shape = shape, skew = skew,
                  dfun = dfun1, pow = pows[[2]],
                  subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                  stop.on.error = FALSE)[[1]]

  sigt <- c(sigt_egarch_longCpp(
    x = x,
    mu = mu,
    coef_inf = coef_inf,
    kappa = kappa,
    gamma = gamma,
    E_asy = E_asy,
    E_mag = E_mag,
    Elnsig2 = Elnsig2,
    powers = pows,
    modulus = mods,
    mode = mode
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape)

}

goal_fun_creator_egarch <- function(
    theta,
    x,
    dfun1,
    dfun2,
    p,
    q,
    incl_mean,
    extra_par,
    skew_par,
    int_fun_asy,
    int_fun_mag,
    sigt_fun,
    pows, mods,
    mode
) {

  res_list <- sigt_fun(theta, x, dfun1, p, q,
              incl_mean, extra_par, skew_par,
              int_fun_asy, int_fun_mag, pows, mods,
              mode)

  sigt <- res_list$sigt
  mu <- res_list$mu
  skew <- res_list$skew
  shape <- res_list$shape

  out <- -sum(log(dfun2(x, mu, sigt, shape, skew)))

  out

}



setMethod("fEGarch_fit", "egarch_type_spec",
  function(spec = egarch_spec(), rt, start_pars = NULL, LB = NULL,
    UB = NULL) {

  rt_core <- zoo::coredata(rt)

  order <- orders(spec)
  p <- order[[1]]
  q <- order[[2]]
  cond_d <- cond_dist(spec)

  dfun1 <- lookup_table[["fun1"]][[cond_d]]
  dfun2 <- lookup_table[["fun2"]][[cond_d]]

  skew_par <- lookup_table[["skew_par"]][[cond_d]]
  extra_par <- lookup_table[["extra_par"]][[cond_d]]

  lm <- long_memo(spec)
  sigt_fun <- sigt_funs_egarch[[lm + 1]]

  mod <- modulus(spec)
  mod_asy <- mod[[1]]
  mod_mag <- mod[[2]]

  pow <- powers(spec)
  pow_asy <- pow[[1]]
  pow_mag <- pow[[2]]
  pow_asy_check <- (pow_asy != 0) + 1
  pow_mag_check <- (pow_mag != 0) + 1

  int_fun_asy <- exp_funs[["asy"]][[mod_asy + 1]][[pow_asy_check]]
  int_fun_mag <- exp_funs[["mag"]][[mod_mag + 1]][[pow_mag_check]]

  mode <- mode_mat[pow_asy_check, pow_mag_check]

  incl_mean <- include_mean(spec)

  goal_fun <- function(theta, rt) {
    goal_fun_creator_egarch(
      theta = theta,
      x = rt,
      dfun1 = dfun1,
      dfun2 = dfun2,
      p = p,
      q = q,
      incl_mean = incl_mean,
      extra_par = extra_par,
      skew_par = skew_par,
      int_fun_asy = int_fun_asy,
      int_fun_mag = int_fun_mag,
      sigt_fun = sigt_fun,
      pows = pow,
      mods = mod,
      mode = mode
    )
  }

  vals <- lookup_table[["vals"]][[cond_d]]

  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  if (lm) {
    d_spar <- 0.1
    d_low <- 0
    d_up <- 0.4999
    d_name <- "d"
  } else {
    d_spar <- d_low <- d_up <- numeric(0)
    d_name <- character(0)
  }

  if (is.null(start_pars)) {

    start_pars <- c(
      list(NULL, 0)[[incl_mean + 1]],
      -9,
      rep(0, p),
      rep(0, q - 1),
      0, 0,
      d_spar,
      sval,
      list(NULL, 1)[[skew_par + 1]]
    )

  }

  if (is.null(LB)) {

    LB <- c(
      list(NULL, -2)[[incl_mean + 1]],
      -15,
      rep(-2, p),
      rep(-2, q - 1),
      rep(-2, 2),
      d_low,
      lval,
      list(NULL, 1e-15)[[skew_par + 1]]
    )

  }

  if (is.null(UB)) {

    UB <- c(
      list(NULL, 2)[[incl_mean + 1]],
      -3,
      rep(2, p),
      rep(2, q - 1),
      rep(2, 2),
      d_up,
      uval,
      list(NULL, Inf)[[skew_par + 1]]
    )

  }

  result <- suppressWarnings(Rsolnp::solnp(
    pars = start_pars,
    fun = goal_fun,
    LB = LB,
    UB = UB,
    control = list(trace = FALSE),
    rt = rt_core
  ))

  pars <- result$pars
  llhood <- -tail(result$values, 1)

  sigt <- sigt_fun(unname(pars), rt_core, dfun1, p, q,
              incl_mean, extra_par, skew_par,
              int_fun_asy, int_fun_mag, pow, mod,
              mode)$sigt

  mu <- c(0, pars[[1]])[[incl_mean + 1]]
  etat <- (rt - mu) / sigt

  if (q > 1) {
    psi_names <- paste0("psi", 2:q)
  } else {
    psi_names <- NULL
  }

  par_names <- c(
    list(NULL, "mu")[[incl_mean + 1]],
    "omega_sig",
    paste0("phi", 1:p),
    psi_names,
    "kappa",
    "gamma",
    d_name,
    list(NULL, c("shape", "df")[[(cond_d %in% c("std", "sstd")) + 1]])[[extra_par + 1]],
    list(NULL, "skew")[[skew_par + 1]]
  )

  serrors <- sqrt(diag(solve(result$hessian)))
  names(pars) <- par_names
  names(serrors) <- par_names

  k <- length(pars)
  n <- length(rt)

  m2llhood <- -2 * llhood
  bic <- (m2llhood + k * log(n)) / n
  aic <- (m2llhood + 2 * k) / n

  if (inherits(rt, what = "zoo")) {
    tp <- as.Date(time(rt))
    sigt <- zoo::zoo(sigt, order.by = tp)
    etat <- zoo::zoo(etat, order.by = tp)
  }

  list_out <- fEGarch_fit_egarch(
    pars = pars,
    se = serrors,
    rt = rt,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = c("aic" = aic, "bic" = bic),
    cond_dist = cond_d,
    orders = order,
    powers = pow,
    modulus = mod,
    long_memo = lm
  )

  list_out

})

sigt_loggarch_short_R <- function(theta, x, dfun1, p, q,
              incl_mean, extra_par, skew_par,
              int_fun) {

  imean <- !incl_mean
  mu <- c(0, theta[[1]])[[incl_mean + 1]]
  step0 <- 2 - imean
  omega_sig <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q
  phi <- theta[(step0 + 1):step1]
  psi <- theta[(step1 + 1):step2]

  if (extra_par) {
    shape <- theta[[step2 + 1]]
  } else {
    shape <- 0
  }
  if (skew_par) {
    skew <- theta[[step2 + 1 + extra_par]]
  } else {
    skew <- 0
  }

  Elneta2 <- integrate(f = int_fun, lower = -Inf, upper = Inf, shape = shape, skew = skew,
                  dfun = dfun1,
                  subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                  stop.on.error = FALSE)[[1]]

  Elnsig2 <- omega_sig
  omega <- omega_sig * (1 - sum(phi))
  l <- max(p, q)
  phi_s <- psi_s <- rep(0, l)
  phi_s[1:p] <- phi
  psi_s[1:q] <- psi
  alpha <- phi_s + psi_s

  sigt <- c(sigt_loggarch_short(
    x = x,
    mu = mu,
    omega = omega,
    phi = phi,
    psi = alpha,
    Elneta2 = Elneta2,
    Elnsig2 = Elnsig2
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape)

}

sigt_loggarch_long_R <- function(theta, x, dfun1, p, q,
              incl_mean, extra_par, skew_par,
              int_fun) {


  imean <- !incl_mean
  mu <- c(0, theta[[1]])[[incl_mean + 1]]
  step0 <- 2 - imean
  omega_sig <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q
  phi <- theta[(step0 + 1):step1]
  psi <- theta[(step1 + 1):step2]
  d <- theta[[step2 + 1]]
  if (extra_par) {
    shape <- theta[[step2 + 2]]
  } else {
    shape <- 0
  }
  if (skew_par) {
    skew <- theta[[step2 + 2 + extra_par]]
  } else {
    skew <- 0
  }

  Elneta2 <- integrate(f = int_fun, lower = -Inf, upper = Inf, shape = shape, skew = skew,
                  dfun = dfun1,
                  subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                  stop.on.error = FALSE)[[1]]

  n <- length(x)
  coef_inf <- esemifar::farima_to_ma(ar = phi, ma = psi, d = d,
                                     max_i = n)[-1]

  Elnsig2 <- omega_sig

  sigt <- c(sigt_loggarch_long(
    x = x,
    mu = mu,
    coef_inf = coef_inf,
    Elneta2 = Elneta2,
    Elnsig2 = Elnsig2
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape)

}

goal_fun_creator_loggarch <- function(
    theta,
    x,
    dfun1,
    dfun2,
    p,
    q,
    incl_mean,
    extra_par,
    skew_par,
    int_fun,
    sigt_fun
) {

  res_list <- sigt_fun(theta, x, dfun1, p, q,
              incl_mean, extra_par, skew_par,
              int_fun)

  sigt <- res_list$sigt
  mu <- res_list$mu
  skew <- res_list$skew
  shape <- res_list$shape

  out <- -sum(log(dfun2(x, mu, sigt, shape, skew)))

  out

}

setMethod("fEGarch_fit", "loggarch_type_spec",
  function(spec = loggarch_spec(), rt, start_pars = NULL, LB = NULL, UB = NULL) {
  rt_core <- zoo::coredata(rt)

  order <- orders(spec)
  p <- order[[1]]
  q <- order[[2]]
  cond_d <- cond_dist(spec)

  dfun1 <- lookup_table[["fun1"]][[cond_d]]
  dfun2 <- lookup_table[["fun2"]][[cond_d]]

  skew_par <- lookup_table[["skew_par"]][[cond_d]]
  extra_par <- lookup_table[["extra_par"]][[cond_d]]

  lm <- long_memo(spec)
  sigt_fun <- sigt_funs_loggarch[[lm + 1]]

  int_fun <- function(x, shape, skew, dfun) {
    log(x^2) * dfun(x, shape = shape, skew = skew)
  }

  incl_mean <- include_mean(spec)

  goal_fun <- function(theta, rt) {
    goal_fun_creator_loggarch(
      theta = theta,
      x = rt,
      dfun1 = dfun1,
      dfun2 = dfun2,
      p = p,
      q = q,
      incl_mean = incl_mean,
      extra_par = extra_par,
      skew_par = skew_par,
      int_fun = int_fun,
      sigt_fun = sigt_fun
    )
  }

  vals <- lookup_table[["vals"]][[cond_d]]

  sval <- vals[[1]]
  lval <- vals[[2]]
  uval <- vals[[3]]

  if (lm) {
    d_spar <- 0.1
    d_low <- 0
    d_up <- 0.4999
    d_name <- "d"
  } else {
    d_spar <- d_low <- d_up <- numeric(0)
    d_name <- character(0)
  }

  if (is.null(start_pars)) {

    start_pars <- c(
      list(NULL, 0.0001)[[incl_mean + 1]],
      -9,
      rep(0.2, p),
      rep(-0.1, q),
      d_spar,
      sval,
      list(NULL, 1)[[skew_par + 1]]
    )

  }

  if (is.null(LB)) {

    LB <- c(
      list(NULL, -0.1)[[incl_mean + 1]],
      -15,
      rep(-2, p),
      rep(-2, q),
      d_low,
      lval,
      list(NULL, 1e-15)[[skew_par + 1]]
    )

  }

  if (is.null(UB)) {

    UB <- c(
      list(NULL, 0.1)[[incl_mean + 1]],
      -3,
      rep(2, p),
      rep(2, q),
      d_up,
      uval,
      list(NULL, Inf)[[skew_par + 1]]
    )

  }

  result <- suppressWarnings(Rsolnp::solnp(
    pars = start_pars,
    fun = goal_fun,
    LB = LB,
    UB = UB,
    control = list(trace = FALSE),
    rt = rt_core
  ))

  pars <- result$pars
  llhood <- -tail(result$values, 1)

  sigt <- sigt_fun(unname(pars), rt_core, dfun1, p, q,
              incl_mean, extra_par, skew_par,
              int_fun)$sigt

  mu <- c(0, pars[[1]])[[incl_mean + 1]]
  etat <- (rt - mu) / sigt

  par_names <- c(
    list(NULL, "mu")[[incl_mean + 1]],
    "omega_sig",
    paste0("phi", 1:p),
    paste0("psi", 1:q),
    d_name,
    list(NULL, c("shape", "df")[[(cond_d %in% c("std", "sstd")) + 1]])[[extra_par + 1]],
    list(NULL, "skew")[[skew_par + 1]]
  )

  serrors <- sqrt(diag(solve(result$hessian)))
  names(pars) <- par_names
  names(serrors) <- par_names

  k <- length(pars)
  n <- length(rt)

  m2llhood <- -2 * llhood
  bic <- (m2llhood + k * log(n)) / n
  aic <- (m2llhood + 2 * k) / n

  if (inherits(rt, what = "zoo")) {
    tp <- as.Date(time(rt))
    sigt <- zoo::zoo(sigt, order.by = tp)
    etat <- zoo::zoo(etat, order.by = tp)
  }

  list_out <- fEGarch_fit_loggarch(
    pars = pars,
    se = serrors,
    rt = rt,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = c("aic" = aic, "bic" = bic),
    cond_dist = cond_d,
    orders = order,
    long_memo = lm
  )

  list_out
  })

#'@export
setMethod(
  "show",
  "fEGarch_fit_egarch",
  function(object) {
  x <- object
  par_names <- names(x@pars)
  se <- unname(x@se)
  pars <- unname(x@pars)

  tvals <- pars / se
  atvals <- abs(tvals)

  pvals <- 2 * (1 - pnorm(atvals))

  df <- data.frame(
    par = sprintf("%.4f", pars),
    se = sprintf("%.4f", se),
    tval = sprintf("%.4f", tvals),
    pval = sprintf("%.4f", pvals)
  )
  row.names(df) <- par_names

  part1 <- paste0(
    "*************************************\n",
    "*     Fitted EGARCH Family Model    *\n",
    "*************************************\n",
    " \n",
    "Type: eGARCH\n",
    "Orders: (", x@orders[[1]], ", ", x@orders[[2]], ")\n",
    "Modulus: (", x@modulus[[1]], ", ", x@modulus[[2]], ")\n",
    "Powers: (", x@powers[[1]], ", ", x@powers[[2]], ")\n",
    "Long memory: ", x@long_memo, "\n",
    "Cond. distribution: ", x@cond_dist, "\n",
    " \n",
    "Fitted parameters:\n"
  )

  part2 <- paste0(
    " \n",
    "Information criteria:\n",
    "AIC: ", sprintf("%.4f", unname(x@inf_criteria[[1]])),
    ", BIC: ", sprintf("%.4f", unname(x@inf_criteria[[2]]))
  )

  cat(part1, fill = TRUE)
  print(df)
  cat(part2, fill = TRUE)
  }
)

#'@export
setMethod(
  "show",
  "fEGarch_fit_loggarch",
  function(object) {
  x <- object
  par_names <- names(x@pars)
  se <- unname(x@se)
  pars <- unname(x@pars)

  tvals <- pars / se
  atvals <- abs(tvals)

  pvals <- 2 * (1 - pnorm(atvals))

  df <- data.frame(
    par = sprintf("%.4f", pars),
    se = sprintf("%.4f", se),
    tval = sprintf("%.4f", tvals),
    pval = sprintf("%.4f", pvals)
  )
  row.names(df) <- par_names

  part1 <- paste0(
    "*************************************\n",
    "*     Fitted EGARCH Family Model    *\n",
    "*************************************\n",
    " \n",
    "Type: logGARCH\n",
    "Orders: (", x@orders[[1]], ", ", x@orders[[2]], ")\n",
    "Long memory: ", x@long_memo, "\n",
    "Cond. distribution: ", x@cond_dist, "\n",
    " \n",
    "Fitted parameters:\n"
  )

  part2 <- paste0(
    " \n",
    "Information criteria:\n",
    "AIC: ", sprintf("%.4f", unname(x@inf_criteria[[1]])),
    ", BIC: ", sprintf("%.4f", unname(x@inf_criteria[[2]]))
  )

  cat(part1, fill = TRUE)
  print(df)
  cat(part2, fill = TRUE)
  }
)

#'@export
setMethod("sigt", "fEGarch_fit", function(x) {x@sigt})
#'@export
setMethod("etat", "fEGarch_fit", function(x) {x@etat})
#'@export
setMethod("inf_criteria", "fEGarch_fit", function(x) {x@inf_criteria})
#'@export
setMethod("llhood", "fEGarch_fit", function(x) {x@llhood})
#'@export
setMethod("pars", "fEGarch_fit", function(x) {x@pars})
#'@export
setMethod("se", "fEGarch_fit", function(x) {x@se})
