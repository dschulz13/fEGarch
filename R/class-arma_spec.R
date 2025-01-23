setClass("arma_spec",
  slots = c(
    orders = "numeric",
    long_memo = "logical"
  ),
  prototype = list(
    orders = c(0, 0),
    long_memo = FALSE
  )
)

setValidity("arma_spec", function(object){
  if (order_checker_arma(object@orders)) {
    "@orders must be of length two and both elements must be at least zero"
  } else if (lmemo_checker(object@long_memo)) {
    "@long_memo must be of length one"
  } else {
    TRUE
  }
})

#'Specification of ARMA Models
#'
#'Specify the model for the conditional mean in a dual model, where the
#'conditional mean is modelled through an ARMA or a FARIMA model and the
#'conditional standard deviations through a GARCH-type model simultaneously.
#'
#'@param orders a two-element numeric vector with the model
#'orders; the first element is the autoregressive order \eqn{p^{*}}, while
#'the second element is the moving-average order \eqn{q^{*}}.
#'@param long_memo a logical value that indicates whether the
#'long-memory version of the model should be considered or not.
#'
#'@return
#'An object of class \code{"arma_spec"} is returned.
#'
#'@details
#'Let \eqn{\left\{y_t\right\}}, with \eqn{t \in \mathbb{Z}} as the time index,
#'be a theoretical time series that follows
#'\deqn{\beta(B)(1- B)^{D}(y_t - \mu)=\alpha(B)r_t,}
#'where \eqn{\beta(B) = 1 - \sum_{i=1}^{p^{*}}\beta_i B^{i}} and
#'\eqn{\alpha(B) = 1 + \sum_{j=1}^{q^{*}}\alpha_j B^{j}} are the AR- and MA-polynomials
#'of orders \eqn{p^{*}} and \eqn{q^{*}}, respectively, with real coefficients
#'\eqn{\beta_i}, \eqn{i=1,\dots,p^{*}}, and \eqn{\alpha_j}, \eqn{j=1,\dots,q^{*}}.
#'\eqn{B} is the backshift operator. \eqn{\beta(B)} and \eqn{\alpha(B)} are
#'commonly assumed to be without common roots and to have roots outside of
#'the unit circle.
#'Furthermore, \mu is a real-valued coefficient representing the unconditional
#'mean in \eqn{\left\{y_t\right\}}. \eqn{D \in [0, 0.5)} is the fractional
#'differencing parameter. \eqn{\left\{r_t\right\}} is a zero-mean (weak) white
#'noise process, for example a member of the GARCH-models (with mean set to zero)
#'presented in this package (see the
#'descriptions in \link{fEGarch_spec}, \link{fiaparch}, \link{lmgarch}).
#'
#'The for \eqn{D=0}, which can be achieved through \code{long_memo = FALSE},
#'the formulas above describe an autoregressive moving-average (ARMA) model.
#'For \eqn{D \in (0, 0.5)}, they describe a fractionally integrated ARMA (FARIMA)
#'model.
#'
#'@export
#'
#'@examples
#'arma_spec()
#'arma_spec(orders = c(1, 1))
#'
arma_spec <- function(
    orders = c(0, 0),
    long_memo = FALSE
) {

  new("arma_spec",
    orders = orders,
    long_memo = long_memo
  )

}

#'Accessors for Class \code{"arma_spec"}
#'
#'Access and change elements in objects of class \code{"arma_spec"}.
#'The method names represent the name of the element to access /
#'manipulate.
#'
#'@param x the input object or object to modify.
#'@param value the value to modify the object \code{x} with.
#'
#'@details
#'These methods are intended to be used for accessing or manipulating
#'individual elements of objects of class \code{"arma_spec"}.
#'
#'@name arma_spec_methods
#'@aliases orders,arma_spec-method
#'
#'@return
#'These methods return an object of class \code{"arma_spec"}.
#'
#'@export
#'
#'@examples
#'test_obj <- arma_spec()
#'orders(test_obj)
#'orders(test_obj) <- c(1, 1)
#'orders(test_obj)
#'
setMethod("orders", "arma_spec", function(x) {x@orders})
#'@export
#'@rdname arma_spec_methods
#'@aliases long_memo,arma_spec-method
setMethod("long_memo", "arma_spec", function(x) {x@long_memo})
#'@export
#'@rdname arma_spec_methods
setMethod("orders<-", "arma_spec", function(x, value) {x@orders <- value; validObject(x); x})
#'@export
#'@rdname arma_spec_methods
setMethod("long_memo<-", "arma_spec", function(x, value) {x@long_memo <- value; validObject(x); x})
#'@export
