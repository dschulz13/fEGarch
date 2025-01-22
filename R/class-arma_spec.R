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
#'Specify the ARMA part in a dual model, where the conditional
#'mean is modelled through an ARMA model and the conditional
#'standard deviations through a GARCH-type model simultaneously.
#'
#'@param orders a two-element numeric vector with the model
#'orders; the first element is the autoregressive order \eqn{p}, while
#'the second element is the moving-average order \eqn{q}.
#'@param long_memo a logical value that indicates whether the
#'long-memory version of the model should be considered or not.
#'
#'
#'@return
#'An object of class \code{"arma_spec"} is returned.
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
