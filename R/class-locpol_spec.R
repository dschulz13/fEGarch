setClass("locpol_spec",
  slots = c(
    poly_order = "numeric",
    kernel_order = "numeric",
    boundary_method = "character"
  ),
  prototype = list(
    poly_order = 3,
    kernel_order = 1,
    boundary_method = "extend"
  )
)

setValidity("locpol_spec", function(object){
  if (poly_checker_locpol(object@poly_order)) {
    "@poly_order must be of length one and either 1 or 3"
  } else if (kernel_checker_locpol(object@kernel_order)) {
    "@kernel_order must be of length one and one of the integers from 0, 1, 2, 3"
  } else if (bmethod_checker_locpol(object@boundary_method)) {
    '@boundary_method must be of length one and either "shorten" or "extend"'
  } else {
    TRUE
  }
})

#'Specification of Nonparametric Local Polynomial Models
#'
#'Specify the nonparametric local polynomial model part in a
#'semiparametric volatility model.
#'
#'@param poly_order a single numeric value, in detail either \code{1}
#'or \code{3}, that represents the local polynomial order.
#'@param kernel_order a single numeric value representing
#'the smoothness of the underlying kernel function;
#'available are \code{0} (uniform), \code{1} (epanechnikov),
#'\code{2} (bisquare), and \code{3} (triweight).
#'@param boundary_method a single character value indicating
#'the smoothing concept to use at boundary points; for
#'\code{"extend"}, the smoothing window is extended toward the
#'interior by the amount that is lost toward the boundary; for
#'\code{"shorten"}, there is no compensation in the smoothing window
#'toward the interior for the loss of window width toward the boundary,
#'i.e. the total smoothing window width reduces more and more when
#'getting closer to the first and last time points.
#'
#'@return
#'An object of class \code{"locpol_spec"} is returned.
#'
#'@details
#'Test
#'
#'@export
#'
#'@examples
#'locpol_spec()
#'locpol_spec(poly_order = 1)
#'locpol_spec(kernel_order = 2)
#'
locpol_spec <- function(
    poly_order = c(3, 1),
    kernel_order = c(1, 0, 2, 3),
    boundary_method = c("extend", "shorten")
) {

  poly_order <- match.arg(poly_order)
  kernel_order <- match.arg(kernel_order)
  boundary_method <- match.arg(boundary_method)

  new("locpol_spec",
    poly_order = poly_order,
    kernel_order = kernel_order,
    boundary_method = boundary_method
  )

}

#'Accessors for Class \code{"locpol_spec"}
#'
#'Access and change elements in objects of class \code{"locpol_spec"}.
#'The method names represent the name of the element to access /
#'manipulate.
#'
#'@param x the input object or object to modify.
#'@param value the value to modify the object \code{x} with.
#'
#'@details
#'These methods are intended to be used for accessing or manipulating
#'individual elements of objects of class \code{"locpol_spec"}.
#'
#'@name locpol_spec_methods
#'@aliases poly_order,locpol_spec-method
#'
#'@return
#'These methods return an object of class \code{"locpol_spec"}.
#'
#'@export
#'
#'@examples
#'test_obj <- locpol_spec()
#'poly_order(test_obj)
#'poly_order(test_obj) <- 1
#'poly_order(test_obj)
#'
setMethod("poly_order", "locpol_spec", function(x) {x@poly_order})
#'@export
#'@rdname locpol_spec_methods
#'@aliases kernel_order,locpol_spec-method
setMethod("kernel_order", "locpol_spec", function(x) {x@kernel_order})
#'@export
#'@rdname locpol_spec_methods
#'@aliases boundary_method,locpol_spec-method
setMethod("boundary_method", "locpol_spec", function(x) {x@boundary_method})
#'@export
#'@rdname locpol_spec_methods
setMethod("poly_order<-", "mean_spec", function(x, value) {x@poly_order <- value; validObject(x); x})
#'@export
#'@rdname locpol_spec_methods
setMethod("kernel_order<-", "mean_spec", function(x, value) {x@kernel_order <- value; validObject(x); x})
#'@export
#'@rdname locpol_spec_methods
setMethod("boundary_method<-", "mean_spec", function(x, value) {x@boundary_method <- value; validObject(x); x})
