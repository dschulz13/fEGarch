#'Generics for Model Specification Accessors
#'
#'A collection of generics for accessors of the model specification
#'objects.
#'
#'@param x the input object or object to modify.
#'@param value the value to modify the object \code{x} with.
#'
#'@details
#'These generics are intended to provide a basis to construct methods
#'for adjusting the the output of model specification objects like
#'\code{egarch_spec()}.
#'
#'@name spec_generics
#'
#'@return
#'These generics themselves do not return anything and are just the
#'foundation for more sophisticated methods.
#'
#'@export
#'@rdname spec_generics
setGeneric("orders", function(x) {standardGeneric("orders")})
#'@export
#'@rdname spec_generics
setGeneric("powers", function(x) {standardGeneric("powers")})
#'@export
#'@rdname spec_generics
setGeneric("long_memo", function(x) {standardGeneric("long_memo")})
#'@export
#'@rdname spec_generics
setGeneric("modulus", function(x) {standardGeneric("modulus")})
#'@export
#'@rdname spec_generics
setGeneric("cond_dist", function(x) {standardGeneric("cond_dist")})
#'@export
#'@rdname spec_generics
setGeneric("include_mean", function(x) {standardGeneric("include_mean")})

#'@export
#'@rdname spec_generics
setGeneric("orders<-", function(x, value) {standardGeneric("orders<-")})
#'@export
#'@rdname spec_generics
setGeneric("powers<-", function(x, value) {standardGeneric("powers<-")})
#'@export
#'@rdname spec_generics
setGeneric("long_memo<-", function(x, value) {standardGeneric("long_memo<-")})
#'@export
#'@rdname spec_generics
setGeneric("modulus<-", function(x, value) {standardGeneric("modulus<-")})
#'@export
#'@rdname spec_generics
setGeneric("cond_dist<-", function(x, value) {standardGeneric("cond_dist<-")})
#'@export
#'@rdname spec_generics
setGeneric("include_mean<-", function(x, value) {standardGeneric("include_mean<-")})


#====================================================#

#'An EGARCH Family Fitting Generic
#'
#'A base generic to build further EGARCH fitting methods.
#'It is in itself without use.
#'
#'@param spec the model specification object.
#'@param rt the (log-)return series.
#'@param start_pars a vector of the length of the output parameter
#'vector with suitable, corresponding starting values for the
#'optimization procedure.
#'@param LB a vector of the length of the output parameter
#'vector with suitable, corresponding lower boundary values for the
#'optimization procedure.
#'@param UB a vector of the length of the output parameter
#'vector with suitable, corresponding upper boundary values for the
#'optimization procedure.
#'
#'@details
#'This generic itself is without use.
#'
#'
#'
setGeneric("fEGarch_fit", function(spec, rt, start_pars = NULL, LB = NULL, UB = NULL) {standardGeneric("fEGarch_fit")})
