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

setGeneric("fEGarch_fit", function(spec, rt, start_pars = NULL, LB = NULL, UB = NULL) {standardGeneric("fEGarch_fit")})

#====================================================#

#'@export
setGeneric("sigt", function(x) {standardGeneric("sigt")})
#'@export
setGeneric("etat", function(x) {standardGeneric("etat")})
#'@export
setGeneric("llhood", function(x) {standardGeneric("llhood")})
#'@export
setGeneric("inf_criteria", function(x) {standardGeneric("inf_criteria")})
#'@export
setGeneric("pars", function(x) {standardGeneric("pars")})
#'@export
setGeneric("se", function(x) {standardGeneric("se")})
#'@export
setGeneric("vcov_mat", function(x) {standardGeneric("vcov_mat")})
