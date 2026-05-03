setClass("fAParch_spec",
  contains = "base_garch_spec",
  slots = c(
    model_type = "character",
    adjust = "logical"
  ),
  prototype = list(
    model_type = character(1),
    adjust = logical(1)
  )
)

setValidity("fAParch_spec", function(object){
  if (length(object@adjust) != 1) {
    "@adjust must be of length one"
  } else if (!(all(object@model_type %in% c("garch", "gjrgarch", "tgarch", "ngarch", "aparch")))) {
    "@model_type must be one of the allowed options"
  } else {
    TRUE
  }
})

#'APARCH-Family Model Specification
#'
#'Specify a model from the APARCH class, including short-memory models,
#'unadjusted long-memory models, and adjusted long-memory models.
#'
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
#'(\code{"ged"}), an average Laplace distribution (\code{"ald"})
#'and the skewed versions of them
#'(\code{"snorm"}, \code{"sstd"}, \code{"sged"}, \code{"sald"}).
#'@param model_type a character value indicating the subtype to implement;
#'available are \code{"garch"}, \code{"gjrgarch"}, \code{"tgarch"},
#'\code{"ngarch"} and \code{"aparch"}.
#'@param adjust a logical value indicating, whether the model should be adjusted
#'in the sense of Feng et al. (2026), if \code{long_memo = TRUE}.
#'For \code{model_type = "garch"} or \code{long_memo = FALSE}, \code{adjust}
#'will automatically be set to \code{FALSE}.
#'
#'@details
#'This is a model specification function to then use in conjunction with
#'\code{\link{fAParch}} for model fitting.
#'
#'@return
#'An object of class \code{"fAParch_spec"} is returned.
#'
#'@references
#'\itemize{
#'\item{Feng, Y., Ayensu, O. K., & Schulz, D. (2026). An adjusted FIAPARCH with finite \eqn{2\delta}th moment
#'and long memory in volatility implemented in R. Unpublished manuscript at Paderborn University.}
#'}
#'
#'@export
#'
#'@examples
#'fAParch_spec()
#'
#'fAParch_spec(cond_dist = "std")
#'
#'fAParch_spec(cond_dist = "sged", long_memo = TRUE, model_type = "tgarch")
#'

fAParch_spec <- function(model_type = c("garch", "gjrgarch", "tgarch", "ngarch", "aparch"), orders = c(1, 1), long_memo = FALSE, cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"), adjust = TRUE) {

  mtype <- match.arg(model_type)
  cond_d <- match.arg(cond_dist)

  if (!long_memo || model_type == "garch") adjust <- FALSE

  new("fAParch_spec",
    model_type = mtype,
    orders = orders,
    long_memo = long_memo,
    cond_dist = cond_d,
    adjust = adjust
  )

}


