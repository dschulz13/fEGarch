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
#'orders.
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
  if (long_memo && any(orders > 1)) {   # Adjust orders for long-memory models
    orders <- c(1, 1)
    message("Only orders p=1 with q=1 are currently supported for the models with fractional differencing by fAParch_spec.")
  }

  new("fAParch_spec",
    model_type = mtype,
    orders = orders,
    long_memo = long_memo,
    cond_dist = cond_d,
    adjust = adjust
  )

}

##########################################

#'Adjusted FIAPARCH Family Submodel Specification

#'Wrappers of \code{\link{fAParch_spec}} that create
#'specifications of specific submodels of the broad
#'adjusted FIAPARCH class.
#'
#'@param orders a two-element numeric vector with the model orders.
#'@param cond_dist a character value stating the underlying
#'conditional distribution to consider; available are a normal
#'distribution (\code{"norm"}), a \eqn{t}-distribution
#'(\code{"std"}), a generalized error distribution
#'(\code{"ged"}), an average Laplace distribution (\code{"ald"})
#'and the skewed versions of them
#'(\code{"snorm"}, \code{"sstd"}, \code{"sged"}, \code{"sald"}).
#'
#'@details
#'\itemize{
#'\item{\code{aparch_spec()} specifies an APARCH model; see \code{\link{aparch}.}}
#'\item{\code{garch_spec()} specifies a GARCH model; see \code{\link{garch}.}}
#'\item{\code{tgarch_spec()} specifies a TGARCH model; see \code{\link{tgarch}.}}
#'\item{\code{gjrgarch_spec()} specifies a GJR-GARCH model; see
#'\code{\link{gjrgarch}.}}
#'\item{\code{ngarch_spec()} specifies an NGARCH model; see \code{\link{ngarch}.}}
#'\item{\code{fiaparch_spec()} specifies a FIAPARCH model; see
#'\code{\link{fiaparch}.}}
#'\item{\code{figarch_spec()} specifies a FIGARCH model; see
#'\code{\link{figarch}.}}
#'\item{\code{fitgarch_spec()} specifies a FITGARCH model; see
#'\code{\link{fitgarch}.}}
#'\item{\code{figjrgarch_spec()} specifies a FIGJR-GARCH model; see
#'\code{\link{figjrgarch}.}}
#'\item{\code{fingarch_spec()} specifies a FINGARCH model; see
#'\code{\link{fingarch}.}}
#'\item{\code{afiaparch_spec()} specifies a AFIAPARCH model; see
#'\code{\link{afiaparch}.}}
#'\item{\code{afitgarch_spec()} specifies a AFITGARCH model; see
#'\code{\link{afitgarch}.}}
#'\item{\code{afigjrgarch_spec()} specifies a AFIGJR-GARCH model; see
#'\code{\link{afigjrgarch}.}}
#'\item{\code{afingarch_spec()} specifies a AFINGARCH model; see
#'\code{\link{afingarch}.}}
#'}
#'
#'These are wrappers for \code{\link{fAParch_spec}} and their results should
#'thus be passed to \code{\link{fAParch}} for model fitting.
#'
#'@name aparch-submodel-specs
#'
#'@export
#'
#'@return
#'An object of class \code{"fAParch_spec"} is returned.
#'
#'@examples
#'spec1 <- afiaparch_spec(orders = c(1, 1), cond_dist = "sged")
#'spec2 <- tgarch_spec(orders = c(2, 1), cond_dist = "std")
#'# The other functions work analogously.
#'

aparch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "aparch", orders = orders, long_memo = FALSE,
               cond_dist = cond_dist, adjust = FALSE)
}

#'@rdname aparch-submodel-specs
#'@export

garch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "garch", orders = orders, long_memo = FALSE,
               cond_dist = cond_dist, adjust = FALSE)
}

#'@rdname aparch-submodel-specs
#'@export

tgarch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "tgarch", orders = orders, long_memo = FALSE,
               cond_dist = cond_dist, adjust = FALSE)
}

#'@rdname aparch-submodel-specs
#'@export

gjrgarch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "gjrgarch", orders = orders, long_memo = FALSE,
               cond_dist = cond_dist, adjust = FALSE)
}

#'@rdname aparch-submodel-specs
#'@export

ngarch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "ngarch", orders = orders, long_memo = FALSE,
               cond_dist = cond_dist, adjust = FALSE)
}

#'@rdname aparch-submodel-specs
#'@export

fiaparch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "aparch", orders = orders, long_memo = TRUE,
               cond_dist = cond_dist, adjust = FALSE)
}

#'@rdname aparch-submodel-specs
#'@export

figarch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "garch", orders = orders, long_memo = TRUE,
               cond_dist = cond_dist, adjust = FALSE)
}

#'@rdname aparch-submodel-specs
#'@export

fitgarch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "tgarch", orders = orders, long_memo = TRUE,
               cond_dist = cond_dist, adjust = FALSE)
}

#'@rdname aparch-submodel-specs
#'@export

figjrgarch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "gjrgarch", orders = orders, long_memo = TRUE,
               cond_dist = cond_dist, adjust = FALSE)
}

#'@rdname aparch-submodel-specs
#'@export

fingarch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "ngarch", orders = orders, long_memo = TRUE,
               cond_dist = cond_dist, adjust = FALSE)
}

#'@rdname aparch-submodel-specs
#'@export

afiaparch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "aparch", orders = orders, long_memo = TRUE,
               cond_dist = cond_dist, adjust = TRUE)
}

#'@rdname aparch-submodel-specs
#'@export

afitgarch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "tgarch", orders = orders, long_memo = TRUE,
               cond_dist = cond_dist, adjust = TRUE)
}

#'@rdname aparch-submodel-specs
#'@export

afigjrgarch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "gjrgarch", orders = orders, long_memo = TRUE,
               cond_dist = cond_dist, adjust = TRUE)
}

#'@rdname aparch-submodel-specs
#'@export

afingarch_spec <- function(orders = c(1, 1), cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")) {
  fAParch_spec(model_type = "ngarch", orders = orders, long_memo = TRUE,
               cond_dist = cond_dist, adjust = TRUE)
}
