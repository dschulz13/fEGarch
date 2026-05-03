#'Deprecated Functions in Package fEGarch
#'
#'A list of deprecated functions in the fEGarch package.
#'
#'@usage garchm_estim(rt, model = c("garch", "gjrgarch", "tgarch", "ngarch",
#'                                  "aparch", "figarch", "figjrgarch",
#'                                  "fitgarch", "fingarch", "fiaparch",
#'                                  "afigjrgarch", "afitgarch", "afingarch",
#'                                  "afiaparch"),
#'                    orders = c(1, 1), cond_dist = c("norm", "std", "ged",
#'                                                    "ald", "snorm", "sstd",
#'                                                    "sged", "sald"),
#'                    drange = c(0, 1), meanspec = mean_spec(), Drange = c(0, 1),
#'                    nonparspec = locpol_spec(), use_nonpar = FALSE, n_test = 0,
#'                    start_pars = NULL, LB = NULL, UB = NULL, control = list(),
#'                    control_nonpar = list(), mean_after_nonpar = FALSE,
#'                    parallel = TRUE,
#'                    ncores = max(1, future::availableCores() - 1),
#'                    trunc = "none", presample = 50, Prange = c(1, 5),
#'                    skip_vcov = FALSE)
#'
#'@param rt the observed series ordered from past to present; can be
#'a numeric vector, a \code{"zoo"} class time series object, or a
#'\code{"ts"} class time series object.
#'@param model any character object among \code{"garch"}, \code{"gjrgarch"},
#'\code{"aparch"}, \code{"tgarch"}, \code{"ngarch"}, \code{"figarch"}, \code{"figjrgarch"},
#'\code{"fitgarch"}, \code{"fiaparch"}, \code{"fingarch"}, \code{"afigjrgarch"},
#'\code{"afitgarch"}, \code{"afiaparch"} and \code{"afingarch"}.
#'@param orders a two-element numeric vector containing the two model
#'orders \eqn{p} and \eqn{q} (see Details for more information); currently,
#'only the default \code{orders = c(1, 1)} is supported for long-memory
#'models; other specifications
#'of a two-element numeric vector will lead to \code{orders = c(1, 1)} being
#'run and a warning message being returned for long-memory models.
#'@param cond_dist the conditional distribution to consider as a
#'character object; the default is a conditional normal distribution
#'\code{"norm"}; available are also, however, a \eqn{t}-distribution
#'(\code{"std"}), a generalized error distribution (\code{"ged"}), an
#'average Laplace distribution (\code{"ald"}),
#'and their four skewed variants (\code{"snorm"}, \code{"sstd"},
#'\code{"sged"}, \code{"sald"}).
#'@param start_pars the starting parameters for the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used; elements should be set with respect to a series rescaled
#'to have sample variance one.
#'@param drange a two-element numeric vector that gives the boundaries of the
#'search interval for the fractional differencing parameter \eqn{d} in the
#'conditional volatility model part of a long-memory model; is
#'overwritten by the settings of the arguments \code{LB} and \code{UB}.
#'@param meanspec an object of class "mean_spec"; indicates the
#'specifications for the model in the conditional mean.
#'@param Drange a two-element numeric vector that indicates the boundaries
#'of the interval over which to search for the fractional differencing
#'parameter \eqn{D} in a long-memory ARMA-type model in the conditional
#'mean model part; by default,
#'\eqn{D} being searched for on the
#'interval from 0 to \eqn{0.5 - 1\times 10^{-6}}; note that specific
#'settings in the arguments
#'\code{LB} and \code{UB} overwrite this argument.
#'@param nonparspec an object of class \code{"locpol_spec"} returned
#'by \code{\link{locpol_spec}}; defines the settings of the nonparametric
#'smoothing technique for \code{use_nonpar = TRUE}.
#'@param use_nonpar a logical indicating whether or not to implement a
#'semiparametric extension of the volatility model defined through \code{spec};
#'see "Details" for more information.
#'@param n_test a single numerical value indicating, how many observations
#'at the end of \code{rt} not to include in the fitting process and to
#'reserve for backtesting.
#'@param LB the lower boundaries of the parameters in the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used; elements should be set with respect to a series rescaled
#'to have sample variance one.
#'@param UB the upper boundaries of the parameters in the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used; elements should be set with respect to a series rescaled
#'to have sample variance one.
#'@param control a list that is passed to \code{control} of the
#'function \code{solnp} of the package \code{Rsolnp}.
#'@param control_nonpar a list containing changes to the arguments
#'for the hyperparameter estimation algorithm in the nonparametric
#'scale function estimation for
#'\code{use_nonpar = TRUE}; see "Details" for more information.
#'@param mean_after_nonpar only for \code{use_nonpar = TRUE}; considers the unconditional mean
#'of the parametric model part in the QMLE step in a semiparametric model; by default, a zero-mean
#'model is considered for the parametric part in a semiparametric model.
#'@param parallel only relevant for a (skewed) average Laplace (AL)
#'distribution, i.e.
#'if \code{cond_dist} in \code{spec} is set to \code{cond_dist = "ald"} or
#'\code{cond_dist = "sald"}; \code{parallel} is a logical value indicating whether
#'or not the slices for the positive integer-valued parameter of the SM
#'distribution should be fitted in parallel for a potential speed boost; on some
#'laptops, \code{parallel = FALSE} may work faster.
#'@param ncores only relevant for a (skewed) average Laplace (AL)
#'distribution, i.e.
#'if \code{cond_dist} in \code{spec} is set to \code{cond_dist = "ald"} or
#'\code{cond_dist = "sald"}, and if simultaneously \code{parallel = TRUE};
#'\code{ncores} is a single numeric value indicating the number of cores to
#'use for parallel computations.
#'@param trunc a positive integer indicating the finite truncation length of the
#'infinite-order polynomials of the infinite-order representations of the
#'long-memory model parts; the character \code{"none"} is an optional input
#'that specifies that truncation should always be applied back to the first (presample) observation
#'time point, i.e. that maximum length filters should be applied at all times.
#'@param presample the presample length for initialization (for extended EGARCH- / Log-GARCH-type
#'models only relevant for the FARIMA-part, as series in log-transformed
#'conditional variance are initialized by zero).
#'@param Prange a two-element vector that indicates the search boundaries for
#'the parameter \eqn{P} in a (skewed) average Laplace distribution.
#'@param skip_vcov a logical indicating whether or not to skip the computation
#'of the variance-covariance matrix of the parameter estimators and therefore
#'also standard error computation.
#'
#'@details
#'These are functions that may be removed in future versions of the package. In particular
#'users should switch from \code{\link{garchm_estim}} to the \code{\link{fAParch_spec}} / \code{\link{fAParch}}
#'pipeline. Both give the same results, but the new pipeline is more in alignment with
#'the overall workflow idea behind the package.
#'
#'
#'@name fEGarch-deprecated
#'
#'@rdname fEGarch-deprecated
#'
NULL
