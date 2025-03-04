# hessCalc <- function(pars, func, control.args = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
# 				hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2), rt) {
#
#   numDeriv::hessian(func = func, x = pars, method = "Richardson",
#                     method.args = control.args, rt = rt)
#
# }

hessCalc <- function(pars, func, control.args = list(
				eps = 1e-8,
				d = 1e-4,
				zero.tol = sqrt(.Machine$double.eps / 7e-7),
				r = 4,
				v = 2,
				show.details = FALSE), rt) {


  out <- numDeriv::hessian(func = func, x = pars, method = "Richardson",
                    method.args = control.args, rt = rt)

  if (any(is.nan(out)) || any(diag(solve(out)) < 0)) {
    control.args$d <- 1e-5
    out <- numDeriv::hessian(func = func, x = pars, method = "Richardson",
              method.args = control.args, rt = rt)
  }

  if (any(is.nan(out)) || any(diag(solve(out)) < 0)) {
    control.args$d <- 1e-6
    out <- numDeriv::hessian(func = func, x = pars, method = "Richardson",
              method.args = control.args, rt = rt)
  }

  if (any(is.nan(out)) || any(diag(solve(out)) < 0)) {
    control.args$d <- 1e-3
    out <- numDeriv::hessian(func = func, x = pars, method = "Richardson",
              method.args = control.args, rt = rt)
  }

  if (any(is.nan(out)) || any(diag(solve(out)) < 0)) {
    control.args$d <- 1e-2
    out <- numDeriv::hessian(func = func, x = pars, method = "Richardson",
              method.args = control.args, rt = rt)
  }

  if (any(is.nan(out))) {
    out <- matrix(NA, nrow = length(pars), ncol = length(pars))
    warning("Unable to compute Hessian matrix. No standard errors available as consequence.")
  }

  out
}
