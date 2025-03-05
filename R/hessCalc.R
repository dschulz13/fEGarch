# Calculation of Hessian matrix using numDeriv::hessian;
# usually produces more stable hessian matrix.
hessCalc <- function(pars, func, control.args = list(
				eps = 1e-4,
				d = 1e-3,
				zero.tol = sqrt(.Machine$double.eps / 7e-7),
				r = 4,
				v = 2,
				show.details = FALSE), rt) {


  out <- numDeriv::hessian(func = func, x = pars, method = "Richardson",
                    method.args = control.args, rt = rt)

  # Try different settings
  while (control.args$eps > 1e-14 && (any(is.nan(out)) || any(diag(solve(out)) < 0))) {

    control.args$eps <- control.args$eps / 10
    out <- numDeriv::hessian(func = func, x = pars, method = "Richardson",
              method.args = control.args, rt = rt)

  }

  control.args$eps <- 1e-4

  # Try different settings
  while (control.args$d > 1e-10 && (any(is.nan(out)) || any(diag(solve(out)) < 0))) {

    control.args$d <- control.args$d / 10
    out <- numDeriv::hessian(func = func, x = pars, method = "Richardson",
              method.args = control.args, rt = rt)

  }

  if (any(is.nan(out))) {
    out <- matrix(NA, nrow = length(pars), ncol = length(pars))
    warning("Unable to compute Hessian matrix. No standard errors available as consequence.")
  }

  out
}
