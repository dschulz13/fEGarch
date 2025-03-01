# hessCalc <- function(pars, func, control.args = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
# 				hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2), rt) {
#
#   numDeriv::hessian(func = func, x = pars, method = "Richardson",
#                     method.args = control.args, rt = rt)
#
# }

hessCalc <- function(pars, func, control.args = list(
				eps = 1e-4,
				d = 0.001,
				zero.tol = sqrt(.Machine$double.eps / 7e-7),
				r = 4,
				v = 2,
				show.details = TRUE), rt) {

  numDeriv::hessian(func = func, x = pars, method = "Richardson",
                    method.args = control.args, rt = rt)

}
