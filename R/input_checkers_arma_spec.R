# For class "arma_spec"

order_checker_arma <- function(arma_order) {
  !((length(arma_order) == 2) && (arma_order[[1]] >= 0 && arma_order[[2]] >= 0))
}


# Function "lmemo_checker" can be reused from "inpout_checkers_egarch_spec.R".
