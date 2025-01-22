# For class "arma_spec"

order_checker_arma <- function(arma_order) {
  !((length(arma_order) == 2) && (arma_order[[1]] >= 0 && arma_order[[2]] >= 0))
}

lmemo_checker<- function(long_memo) {
  !(length(long_memo) == 1)
}
