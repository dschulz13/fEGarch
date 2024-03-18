#'Simulate from a Standard Normal Distribution
#'
#'@export
#'

rnorm_st <- function(n) {
  rnorm(n, mean = 0, sd = 1)
}

#'Simulate from a Standardized Student's t Distribution
#'
#'@export
#'

rstd_st <- function(n, df) {
  rt(n, df = df) * sqrt((df - 2) / df)
}

#'Simulate from a Standardized Generalized Error Distribution
#'
#'@export
#'

rged_st <- function(n, shape) {
  # Work on this
}



