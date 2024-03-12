### Constant E(|Zt|) for different error distributions

### Standard normal distribution

E_zt_ <- sqrt(2 / pi)

### Standardized t-distribution

E_zt_ <- function(df) {
  2 * sqrt((df - 2) / pi) * (gamma((df + 1) / 2) / (gamma(df / 2) * (df - 1)))
}

# Skew-normal distribution

E_zt_ <- function(skew_par) {

  delta <- skew_par / sqrt(1 + skew_par^2)
  scale <- 1 / sqrt(1 - (2 * delta^2) / pi)
  location <- -scale * delta * sqrt(2 / pi)


}

# Skew-t distribution

E_zt_ <- function(df, skew_par) {

}

# GED

E_zt_ <- function(shape) {
  scale <- sqrt(1 / (gamma(3 / shape) / gamma(1 / shape)))
  scale * gamma(2 / shape) / gamma(1 / shape)
}

# Skew-GED distribution


