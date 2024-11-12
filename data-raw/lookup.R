lookup_table <- list(
  fun1 = list(
    "norm" = function(x, shape, skew) {fEGarch:::pdf_norm_v1(x)},
    "std" = function(x, shape, skew) {fEGarch:::pdf_std_v1(x, df = shape)},
    "ged" = function(x, shape, skew) {fEGarch:::pdf_ged_v1(x, shape = shape)},
    "snorm" = function(x, shape, skew) {fEGarch:::pdf_skew_norm_s(x, skew = skew)},
    "sstd" = function(x, shape, skew) {fEGarch:::pdf_skew_sstd_s(x, df = shape, skew = skew)},
    "sged" = function(x, shape, skew) {fEGarch:::pdf_skew_sged_s(x, shape = shape, skew = skew)}
  ),
  fun2 = list(
    "norm" = function(x, mu, sigt, shape, skew) {
      fEGarch:::pdf_norm(x, mu = mu, sigt = sigt)
    },
    "std" = function(x, mu, sigt, shape, skew) {
      fEGarch:::pdf_std(x, mu = mu, sigt = sigt, df = shape)
    },
    "ged" = function(x, mu, sigt, shape, skew) {
      fEGarch:::pdf_ged(x, mu = mu, sigt = sigt, shape = shape)
    },
    "snorm" = function(x, mu, sigt, shape, skew) {
      fEGarch:::pdf_skew_norm_final(x, mu = mu, sigt = sigt, skew = skew)
    },
    "sstd" = function(x, mu, sigt, shape, skew) {
      fEGarch:::pdf_skew_sstd_final(x, mu = mu, sigt = sigt, df = shape, skew = skew)
    },
    "sged" = function(x, mu, sigt, shape, skew) {
      fEGarch:::pdf_skew_sged_final(x, mu = mu, sigt = sigt, shape = shape, skew = skew)
    }
  ),
  # starting_par, lower bound, upper bound
  vals = list(
    "norm" = list(NULL, NULL, NULL),
    "std" = list(10, 4, Inf),
    "ged" = list(2, 1e-15, Inf),
    "snorm" = list(NULL, NULL, NULL),
    "sstd" = list(10, 4, Inf),
    "sged" = list(2, 1e-15, Inf)
  ),
  extra_par = list(
    "norm" = FALSE,
    "std" = TRUE,
    "ged" = TRUE,
    "snorm" = FALSE,
    "sstd" = TRUE,
    "sged" = TRUE
  ),
  skew_par = list(
    "norm" = FALSE,
    "std" = FALSE,
    "ged" = FALSE,
    "snorm" = TRUE,
    "sstd" = TRUE,
    "sged" = TRUE
  )
)

fitFUNS <- list(
  "mloggarch" = fEGarch:::mloggarchFit,
  "megarch" = fEGarch:::megarchFit,
  "fiegarch" = fEGarch:::fiegarchFit,
  "filoggarch" = fEGarch:::filoggarchFit,
  "mafiloggarch" = fEGarch:::mafiloggarchFit
)

#usethis::use_data(lookup_table, fitFUNS, internal = TRUE, overwrite = TRUE)
