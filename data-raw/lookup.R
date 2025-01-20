lookup_table <- list(
  fun1 = list(
    "norm" = function(x, shape, skew) {c(get("pdf_norm_v1", envir = asNamespace("fEGarch"))(x))},
    "std" = function(x, shape, skew) {c(get("pdf_std_v1", envir = asNamespace("fEGarch"))(x, df = shape))},
    "ged" = function(x, shape, skew) {c(get("pdf_ged_v1", envir = asNamespace("fEGarch"))(x, shape = shape))},
    "ald" = function(x, shape, skew) {c(get("pdf_ald_v1", envir = asNamespace("fEGarch"))(x, P = shape))},
    "snorm" = function(x, shape, skew) {get("pdf_skew_norm_s", envir = asNamespace("fEGarch"))(x, skew = skew)},
    "sstd" = function(x, shape, skew) {get("pdf_skew_sstd_s", envir = asNamespace("fEGarch"))(x, df = shape, skew = skew)},
    "sged" = function(x, shape, skew) {get("pdf_skew_sged_s", envir = asNamespace("fEGarch"))(x, shape = shape, skew = skew)},
    "sald" = function(x, shape, skew) {get("pdf_skew_sald_s", envir = asNamespace("fEGarch"))(x, P = shape, skew = skew)}
  ),
  fun2 = list(
    "norm" = function(x, mu, sigt, shape, skew) {
      c(get("pdf_norm", envir = asNamespace("fEGarch"))(x, mu = mu, sigt = sigt))
    },
    "std" = function(x, mu, sigt, shape, skew) {
      c(get("pdf_std", envir = asNamespace("fEGarch"))(x, mu = mu, sigt = sigt, df = shape))
    },
    "ged" = function(x, mu, sigt, shape, skew) {
      c(get("pdf_ged", envir = asNamespace("fEGarch"))(x, mu = mu, sigt = sigt, shape = shape))
    },
    "ald" = function(x, mu, sigt, shape, skew) {
      c(get("pdf_ald", envir = asNamespace("fEGarch"))(x, mu = mu, sigt = sigt, P = shape))
    },
    "snorm" = function(x, mu, sigt, shape, skew) {
      get("pdf_skew_norm_final", envir = asNamespace("fEGarch"))(x, mu = mu, sigt = sigt, skew = skew)
    },
    "sstd" = function(x, mu, sigt, shape, skew) {
      get("pdf_skew_sstd_final", envir = asNamespace("fEGarch"))(x, mu = mu, sigt = sigt, df = shape, skew = skew)
    },
    "sged" = function(x, mu, sigt, shape, skew) {
      get("pdf_skew_sged_final", envir = asNamespace("fEGarch"))(x, mu = mu, sigt = sigt, shape = shape, skew = skew)
    },
    "sald" = function(x, mu, sigt, shape, skew) {
      get("pdf_skew_sald_final", envir = asNamespace("fEGarch"))(x, mu = mu, sigt = sigt, P = shape, skew = skew)
    }
  ),
  # starting_par, lower bound, upper bound
  vals = list(
    "norm" = list(NULL, NULL, NULL),
    "std" = list(10, 4, Inf),
    "ged" = list(2, 1, Inf),
    "ald" = list(NULL, NULL, NULL),   # even though SM-distribution has an additional parameter, it is a non-negative integer that we circle through
    "snorm" = list(NULL, NULL, NULL),
    "sstd" = list(10, 4, Inf),
    "sged" = list(2, 1, Inf),
    "sald" = list(NULL, NULL, NULL)
  ),
  extra_par = list(
    "norm" = FALSE,
    "std" = TRUE,
    "ged" = TRUE,
    "ald" = FALSE,
    "snorm" = FALSE,
    "sstd" = TRUE,
    "sged" = TRUE,
    "sald" = FALSE
  ),
  skew_par = list(
    "norm" = FALSE,
    "std" = FALSE,
    "ged" = FALSE,
    "ald" = FALSE,
    "snorm" = TRUE,
    "sstd" = TRUE,
    "sged" = TRUE,
    "sald" = TRUE
  ),
  sim_funs = list(
    "norm" = function(x, shape, skew) {get("rnorm_s", envir = asNamespace("fEGarch"))(x)},
    "std" = function(x, shape, skew) {get("rstd_s", envir = asNamespace("fEGarch"))(x, df = shape)},
    "ged" = function(x, shape, skew) {get("rged_s", envir = asNamespace("fEGarch"))(x, shape = shape)},
    "ald" = function(x, shape, skew) {get("rald_s", envir = asNamespace("fEGarch"))(x, P = shape)},
    "snorm" = function(x, shape, skew) {get("rsnorm_s", envir = asNamespace("fEGarch"))(x, skew = skew)},
    "sstd" = function(x, shape, skew) {get("rsstd_s", envir = asNamespace("fEGarch"))(x, df = shape, skew = skew)},
    "sged" = function(x, shape, skew) {get("rsged_s", envir = asNamespace("fEGarch"))(x, shape = shape, skew = skew)},
    "sald" = function(x, shape, skew) {get("rsald_s", envir = asNamespace("fEGarch"))(x, P = shape, skew = skew)}
  ),
  sim_pars_default = list(
     mu = 0,
     omega_sig = -9,
     phi = 0.8,
     psi = numeric(0),
     kappa = -0.2,
     gamma = 0.3,
     d = 0.25,
     df = 10,
     shape = 2,
     P = 3,
     skew = 1
   ),
   sim_pars_default_fiaparch = list(
     mu = 0,
     omega = 0.0004,
     phi = 0.2,
     beta = 0.4,
     gamma = 0.1,
     delta = 2,
     d = 0.25,
     df = 10,
     P = 3,
     shape = 2,
     skew = 1
   )
)

# fitFUNS <- list(
#   "mloggarch" = fEGarch:::mloggarchFit,
#   "megarch" = fEGarch:::megarchFit,
#   "fiegarch" = fEGarch:::fiegarchFit,
#   "filoggarch" = fEGarch:::filoggarchFit,
#   "mafiloggarch" = fEGarch:::mafiloggarchFit
# )
#
# model_defaults <- data.frame(
#     Model = c("eGARCH", "meGARCH", "mlogGARCH", "filogGARCH",
#               "fieGARCH", "mafilogGARCH", "feGARCH", "logGARCH"),
#     powerTransf1 = c(1, 0, 0, 0, 1, 0, NA_real_, 0),
#     powerTransf2 = c(1, 1, 0, 0, 1, 0, NA_real_, 0),
#     longMemo = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
#     modulus1 = c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
#     modulus2 = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
#     model_type = c(2, 2, 2, 1, 2, 2, 2, 1)
#   )

# rnames <- model_defaults[, 1]
# model_defaults <- model_defaults[, -1]
# row.names(model_defaults) <- rnames

exp_funs <- list(
  "asy" = list(
    list(
      function(x, pow, shape, skew, dfun) {sign(x) * log(abs(x)) * dfun(x, shape = shape, skew = skew)},
      function(x, pow, shape, skew, dfun) {sign(x) * abs(x)^pow * dfun(x, shape = shape, skew = skew)}
    ),
    list(
      function(x, pow, shape, skew, dfun) {sign(x) * log(abs(x) + 1) * dfun(x, shape = shape, skew = skew)},
      function(x, pow, shape, skew, dfun) {sign(x) * (abs(x) + 1)^pow * dfun(x, shape = shape, skew = skew)}
    )
  ),
  "mag" = list(
    list(
      function(x, pow, shape, skew, dfun) {log(abs(x)) * dfun(x, shape = shape, skew = skew)},
      function(x, pow, shape, skew, dfun) {abs(x)^pow * dfun(x, shape = shape, skew = skew)}
    ),
    list(
      function(x, pow, shape, skew, dfun) {log(abs(x) + 1) * dfun(x, shape = shape, skew = skew)},
      function(x, pow, shape, skew, dfun) {(abs(x) + 1)^pow * dfun(x, shape = shape, skew = skew)}
    )
  )
)

mode_mat <- matrix(NA, nrow = 2, ncol = 2)
mode_mat[1, 1] <- "log_log"
mode_mat[2, 2] <- "pow_pow"
mode_mat[1, 2] <- "log_pow"
mode_mat[2, 1] <- "pow_log"

sigt_funs_egarch <- list(
  get("sigt_egarch_short_R", envir = asNamespace("fEGarch")),
  get("sigt_egarch_long_R", envir = asNamespace("fEGarch"))
)

sigt_funs_loggarch <- list(
  get("sigt_loggarch_short_R", envir = asNamespace("fEGarch")),
  get("sigt_loggarch_long_R", envir = asNamespace("fEGarch"))
)

#usethis::use_data(lookup_table, exp_funs, mode_mat, sigt_funs_egarch, sigt_funs_loggarch, internal = TRUE, overwrite = TRUE)
