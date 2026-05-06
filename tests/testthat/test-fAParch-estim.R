window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
rt <- window.zoo(SP500, end = "2002-12-31")

test_that("fAParch_spec functions work as intended", {
  expect_no_error({
    fAParch_spec(model_type = "garch", orders = c(2, 2), long_memo = FALSE,
                 cond_dist = "ald", adjust = FALSE)
  })
  expect_no_error({
    fAParch_spec(model_type = "aparch", orders = c(1, 1), long_memo = TRUE,
                 cond_dist = "ald", adjust = TRUE)
  })
  expect_no_error({
    afiaparch_spec(orders = c(1, 1), cond_dist = "ald")
  })
  expect_message({
    fAParch_spec(orders = c(2, 1), cond_dist = "std", long_memo = TRUE, adjust = TRUE)
  }, "Only orders p=1 with q=1 are currently supported for the models with fractional differencing by fAParch_spec.")
  expect_message({
    fAParch_spec(orders = c(2, 1), cond_dist = "std", long_memo = TRUE, adjust = FALSE)
  }, "Only orders p=1 with q=1 are currently supported for the models with fractional differencing by fAParch_spec.")
  expect_message({
    afiaparch_spec(orders = c(1, 2), cond_dist = "ald")
  }, "Only orders p=1 with q=1 are currently supported for the models with fractional differencing by fAParch_spec.")

  expect_equal({suppressMessages(afiaparch_spec(orders = c(1, 2), cond_dist = "ald"))@orders}, c(1, 1))
  expect_equal({suppressMessages(fAParch_spec(orders = c(3, 2), cond_dist = "ald", long_memo = TRUE, adjust = TRUE))@orders}, c(1, 1))
  expect_equal({suppressMessages(fAParch_spec(orders = c(3, 2), cond_dist = "ald", long_memo = TRUE, adjust = FALSE))@orders}, c(1, 1))

})

test_that("fAParch function works as intended", {

  expect_no_error({
   garch_spec() %>% fAParch(rt)
  })

  expect_no_error({
   suppressWarnings(aparch_spec() %>% fAParch(rt))
  })

  expect_no_error({
   gjrgarch_spec() %>% fAParch(rt)
  })

  expect_no_error({
   figarch_spec() %>% fAParch(rt)
  })

  expect_no_error({
   suppressWarnings(fiaparch_spec() %>% fAParch(rt))
  })

  expect_no_error({
   figjrgarch_spec() %>% fAParch(rt)
  })

  expect_no_error({
   fitgarch_spec() %>% fAParch(rt)
  })

  expect_no_error({
   fingarch_spec() %>% fAParch(rt)
  })

  expect_no_error({
   tgarch_spec() %>% fAParch(rt)
  })

  expect_no_error({
   ngarch_spec() %>% fAParch(rt)
  })

  expect_no_error({
   figarch_spec() %>% fAParch(rt, trunc = 50)
  })

  expect_no_error({
   suppressWarnings(fiaparch_spec() %>% fAParch(rt, trunc = 50))
  })

  expect_no_error({
   suppressWarnings(figjrgarch_spec() %>% fAParch(rt, trunc = 50))
  })

  expect_no_error({
   suppressWarnings(fitgarch_spec() %>% fAParch(rt, trunc = 50))
  })

  expect_no_error({
   suppressWarnings(fingarch_spec() %>% fAParch(rt, trunc = 50))
  })

  expect_no_error({
   suppressWarnings(afiaparch_spec() %>% fAParch(rt))
  })

  expect_no_error({
   suppressWarnings(afitgarch_spec() %>% fAParch(rt))
  })

  expect_no_error({
   suppressWarnings(afigjrgarch_spec() %>% fAParch(rt))
  })

  expect_no_error({
   suppressWarnings(afingarch_spec() %>% fAParch(rt))
  })

  expect_no_error({
   suppressWarnings(afiaparch_spec() %>% fAParch(rt, trunc = 50))
  })

  expect_no_error({
   suppressWarnings(afitgarch_spec() %>% fAParch(rt, trunc = 50))
  })

  expect_equal({
   suppressWarnings(afitgarch_spec() %>% fAParch(rt, trunc = 50))@trunc
  }, 50)

  expect_equal({
   suppressWarnings(afitgarch_spec() %>% fAParch(rt))@trunc
  }, "none")

  expect_no_error({
   suppressWarnings(afigjrgarch_spec() %>% fAParch(rt, trunc = 50))
  })

  expect_no_error({
   suppressWarnings(afingarch_spec() %>% fAParch(rt, trunc = 50))
  })

  expect_no_error({
   garch_spec(cond_dist = "norm") %>% fAParch(rt)
  })

  expect_no_error({
   garch_spec(cond_dist = "std") %>% fAParch(rt)
  })

  expect_no_error({
   garch_spec(cond_dist = "ged") %>% fAParch(rt)
  })

  expect_no_error({
   suppressWarnings(garch_spec(cond_dist = "ald") %>% fAParch(rt, parallel = FALSE, Prange = c(1, 3)))
  })

  expect_no_error({
   suppressWarnings(garch_spec(cond_dist = "sald") %>% fAParch(rt, parallel = FALSE, Prange = c(1, 3)))
  })

  expect_no_error({
   garch_spec(cond_dist = "snorm") %>% fAParch(rt)
  })

  expect_no_error({
   garch_spec(cond_dist = "sstd") %>% fAParch(rt)
  })

  expect_no_error({
   garch_spec(cond_dist = "sged") %>% fAParch(rt)
  })

  expect_no_error({
   garch_spec(cond_dist = "norm", orders = c(1, 1)) %>% fAParch(rt)
  })

  expect_no_error({
   garch_spec(cond_dist = "norm", orders = c(2, 1)) %>% fAParch(rt)
  })

  expect_no_error({
   garch_spec(cond_dist = "norm", orders = c(1, 2)) %>% fAParch(rt)
  })

  expect_equal({length(is(afitgarch_spec(cond_dist = "norm", orders = c(1, 1)) %>% fAParch(rt)))}, 2)

  expect_equal({is(afitgarch_spec(cond_dist = "norm", orders = c(1, 1)) %>% fAParch(rt))}, c("fEGarch_fit_afitgarch", "fEGarch_fit"))

  expect_equal({length(is(afiaparch_spec(cond_dist = "norm", orders = c(1, 1)) %>% fAParch(rt)))}, 2)

  expect_equal({is(afiaparch_spec(cond_dist = "norm", orders = c(1, 1)) %>% fAParch(rt))}, c("fEGarch_fit_afiaparch", "fEGarch_fit"))

  expect_equal({length(is(afigjrgarch_spec(cond_dist = "norm", orders = c(1, 1)) %>% fAParch(rt)))}, 2)

  expect_equal({is(afigjrgarch_spec(cond_dist = "norm", orders = c(1, 1)) %>% fAParch(rt))}, c("fEGarch_fit_afigjrgarch", "fEGarch_fit"))

  expect_equal({length(is(afingarch_spec(cond_dist = "norm", orders = c(1, 1)) %>% fAParch(rt)))}, 2)

  expect_equal({is(afingarch_spec(cond_dist = "norm", orders = c(1, 1)) %>% fAParch(rt))}, c("fEGarch_fit_afingarch", "fEGarch_fit"))


})
