window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
rt <- window.zoo(SP500, end = "2002-12-31")

test_that("fAParch forecasting works as intended", {
  expect_no_error({
    fAParch_spec(model_type = "garch", orders = c(2, 2), long_memo = FALSE,
                 cond_dist = "ald", adjust = FALSE)
  })

})
