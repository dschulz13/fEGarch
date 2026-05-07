window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
rt <- window.zoo(SP500, end = "2002-12-31")
.debug_mode_fEGarch <- TRUE

test_that("Bandwidth selection works as intended", {

  expect_equal({(egarch_spec() %>% fEGarch(rt, use_nonpar = TRUE))@nonpar_model$b0}, 0.189582, tolerance = 1e-04)

  expect_equal({(fiegarch_spec() %>% fEGarch(rt, use_nonpar = TRUE))@nonpar_model$b0}, 0.1893615, tolerance = 1e-04)

  expect_equal({(egarch_spec() %>% fEGarch(rt, use_nonpar = TRUE, nonparspec = locpol_spec(poly_order = 1)))@nonpar_model$b0}, 0.07757841, tolerance = 1e-04)

  expect_equal({(fiegarch_spec() %>% fEGarch(rt, use_nonpar = TRUE, nonparspec = locpol_spec(poly_order = 1)))@nonpar_model$b0}, 0.07789738, tolerance = 1e-04)


  expect_equal({(tgarch_spec() %>% fAParch(rt, use_nonpar = TRUE))@nonpar_model$b0}, 0.189582, tolerance = 1e-04)

  expect_equal({(afitgarch_spec() %>% fAParch(rt, use_nonpar = TRUE))@nonpar_model$b0}, 0.1893615, tolerance = 1e-04)

  expect_equal({(tgarch_spec() %>% fAParch(rt, use_nonpar = TRUE, nonparspec = locpol_spec(poly_order = 1)))@nonpar_model$b0}, 0.07757841, tolerance = 1e-04)

  expect_equal({(afitgarch_spec() %>% fAParch(rt, use_nonpar = TRUE, nonparspec = locpol_spec(poly_order = 1)))@nonpar_model$b0}, 0.07789738, tolerance = 1e-04)


})
