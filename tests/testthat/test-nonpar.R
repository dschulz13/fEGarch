window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
rt <- window.zoo(SP500, end = "2002-12-31")
spec <- egarch_spec()
spec2 <- fiegarch_spec()

test_that("Bandwidth selections are employed for short and long memory correctly", {
  expect_equal(
   fEGarch(spec, rt, use_nonpar = TRUE)@nonpar_model$b0, 0.189582019547447
  )
  expect_equal(
   fEGarch(spec2, rt, use_nonpar = TRUE)@nonpar_model$b0, 0.189361545957288
  )


})
