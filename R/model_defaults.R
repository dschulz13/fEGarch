model_defaults <- data.frame(
  powerTransf1 = c(FALSE, FALSE),
  powerTransf2 = c(FALSE, FALSE),
  longMemo = c(FALSE, TRUE),
  modulus = c(FALSE, FALSE)
)

row.names(model_defaults) <- c("eGARCH", "fieGARCH")
