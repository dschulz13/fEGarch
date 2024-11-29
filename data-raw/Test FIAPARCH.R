library(zoo) # for time series handling
library(yfR) # for downloading data
library(magrittr) # for pipe operator "%>%"
library(fEGarch) # for model fitting
data0 <- yf_get("^GDAXI", first_date = "2010-01-01",
last_date = "2024-10-31")
rt <- data0$price_close %>%
zoo(order.by = as.Date(data0$ref_date)) %>%
log() %>%
diff()

#test <- lmgarchFit(rt)


test <- fiaparch(rt, cond_dist = "std")
library(rugarch)
spec <- ugarchspec(
  variance.model = list(
    model = "apARCH", garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 0)
  ),
  distribution.model = "sged"
)
fit <- ugarchfit(spec, coredata(rt))
