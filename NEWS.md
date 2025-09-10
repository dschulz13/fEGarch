# fEGarch 1.0.2
- a bug was fixed, where for the nonparametric scale estimation always the
  automated bandwidth selection under long-memory errors was employed
- a safety measure was included to avoid that the log-likelihood can
  become -Inf
- math formulas were fully removed from the README, because they were not
  displayed properly on CRAN

# fEGarch 1.0.1
- a bug was fixed, where for most models except for EGARCH-type models the
  mean estimate was not considered properly in volatility forecasts
- a bug was fixed, where for EGARCH-family models rolling point forecasts
  under the (skewed) ALD were erroneous
- the mathematical formula in the README was not displayed correctly on the
  CRAN servers and was therefore removed entirely
- the package description was adjusted to contain also information on the
  the more advanced models and backtesting capabilities of the package,
  making it easier for users to see what to expect from the package
- the package title was adjusted to contain more information on the
  package's contents, making it easier for users to see what to expect from
  the package
- a bug was fixed in the computation of the hessian matrix, where now
  it is caught if a hessian is not invertible through solve() and
  the fitting functions therefore now don't stop due to an error anymore
- an additional dataset `UKinflation` was added to the package to show
  applications of the package beyond return data
- some examples in the README were adjusted / added
