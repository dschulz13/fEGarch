# Check for new package version patch 1.0.2

## Test environments
* local Windows install, R 4.4.2
* win-builder (oldrelease, release and devel)
* macOS builder
* R-hub (linux, m1-san, macos, macos-arm64, windows, atlas, c23, clang-asan, clang-ubsan, clang16, clang17, clang18, clang19, donttest, gcc-asan, gcc13, gcc14, gcc15, mkl, nold, noremap, nosuggests, ubuntu-clang, ubuntu-gcc12, ubuntu-next, ubuntu-release)

## R CMD check results

There were no errors.

There were no warnings.

There was one INFO:

* checking installed package size ... INFO
  installed size is 10.0Mb
  sub-directories of 1Mb or more:
    libs   8.3Mb

This was caused only by the macOS builder. The package heavily relies on the usage of C++ for time efficiency and therefore has a larger size.

There was one NOTE:

* checking CRAN incoming feasibility ... [112s] NOTE
  Maintainer: 'Dominik Schulz <dominik.schulz@uni-paderborn.de>'
  
Dominik Schulz with that email address is indeed the correct maintainer of the package.

# Checks for old package version 1.0.1

## Test environments
* local Windows install, R 4.4.2
* win-builder (oldrelease, release and devel)
* macOS builder
* R-hub (linux, m1-san, macos, macos-arm64, windows, atlas, c23, clang-asan, clang-ubsan, clang16, clang17, clang18, clang19, donttest, gcc-asan, gcc13, gcc14, gcc15, mkl, nold, noremap, nosuggests, ubuntu-clang, ubuntu-gcc12, ubuntu-next, ubuntu-release)

## R CMD check results

There were no errors.

There were no warnings.

There was two NOTEs:

* checking CRAN incoming feasibility ... [25s] NOTE
  Maintainer: 'Dominik Schulz <dominik.schulz@uni-paderborn.de>'

This is due to the quick patch after the first submission of the
package. An important bug for one of the main parts of the package was found that
affects volatility forecasts. Please excuse this quick but urgent
patch.

* checking installed package size ... NOTE
  installed size is 10.0Mb
  sub-directories of 1Mb or more:
    libs   8.3Mb

This was caused only by the macOS builder. The package heavily relies on the usage of C++ for time efficiency and therefore has a larger size.

# Checks for old package version 1.0.0

## Test environments
* local Windows install, R 4.4.2
* win-builder (oldrelease, release and devel)
* macOS builder
* R-hub (linux, m1-san, macos, macos-arm64, windows, atlas, c23, clang-asan, clang-ubsan, clang16, clang17, clang18, clang19, donttest, gcc-asan, gcc13, gcc14, gcc15, mkl, nold, noremap, nosuggests, ubuntu-clang, ubuntu-gcc12, ubuntu-next, ubuntu-release)

## R CMD check results

There were no errors.

There were no warnings.

There were two NOTEs:

* checking installed package size ... NOTE
  installed size is  9.9Mb
  sub-directories of 1Mb or more:
    libs   8.3Mb

This was caused only by the macOS builder. The package heavily relies on the usage of C++ for time efficiency and therefore has a larger size.

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Dominik Schulz <dominik.schulz@uni-paderborn.de>'

  New submission


This is indeed a new package submission.


Other mentions of the win-builder include a couple of potentially misspelled words, which were rechecked to be spelled correctly, and the mention of (possibly) invalid URLs in the documentation of the object `SP500`, which were also rechecked to be correct.

## Addressing Comments on Initial Submission

The first submission of the package underwent a manual inspection by a CRAN
member and was rejected with further comments. Find here the list of
comments by the CRAN member and how the issues were addressed.

One comment received:

* Please add \value to .Rd files regarding exported methods and explain
  the functions results in the documentation. Please write about the
  structure of the output (class) and also what the output means. (If a
  function does not return a value, please document that too, e.g.
  \value{No return value, called for side effects} or similar)
  For more details:
  <https://contributor.r-project.org/cran-cookbook/docs_issues.html#missing-value-tags-in-.rd-files>
  Missing Rd-tags:
      loss_functions-fEGarch_risk-method.Rd: \value

A "Value" section was added to the documentation of
"loss_functions-fEGarch_risk-method" explaining the output
of the method.
