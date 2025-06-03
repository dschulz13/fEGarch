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
  installed size is  9.8Mb
  sub-directories of 1Mb or more:
    libs   8.3Mb

This was caused only by the macOS builder. The package heavily relies on the usage of C++ for time efficiency and therefore has a larger size.

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Dominik Schulz <dominik.schulz@uni-paderborn.de>'

  New submission


This is indeed a new package submission.


Other mentions of the win-builder include a couple of potentially misspelled words, which were rechecked to be spelled correctly, and the mention of (possibly) invalid URLs in the documentation of the object `SP500`, which were also rechecked to be correct.
