## Test environments
* Windows 10 (local), R-release
* Windows Server 2008, (winbuilder and rhub), R-oldrel, R-release, R-devel
* macOS High Sierra, (rhub), R-release
* macOS Big Sur, (rhub), R-release
* Ubuntu, (rhub), R-release, R-devel
* Fedora, (rhub), R-devel
* Debian, (rhub), R-release, R-patched, R-devel
* Solaris, (rhub), R-release


## R CMD check results
There were no ERRORs or WARNINGs.


There were 2 NOTEs:


* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'rhdf5'

  The BioConductor package 'rhdf5' is an optional dependency. It's not available on all platforms.


* checking installed package size ... NOTE
  installed size is  5.8Mb
  sub-directories of 1Mb or more:
    libs   5.4Mb

  Including Rcpp and RcppParallel sometimes results in compiled libraries of this size.


## Downstream dependencies
There are no reverse dependencies.
