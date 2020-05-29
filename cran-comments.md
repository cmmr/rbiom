## Test environments
* Windows Server 2008 (win-builder and rhub), R devel, 4.0.0, and 3.6.3
* Windows 7 and 10 (local), R 4.0.0
* Ubuntu 16.04 (travis-ci and rhub), R devel and 4.0.0
* Debian 9.3 (rhub), R devel, 4.0.0 patched, and 3.6.3
* Fedora (rhub), R devel
* CentOS 6 (rhub), R 3.5.2
* macOS 10.13.6 (rhub), R 4.0.0
* Solaris 10 (rhub), R 4.0.0


## R CMD check results
There were no ERRORs or WARNINGs.


There were 2 NOTEs:


* checking package dependencies ... NOTE
  Package suggested but not available for checking: 'rhdf5'

  The BioConductor package 'rhdf5' is an optional dependency. It's not available on all platforms.


* checking installed package size ... NOTE
  installed size is  5.5Mb
  sub-directories of 1Mb or more:
    libs   5.2Mb

  Including Rcpp and RcppParallel sometimes results in compiled libraries of this size.


## Downstream dependencies
There are no reverse dependencies.
