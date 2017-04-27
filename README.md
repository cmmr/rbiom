rbiom
=======

[![Travis-CI Build Status](https://travis-ci.org/dansmith01/rbiom.svg?branch=master)](https://travis-ci.org/dansmith01/rbiom) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rbiom)](https://cran.r-project.org/package=rbiom)

rbiom provides utilities for working with biological observation matrix (BIOM) datasets. This includes reading/writing .biom files, rarefying counts, computing alpha/beta diversity metrics (e.g. Chao1, UniFrac), and more. Multithreading is built into the cpu intensive calculations.


Installation
------------

rbiom is not currently available from CRAN, but you can install the development version from github with:

``` r
install.packages(c("devtools", "ape", "doSNOW", "foreach", "h5", "methods", "plyr", "Rcpp", "rjson", "slam"))
devtools::install_github("cmmr/rbiom")
```


Usage
-----

``` r
library(rbiom)

infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
biom <- read.biom(infile)

# Rarefy to 1000 reads per sample
biom <- rarefy(biom, depth=1000)

# Summarize counts by phylum
phyla <- taxa.rollup(biom, 'Phylum')
phyla[1:4,1:6]

# Work with metadata
table(biom$metadata$Sex, biom$metadata$Body.Site)
sprintf("Mean age: %.1f", mean(biom$metadata$Age))

# Draw the phylogenetic tree
plot(biom$phylogeny)

# Get unifrac distance matrix
dm <- beta.div(biom, 'unifrac')
```

Several functions will by default use all available CPU cores. To limit the number of cores used, you can set the rbiom.max.threads option:

``` r
options('rbiom.max.threads') <- 6
```



