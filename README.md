rbiom
=======

[![Travis-CI Build Status](https://travis-ci.org/cmmr/rbiom.svg?branch=master)](https://travis-ci.org/cmmr/rbiom) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rbiom)](https://cran.r-project.org/package=rbiom)

This package is a toolkit for working with Biological Observation Matrix
(BIOM) files. Features include reading/writing all BIOM formats, rarefaction,
alpha diversity, beta diversity (including UniFrac), summarizing counts by 
taxonomic level, and sample subsetting. Standalone functions for reading,
writing, and subsetting phylogenetic trees are also provided. All CPU
intensive operations are encoded in C with multi-thread support.

Reference material is available online at https://cmmr.github.io/rbiom/index.html


Installation
------------

rbiom is not currently available from CRAN, but you can install the development version from github with:

``` r
install.packages("devtools")
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
options('rbiom.max.threads' = 6)
```



