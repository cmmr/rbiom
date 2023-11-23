rbiom
=======

[![Travis-CI Build Status](https://travis-ci.org/cmmr/rbiom.svg?branch=master)](https://travis-ci.org/cmmr/rbiom) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rbiom)](https://cran.r-project.org/package=rbiom)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/r-rbiom/badges/version.svg)](https://anaconda.org/conda-forge/r-rbiom)

This package is a toolkit for working with Biological Observation Matrix
(BIOM) files. Features include reading/writing all BIOM formats, rarefaction,
alpha diversity, beta diversity (including UniFrac), summarizing counts by 
taxonomic level, and sample subsetting. Standalone functions for reading,
writing, and subsetting phylogenetic trees are also provided. All CPU
intensive operations are encoded in C with multi-thread support.

Reference material is available online at https://cmmr.github.io/rbiom/index.html

Source code can be found at https://github.com/cmmr/rbiom


Installation
------------

The latest stable version can be installed from CRAN.

```r
install.packages("rbiom")
```

The development version is available on GitHub.

```r
install.packages("remotes")
remotes::install_github("cmmr/rbiom")
```


Usage
-----

```r
library(rbiom)

infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
biom <- read_biom(infile)

# Rarefy to 1000 reads per sample
biom <- sample_rarefy(biom, depth=1000)

# Summarize counts by phylum
phyla <- taxa_rollup(biom, 'Phylum')
phyla[1:4,1:6]

# Work with metadata
table(biom$metadata$Sex, biom$metadata$Body.Site)
sprintf("Mean age: %.1f", mean(biom$metadata$Age))

# Draw the phylogenetic tree
plot(biom$phylogeny)

# Get unifrac distance matrix
dm <- bdiv_distmat(biom, 'unifrac')
```


Parallel Processing
-------------------

Computation of beta diversity metrics (UniFrac, Bray-Curtis, etc) will use all available CPU cores by default. To limit the number of cores used, you can set the numThreads option:

```r
RcppParallel::setThreadOptions(numThreads = 4)
```


Caching
-------

Caching is enabled by default. rbiom will store a maximum of 200MB in the temporary directory given by `file.path(tempdir(), "rbiom", "cache")`. The following commands can be used to change the cache directory, storage limit (given in bytes), and key hashing function:
```r
options(rbiom.cache_dir="/tmp/rbiom_cache")
options(rbiom.cache_size=5 * 1024 ^ 2) # 5MB
options(rbiom.cache_hash=rlang::hash)

Sys.setenv(RBIOM_CACHE_DIR="/tmp/rbiom_cache")
Sys.setenv(RBIOM_CACHE_SIZE=1024 ^ 2) # 1GB
```

Setting the cache directory to `"FALSE"` will disable caching.

R options will override environment variables.



Building from source
--------------------

rbiom requires the following system libraries which can be installed through your operating system's package manager.

* deb (Debian, Ubuntu): `libudunits2-dev libssl-dev libxml2-dev libcurl4-openssl-dev libgdal-dev`
* rpm (Fedora, CentOS, RHEL): `udunits2-devel openssl-devel libxml2-devel libcurl-devel gdal-devel`
* csw (Solaris): `libssl_dev openssl@1.1 libxml2_dev gdal_dev`
* brew (OSX): `udunits`



