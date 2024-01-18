
<!-- Run `devtools::build_readme(); pkgdown::build_home()` after editing.  -->

# rbiom

<!-- badges: start -->

[![Travis-CI Build
Status](https://travis-ci.org/cmmr/rbiom.svg?branch=master)](https://travis-ci.org/cmmr/rbiom)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rbiom)](https://cran.r-project.org/package=rbiom)
[![Anaconda-Server
Badge](https://anaconda.org/conda-forge/r-rbiom/badges/version.svg)](https://anaconda.org/conda-forge/r-rbiom)
<!-- badges: end -->

This package is a toolkit for working with Biological Observation Matrix
(BIOM) files. Features include reading/writing all ‘BIOM’ formats,
rarefaction, alpha diversity, beta diversity (including ‘UniFrac’),
summarizing counts by taxonomic level, subsetting, visualizations, and
statistical analysis. All CPU intensive operations are written in C.

Reference material is available online at
<https://cmmr.github.io/rbiom/index.html>

Source code can be found at <https://github.com/cmmr/rbiom>

## Installation

The latest stable version can be installed from CRAN.

``` r
install.packages("rbiom")
```

The development version is available on GitHub.

``` r
install.packages("remotes")
remotes::install_github("cmmr/rbiom")
```

## Usage

#### Import and rarefy abundance counts.

``` r
library(rbiom)

infile <- system.file(package = "rbiom", "extdata", "hmp50.bz2")
biom   <- rarefy(infile)
```

#### Explore associations with metadata.

``` r
bdiv_ord_plot(biom, color.by = "Body Site", facet.by = "Sex")
```

![](man/figures/README-bdiv-1.png)<!-- -->

``` r
adiv_boxplot(biom, x = "Sex", adiv = c("otu", "shan"), color.by = "Body Site")
```

![](man/figures/README-bdiv-2.png)<!-- -->

#### Compute statistics for correlation models.

``` r
taxa_corrplot(biom, x = "Age", layers = "stc", taxa = .01, p.top = 4, color.by = "bod")
```

![](man/figures/README-stats-1.png)<!-- -->

``` r
taxa_stats(biom, regr = "Age", stat.by = "Body Site", taxa = 0.01)
#> # Test:     Is each trendline's slope non-zero?
#> # Model:    stats::lm(.abundance ~ Age * `Body Site`)
#> # A tibble: 80 × 10
#>    .taxa      `Body Site`    .trend   .se   .df  .lower .upper  .t.ratio  .p.val
#>    <fct>      <fct>           <dbl> <dbl> <dbl>   <dbl>  <dbl>     <dbl>   <dbl>
#>  1 Lactobaci… Anterior n… -3.47e- 1  1.25    39  -2.88    2.19 -2.77e- 1 0.783  
#>  2 Lactobaci… Buccal muc… -5.09e- 1  2.61    39  -5.79    4.78 -1.95e- 1 0.847  
#>  3 Lactobaci… Mid vagina  -4.29e+ 0  1.25    39  -6.81   -1.76 -3.44e+ 0 0.00142
#>  4 Lactobaci… Saliva      -4.42e- 1  1.68    39  -3.84    2.95 -2.63e- 1 0.794  
#>  5 Lactobaci… Stool        8.33e-16  2.05    39  -4.14    4.14  4.07e-16 1      
#>  6 Streptoco… Anterior n…  3.71e+ 0  2.05    39  -0.436   7.85  1.81e+ 0 0.0780 
#>  7 Streptoco… Buccal muc… -3.01e+ 0  4.27    39 -11.6     5.63 -7.04e- 1 0.486  
#>  8 Streptoco… Mid vagina  -1.60e+ 0  2.04    39  -5.72    2.53 -7.83e- 1 0.438  
#>  9 Streptoco… Saliva      -3.38e+ 0  2.74    39  -8.92    2.17 -1.23e+ 0 0.225  
#> 10 Streptoco… Stool        7.56e- 1  3.35    39  -6.01    7.52  2.26e- 1 0.822  
#> # ℹ 70 more rows
#> # ℹ 1 more variable: .adj.p <dbl>
```

#### Summarize counts by taxonomic rank.

``` r
taxa_heatmap(biom, taxa = 30, color.by = c("body", "age"), limit.by = c(sex = "Male"))
```

![](man/figures/README-taxa-1.png)<!-- -->

``` r
taxa_barplot(biom, rank = "Phylum")
```

![](man/figures/README-taxa-2.png)<!-- -->

``` r
taxa_table(biom, 'Phylum')
#> # A tibble: 637 × 8
#>    .rank  .sample .taxa               .abundance   Age   BMI `Body Site`   Sex  
#>    <fct>  <chr>   <fct>                    <dbl> <dbl> <dbl> <fct>         <fct>
#>  1 Phylum HMP01   Actinobacteria              13    22    20 Buccal mucosa Fema…
#>  2 Phylum HMP01   Bacteroidetes              192    22    20 Buccal mucosa Fema…
#>  3 Phylum HMP01   Cyanobacteria                0    22    20 Buccal mucosa Fema…
#>  4 Phylum HMP01   Deinococcus Thermus          0    22    20 Buccal mucosa Fema…
#>  5 Phylum HMP01   Firmicutes                 854    22    20 Buccal mucosa Fema…
#>  6 Phylum HMP01   Fusobacteria                37    22    20 Buccal mucosa Fema…
#>  7 Phylum HMP01   Gracilibacteria             13    22    20 Buccal mucosa Fema…
#>  8 Phylum HMP01   Proteobacteria              74    22    20 Buccal mucosa Fema…
#>  9 Phylum HMP01   Saccharibacteria             0    22    20 Buccal mucosa Fema…
#> 10 Phylum HMP01   Spirochaetae                 0    22    20 Buccal mucosa Fema…
#> # ℹ 627 more rows
```

## Parallel Processing

Computation of beta diversity metrics (UniFrac, Bray-Curtis, etc) will
use all available CPU cores by default. To limit the number of cores
used, you can set the numThreads option:

``` r
RcppParallel::setThreadOptions(numThreads = 4)
```

## Building from source

rbiom requires the following system libraries which can be installed
through your operating system’s package manager.

- deb (Debian, Ubuntu):
  `libudunits2-dev libssl-dev libxml2-dev libcurl4-openssl-dev libgdal-dev`
- rpm (Fedora, CentOS, RHEL):
  `udunits2-devel openssl-devel libxml2-devel libcurl-devel gdal-devel`
- csw (Solaris): `libssl_dev openssl@1.1 libxml2_dev gdal_dev`
- brew (OSX): `udunits`
