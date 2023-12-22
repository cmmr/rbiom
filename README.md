
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
#> # Test:     Estimated marginal means of linear trends.
#> # Model:    stats::lm(.abundance ~ Age * `Body Site`)
#> # A tibble: 80 × 10
#>    .taxa       `Body Site`    .trend   .se   .df .lower .upper  .t.ratio  .p.val
#>    <fct>       <fct>           <dbl> <dbl> <dbl>  <dbl>  <dbl>     <dbl>   <dbl>
#>  1 Lactobacil… Anterior n… -6.82e- 3  8.45    39  -17.1   17.1 -8.06e- 4 9.99e-1
#>  2 Lactobacil… Buccal muc… -2.68e- 2 17.6     39  -35.7   35.6 -1.52e- 3 9.99e-1
#>  3 Lactobacil… Mid vagina  -3.13e+ 1  8.42    39  -48.4  -14.3 -3.72e+ 0 6.19e-4
#>  4 Lactobacil… Saliva      -2.33e- 2 11.3     39  -22.9   22.9 -2.06e- 3 9.98e-1
#>  5 Lactobacil… Stool        4.75e-15 13.8     39  -27.9   27.9  3.44e-16 1   e+0
#>  6 Streptococ… Anterior n…  1.50e+ 0  6.15    39  -10.9   13.9  2.45e- 1 8.08e-1
#>  7 Streptococ… Buccal muc… -5.33e+ 1 12.8     39  -79.2  -27.3 -4.16e+ 0 1.71e-4
#>  8 Streptococ… Mid vagina  -1.33e- 1  6.12    39  -12.5   12.2 -2.17e- 2 9.83e-1
#>  9 Streptococ… Saliva      -4.56e+ 0  8.23    39  -21.2   12.1 -5.55e- 1 5.82e-1
#> 10 Streptococ… Stool        3.37e- 2 10.0     39  -20.3   20.3  3.35e- 3 9.97e-1
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
