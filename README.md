
<!-- Run `devtools::build_readme()` after editing.  -->

# rbiom <img src="man/figures/logo.png" align="right" width="174" height="200" alt="rbiom logo" />

<!-- badges: start -->

[![cran](https://www.r-pkg.org/badges/version/rbiom)](https://CRAN.R-project.org/package=rbiom)
[![conda](https://anaconda.org/conda-forge/r-rbiom/badges/version.svg)](https://anaconda.org/conda-forge/r-rbiom)
[![downloads](https://cranlogs.r-pkg.org/badges/grand-total/rbiom)](https://cranlogs.r-pkg.org/)
[![covr](https://codecov.io/gh/cmmr/rbiom/graph/badge.svg)](https://app.codecov.io/gh/cmmr/rbiom)
<!-- badges: end -->

`rbiom` is designed for microbiome researchers, providing visualizations
and statistical analyses from Biological Observation Matrix (BIOM)
files.

## Installation

The latest stable version can be installed from CRAN.

``` r
install.packages('rbiom')
```

The development version is available on GitHub.

``` r
install.packages('pak')
pak::pak('cmmr/rbiom')
```

## Usage

#### Import and rarefy abundance counts.

``` r
library(rbiom)

infile <- system.file(package = 'rbiom', 'extdata', 'hmp50.bz2')
biom   <- rarefy(infile)
```

#### Explore associations with metadata.

``` r
bdiv_ord_plot(biom, stat.by = 'Body Site', facet.by = 'Sex')
```

![](man/figures/README-bdiv-1.png)<!-- -->

``` r
adiv_boxplot(biom, x = 'Sex', adiv = c('otu', 'shan'), stat.by = 'Body Site')
```

![](man/figures/README-bdiv-2.png)<!-- -->

``` r
subset(biom, `Body Site` == 'Buccal mucosa') %>% 
  taxa_corrplot('Age', taxa = 2, layers = 'ptc', fit = 'lm', test = 'emtrends')
```

![](man/figures/README-bdiv-3.png)<!-- -->

#### Summarize counts by taxonomic rank.

``` r
taxa_heatmap(biom, taxa = 10, tracks = c('body', 'age'))
```

![](man/figures/README-taxa-1.png)<!-- -->

``` r
taxa_stacked(biom, rank = 'Phylum')
```

![](man/figures/README-taxa-2.png)<!-- -->

``` r
taxa_table(biom, 'Phylum')
#> # A tibble: 294 × 8
#>    .rank  .sample .taxa          .abundance   Age   BMI `Body Site`   Sex   
#>    <fct>  <chr>   <fct>               <dbl> <dbl> <dbl> <fct>         <fct> 
#>  1 Phylum HMP01   Firmicutes            856    22    20 Buccal mucosa Female
#>  2 Phylum HMP01   Bacteroidetes         199    22    20 Buccal mucosa Female
#>  3 Phylum HMP01   Actinobacteria         16    22    20 Buccal mucosa Female
#>  4 Phylum HMP01   Proteobacteria         72    22    20 Buccal mucosa Female
#>  5 Phylum HMP01   Fusobacteria           32    22    20 Buccal mucosa Female
#>  6 Phylum HMP01   Tenericutes             0    22    20 Buccal mucosa Female
#>  7 Phylum HMP02   Firmicutes            803    24    23 Buccal mucosa Male  
#>  8 Phylum HMP02   Bacteroidetes         192    24    23 Buccal mucosa Male  
#>  9 Phylum HMP02   Actinobacteria         52    24    23 Buccal mucosa Male  
#> 10 Phylum HMP02   Proteobacteria         96    24    23 Buccal mucosa Male  
#> # ℹ 284 more rows
```

## Documentation

The online manual for `rbiom` is available at
<https://cmmr.github.io/rbiom/>. It includes a getting started guide,
articles that explore specific use cases, and reference pages for each
function.

## Community guidelines

### Support

Bug reports, feature requests, and general questions can be submitted at
<https://github.com/cmmr/rbiom/issues>.

### Contributing

Pull requests are welcome. Please ensure contributed code is covered by
tests and documentation (add additional tests and documentation as
needed) and that it passes all automated tests.

## Automated tests

The following commands will check if `rbiom` passes the bundled testing
suite.

``` r
install.packages('testthat')
testthat::test_check('rbiom')
```
