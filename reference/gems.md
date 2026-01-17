# Global Enteric Multicenter Study (n = 1,006)

Global Enteric Multicenter Study (n = 1,006)

## Usage

``` r
gems
```

## Format

An rbiom object with 1,006 samples. Includes metadata and taxonomy.

- diarrhea - :

  Case or Control

- age - :

  0 - 4.8 (years old)

- country - :

  Bangladesh, Gambia, Kenya, or Mali

## Source

[doi:10.1186/gb-2014-15-6-r76](https://doi.org/10.1186/gb-2014-15-6-r76)
and [doi:10.1093/nar/gkx1027](https://doi.org/10.1093/nar/gkx1027)

## See also

Other Built-In Datasets:
[`babies`](https://cmmr.github.io/rbiom/reference/babies.md),
[`hmp50`](https://cmmr.github.io/rbiom/reference/hmp50.md)

## Examples

``` r
gems
#> 
#> ══ Global Enteric Multicenter Study ════════════════════════
#> 
#> Stool samples collected from 1006 participants under the
#> age of 5 from Bangladesh, Gambia, Kenya, and Mali. Includes
#> 492 healthy controls and 514 moderate-to-severe diarrhea
#> (MSD) cases. Original study by Pop et al, 2014
#> (<https://doi.org/10.1186/gb-2014-15-6-r76>); dataset
#> retrieved from MicrobiomeDB
#> (<https://doi.org/10.1093/nar/gkx1027>).
#> 
#>    1006 Samples: SRS608640, SRS608279, SRS607780, ... 
#>     767 OTUs:    OTU001, OTU002, OTU003, ..., and OTU767 
#>       8 Ranks:   .otu, Kingdom, Phylum, ..., and Species 
#>       4 Fields:  .sample, diarrhea, age, and country 
#>         Tree:    <absent>
#> 
#> ── 116 - 11k reads/sample ──────────────────── 2021-03-15 ──
#> 
table(gems$metadata$country)
#> 
#> Bangladesh     Gambia      Kenya       Mali 
#>        215        272        307        212 
```
