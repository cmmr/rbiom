# Human Microbiome Project - demo dataset (n = 50)

Human Microbiome Project - demo dataset (n = 50)

## Usage

``` r
hmp50
```

## Format

An rbiom object with 50 samples. Includes metadata, taxonomy, phylogeny,
and sequences.

- Sex - :

  Male or Female

- Body Site - :

  Anterior nares, Buccal mucosa, Mid vagina, Saliva, or Stool

- Age - :

  21 - 40

- BMI - :

  19 - 32

## Source

[doi:10.1101/gr.096651.109](https://doi.org/10.1101/gr.096651.109)

## See also

Other Built-In Datasets:
[`babies`](https://cmmr.github.io/rbiom/reference/babies.md),
[`gems`](https://cmmr.github.io/rbiom/reference/gems.md)

## Examples

``` r
hmp50
#> 
#> ══ Human Microbiome Project - 50 Sample Demo ═══════════════
#> 
#> Oral, nasal, vaginal, and fecal samples from a diverse set
#> of healthy volunteers. Source: Human Microbiome Project
#> (<https://hmpdacc.org>).
#> 
#>      50 Samples: HMP01, HMP02, HMP03, ..., and HMP50 
#>     490 OTUs:    Unc01yki, Unc53100, LtbAci52, ... 
#>       7 Ranks:   .otu, Kingdom, Phylum, ..., and Genus 
#>       5 Fields:  .sample, Age, BMI, Body Site, and Sex 
#>         Tree:    <present>
#> 
#> ── 182 - 22k reads/sample ──────────────────── 2023-09-22 ──
#> 
hmp50$metadata
#> # A tibble: 50 × 5
#>    .sample   Age   BMI `Body Site`    Sex   
#>  * <chr>   <dbl> <dbl> <fct>          <fct> 
#>  1 HMP01      22    20 Buccal mucosa  Female
#>  2 HMP02      24    23 Buccal mucosa  Male  
#>  3 HMP03      28    26 Saliva         Male  
#>  4 HMP04      25    23 Saliva         Male  
#>  5 HMP05      27    24 Buccal mucosa  Female
#>  6 HMP06      32    25 Saliva         Male  
#>  7 HMP07      26    22 Buccal mucosa  Male  
#>  8 HMP08      27    26 Saliva         Female
#>  9 HMP09      33    32 Saliva         Male  
#> 10 HMP10      22    20 Anterior nares Female
#> # ℹ 40 more rows
```
