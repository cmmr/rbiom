# Longitudinal Stool Samples from Infants (n = 2,684)

Longitudinal Stool Samples from Infants (n = 2,684)

## Usage

``` r
babies
```

## Format

An rbiom object with 2,684 samples. Includes metadata and taxonomy.

- Subject ID - :

  ID1, ID2, ..., ID12

- Sex - :

  Male or Female

- Age (days) - :

  1 - 266

- Child's diet - :

  "Breast milk", "Breast milk and formula", or "Formula"

- Sample collection - :

  "Frozen upon collection" or "Stored in alcohol"

- Antibiotic exposure - :

  Yes or No

- Antifungal exposure - :

  Yes or No

- Delivery mode - :

  Cesarean or Vaginal

- Solid food introduced (Age) - :

  116 - 247

## Source

<https://www.nature.com/articles/s41467-018-04641-7> and
[doi:10.1038/s41467-017-01973-8](https://doi.org/10.1038/s41467-017-01973-8)

## See also

Other Built-In Datasets:
[`gems`](https://cmmr.github.io/rbiom/reference/gems.md),
[`hmp50`](https://cmmr.github.io/rbiom/reference/hmp50.md)

## Examples

``` r
babies
#> 
#> ══ DailyBaby ═══════════════════════════════════════════════
#> 
#> Longitudinal stool microbiome from 12 infants (Muinck &
#> Trosvik, 2018
#> (<https://www.nature.com/articles/s41467-018-04641-7>)).
#> Retrieved from MicrobiomeHD
#> (<https://doi.org/10.1038/s41467-017-01973-8>).
#> 
#>    2684 Samples: SRR7044380, SRR7058817, SRR7039319, ... 
#>     506 OTUs:    ASV085, ASV022, ASV251, ..., and ASV607 
#>       8 Ranks:   .otu, Kingdom, Phylum, ..., and Species 
#>      10 Fields:  .sample, Subject ID, Sex, ... 
#>         Tree:    <absent>
#> 
#> ── 500 reads/sample ────────────────────────── 2018-06-08 ──
#> 
head(babies$metadata$Age)
#> Warning: Unknown or uninitialised column: `Age`.
#> NULL
```
