# Combine several rbiom objects into one.

WARNING: It is generally ill-advised to merge BIOM datasets, as OTUs
mappings are dependent on upstream clustering and are not equivalent
between BIOM files.

## Usage

``` r
biom_merge(
  ...,
  metadata = NA,
  taxonomy = NA,
  tree = NULL,
  sequences = NA,
  id = NA,
  comment = NA
)
```

## Arguments

- ...:

  Any number of rbiom objects (e.g. from
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md)),
  lists of rbiom objects, or valid arguments to the `biom` parameter of
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md)
  (for instance file names).

- metadata, taxonomy, tree, sequences, id, comment:

  Replace the corresponding data in the merged rbiom object with these
  values. Set to `NULL` to not inherit a particular component. The
  default, `NA`, will attempt to create the component based on `...`
  values. The merged phylogenetic tree cannot be inferred.

## Value

An [rbiom
object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md).

## See also

Other biom: [`bdply()`](https://cmmr.github.io/rbiom/reference/bdply.md)

## Examples

``` r
    library(rbiom)
    
    b1 <- as_rbiom(hmp50$counts[,1:4])
    b2 <- as_rbiom(hmp50$counts[,5:8])
    
    biom <- biom_merge(b1, b2)
    print(biom)
#> 
#> ══ Merged BIOM ═════════════════════════════════════════════
#> 
#>       8 Samples: HMP01, HMP02, HMP03, ..., and HMP08 
#>     168 OTUs:    Unc53100, PpbAcne6, UncO2012, ... 
#>       1 Ranks:   .otu 
#>       1 Fields:  .sample 
#>         Tree:    <absent>
#> 
#> ── 1.4k - 4.2k reads/sample ────────────────── 2026-01-23 ──
#> 
    
    biom$tree     <- hmp50$tree
    biom$metadata <- hmp50$metadata
#> Warning: ℹ Ignoring metadata for 42 samples not currently in biom object: "HMP09",
#>   "HMP10", "HMP11", "HMP12", "HMP13", "HMP14", "HMP15", "HMP16", "HMP17",
#>   "HMP18", "HMP19", "HMP20", "HMP21", "HMP22", "HMP23", "HMP24", "HMP25",
#>   "HMP26", …, "HMP49", and "HMP50".
    print(biom)
#> 
#> ══ Merged BIOM ═════════════════════════════════════════════
#> 
#>       8 Samples: HMP01, HMP02, HMP03, ..., and HMP08 
#>     168 OTUs:    Unc53100, PpbAcne6, UncO2012, ... 
#>       1 Ranks:   .otu 
#>       5 Fields:  .sample, Age, BMI, Body Site, and Sex 
#>         Tree:    <present>
#> 
#> ── 1.4k - 4.2k reads/sample ────────────────── 2026-01-23 ──
#> 
```
