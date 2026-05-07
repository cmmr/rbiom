# Convert a variety of data types to an rbiom object.

Construct an rbiom object. The returned object is an R6 reference class.
Use `b <- a$clone()` to create copies, not `b <- a`.

## Usage

``` r
as_rbiom(biom, ...)
```

## Arguments

- biom:

  Object which can be coerced to an rbiom-class object. For example:

  *file* -

  :   Filepath or URL to a biom file.

  *matrix* -

  :   An abundance matrix with OTUs in rows and samples in columns.

  `phyloseq`-class object -

  :   From the phyloseq Bioconductor R package.

  *list* -

  :   With `counts` and optionally `metadata`, `taxonomy`, `tree`, etc
      (see details).

- ...:

  Properties to overwrite in biom: `metadata`, `taxonomy`, `tree`, etc
  (see details). Setting `underscores` here will pass it to
  [`read_tree()`](https://cmmr.github.io/rbiom/reference/read_tree.md).

## Value

An [rbiom
object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md).

## Examples

``` r
    library(rbiom)
    
    # create a simple matrix ------------------------
    mtx <- matrix(
      data     = floor(runif(24) * 1000), 
      nrow     = 6, 
      dimnames = list(paste0("OTU", 1:6), paste0("Sample", 1:4)) )
    mtx
#>      Sample1 Sample2 Sample3 Sample4
#> OTU1      80     497      34     388
#> OTU2     834     289     320     975
#> OTU3     600     732     402     289
#> OTU4     157     772     195     678
#> OTU5       7     874     403     735
#> OTU6     466     174      63     195
    
    # and some sample metadata ----------------------
    df <- data.frame(
      .sample   = paste0("Sample", 1:4),
      treatment = c("A", "B", "A", "B"),
      days      = c(12, 3, 7, 8) )
    
    # convert data set to rbiom ---------------------
    biom <- as_rbiom(mtx, metadata = df, id = "My BIOM")
    biom
#> 
#> ══ My BIOM ═════════════════════════════════════════════════
#> 
#>       4 Samples: Sample1, Sample2, Sample3, and Sample4 
#>       6 OTUs:    OTU1, OTU2, OTU3, OTU4, OTU5, and OTU6 
#>       1 Ranks:   .otu 
#>       3 Fields:  .sample, treatment, and days 
#>         Tree:    <absent>
#> 
#> ── 1.4k - 3.3k reads/sample ────────────────── 2026-05-07 ──
#> 
```
