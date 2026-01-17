# Subset an rbiom object by sample names, OTU names, metadata, or taxonomy.

Dropping samples or OTUs will lead to observations being removed from
the OTU matrix (`biom$counts`). OTUs and samples with zero observations
are automatically removed from the rbiom object.

## Usage

``` r
# S3 method for class 'rbiom'
subset(x, subset, clone = TRUE, ...)

# S3 method for class 'rbiom'
x[i, j, ..., clone = TRUE, drop = FALSE]

# S3 method for class 'rbiom'
na.omit(object, fields = ".all", clone = TRUE, ...)

subset_taxa(x, subset, clone = TRUE, ...)
```

## Arguments

- x:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), such
  as from
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- subset:

  Logical expression for rows to keep. See
  [`base::subset()`](https://rdrr.io/r/base/subset.html).

- clone:

  Create a copy of `biom` before modifying. If `FALSE`, `biom` is
  modified in place as a side-effect. See [speed
  ups](https://cmmr.github.io/rbiom/reference/speed.md) for use cases.
  Default: `TRUE`

- ...:

  Not used.

- i, j:

  The sample or OTU names to keep. Or a logical/integer vector
  indicating which sample names from `biom$samples` or `biom$otus` to
  keep. Subsetting with `[i]` takes `i` as samples, whereas `[i,j]`
  takes `i` as otus and `j` as samples (corresponding to `[rows, cols]`
  in the underlying `biom$counts` matrix).

- drop:

  Not used

- object:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), such
  as from
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- fields:

  Which metadata field(s) to check for `NA`s, or `".all"` to check all
  metadata fields.

## Value

An [rbiom
object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md).

## See also

Other transformations:
[`matrix_ops`](https://cmmr.github.io/rbiom/reference/matrix_ops.md),
[`modify_metadata`](https://cmmr.github.io/rbiom/reference/modify_metadata.md),
[`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md),
[`slice_metadata`](https://cmmr.github.io/rbiom/reference/slice_metadata.md),
[`with()`](https://cmmr.github.io/rbiom/reference/with.md)

## Examples

``` r
    library(rbiom)
    library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
    
    # Subset to specific samples
    biom <- hmp50[c('HMP20', 'HMP42', 'HMP12')]
    biom$metadata
#> # A tibble: 3 × 5
#>   .sample   Age   BMI `Body Site` Sex   
#> * <chr>   <dbl> <dbl> <fct>       <fct> 
#> 1 HMP20      27    22 Stool       Female
#> 2 HMP42      34    19 Mid vagina  Female
#> 3 HMP12      35    26 Stool       Male  
    
    # Subset to specific OTUs
    biom <- hmp50[c('LtbAci52', 'UncO2012'),] # <- Trailing ,
    biom$taxonomy
#> # A tibble: 2 × 7
#>   .otu     Kingdom  Phylum        Class       Order           Family       Genus
#> * <chr>    <fct>    <fct>         <fct>       <fct>           <fct>        <fct>
#> 1 LtbAci52 Bacteria Firmicutes    Bacilli     Lactobacillales Lactobacill… Lact…
#> 2 UncO2012 Bacteria Bacteroidetes Bacteroidia Bacteroidales   Bacteroidac… Bact…
    
    # Subset to specific samples and OTUs
    biom <- hmp50[c('LtbAci52', 'UncO2012'), c('HMP20', 'HMP42', 'HMP12')]
    as.matrix(biom)
#>          HMP20 HMP42 HMP12
#> LtbAci52     0  1167     0
#> UncO2012   119     1  1196
    
    # Subset samples according to metadata
    biom <- subset(hmp50, `Body Site` %in% c('Saliva') & Age < 25)
    biom$metadata
#> # A tibble: 3 × 5
#>   .sample   Age   BMI `Body Site` Sex   
#> * <chr>   <dbl> <dbl> <fct>       <fct> 
#> 1 HMP27      24    30 Saliva      Female
#> 2 HMP28      23    19 Saliva      Female
#> 3 HMP30      24    21 Saliva      Female
    
    # Subset OTUs according to taxonomy
    biom <- subset_taxa(hmp50, Phylum == 'Cyanobacteria')
    biom$taxonomy
#> # A tibble: 4 × 7
#>   .otu     Kingdom  Phylum        Class           Order             Family Genus
#> * <chr>    <fct>    <fct>         <fct>           <fct>             <fct>  <fct>
#> 1 Hu4Lup30 Bacteria Cyanobacteria Chloroplast     o                 f      g    
#> 2 Unc02oth Bacteria Cyanobacteria Melainabacteria Obscuribacterales f      g    
#> 3 PinJeffr Bacteria Cyanobacteria Chloroplast     o                 f      g    
#> 4 Unc48787 Bacteria Cyanobacteria Melainabacteria Gastranaerophila… f      g    
    
    # Remove samples with NA metadata values
    biom <- mutate(hmp50, BS2 = na_if(`Body Site`, 'Saliva'))
    biom$metadata
#> # A tibble: 50 × 6
#>    .sample   Age   BMI `Body Site`    Sex    BS2           
#>  * <chr>   <dbl> <dbl> <fct>          <fct>  <fct>         
#>  1 HMP01      22    20 Buccal mucosa  Female Buccal mucosa 
#>  2 HMP02      24    23 Buccal mucosa  Male   Buccal mucosa 
#>  3 HMP03      28    26 Saliva         Male   NA            
#>  4 HMP04      25    23 Saliva         Male   NA            
#>  5 HMP05      27    24 Buccal mucosa  Female Buccal mucosa 
#>  6 HMP06      32    25 Saliva         Male   NA            
#>  7 HMP07      26    22 Buccal mucosa  Male   Buccal mucosa 
#>  8 HMP08      27    26 Saliva         Female NA            
#>  9 HMP09      33    32 Saliva         Male   NA            
#> 10 HMP10      22    20 Anterior nares Female Anterior nares
#> # ℹ 40 more rows
    biom <- na.omit(biom)
    biom$metadata
#> # A tibble: 40 × 6
#>    .sample   Age   BMI `Body Site`    Sex    BS2           
#>  * <chr>   <dbl> <dbl> <fct>          <fct>  <fct>         
#>  1 HMP01      22    20 Buccal mucosa  Female Buccal mucosa 
#>  2 HMP02      24    23 Buccal mucosa  Male   Buccal mucosa 
#>  3 HMP05      27    24 Buccal mucosa  Female Buccal mucosa 
#>  4 HMP07      26    22 Buccal mucosa  Male   Buccal mucosa 
#>  5 HMP10      22    20 Anterior nares Female Anterior nares
#>  6 HMP11      24    23 Buccal mucosa  Female Buccal mucosa 
#>  7 HMP12      35    26 Stool          Male   Stool         
#>  8 HMP13      24    21 Buccal mucosa  Female Buccal mucosa 
#>  9 HMP14      32    26 Buccal mucosa  Male   Buccal mucosa 
#> 10 HMP15      25    21 Anterior nares Female Anterior nares
#> # ℹ 30 more rows
```
