# Working with rbiom Objects.

Rbiom objects make it easy to access and manipulate your BIOM data,
ensuring all the disparate components remain in sync. These objects
behave largely like lists, in that you can access and assign to them
using the `$` operator. The sections below list all the fields which can
be read and/or written, and the helper functions for common tasks like
rarefying and subsetting. To create an rbiom object, see
[`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

Use `$clone()` to create a copy of an rbiom object. This is necessary
because rbiom objects are **passed by reference**. The usual `<-`
assignment operator will simply create a second reference to the same
object - it will not create a second object. See [speed
ups](https://cmmr.github.io/rbiom/reference/speed.md) for more details.

## Readable Fields

Reading from fields will not change the rbiom object.

|                          |                                                        |
|--------------------------|--------------------------------------------------------|
| **Accessor**             | **Content**                                            |
| `$counts`                | Abundance of each OTU in each sample.                  |
| `$metadata`              | Sample mappings to metadata (treatment, patient, etc). |
| `$taxonomy`              | OTU mappings to taxonomic ranks (genus, phylum, etc).  |
| `$otus`, `$n_otus`       | OTU names.                                             |
| `$samples`, `$n_samples` | Sample names.                                          |
| `$fields`, `$n_fields`   | Metadata field names.                                  |
| `$ranks`, `$n_ranks`     | Taxonomic rank names.                                  |
| `$tree`, `$sequences`    | Phylogenetic tree / sequences for the OTUs, or `NULL`. |
| `$id`, `$comment`        | Arbitrary strings for describing the dataset.          |
| `$depth`                 | Rarefaction depth, or `NULL` if unrarefied.            |
| `$date`                  | Date from BIOM file.                                   |

## Writable Fields

Assigning new values to these components will trigger validation checks
and inter-component synchronization.

|                   |                                                         |
|-------------------|---------------------------------------------------------|
| **Component**     | **What can be assigned.**                               |
| `$counts`         | Matrix of abundances; OTUs (rows) by samples (columns). |
| `$metadata`       | Data.frame with `'.sample'` column, or a file name.     |
| `$taxonomy`       | Data.frame with `'.otu'` as the first column.           |
| `$otus`           | Character vector with new names for the OTUs.           |
| `$samples`        | Character vector with new names for the samples.        |
| `$tree`           | Phylo object with the phylogenetic tree for the OTUs.   |
| `$sequences`      | Named character vector of OTU reference sequences.      |
| `$id`, `$comment` | String with dataset's title or comment.                 |
| `$date`           | Date-like object, or `"%Y-%m-%dT%H:%M:%SZ"` string.     |

## Transformations

All functions return an rbiom object.

|                                                                 |                                                  |
|-----------------------------------------------------------------|--------------------------------------------------|
| **Function**                                                    | **Transformation**                               |
| `<rbiom>$clone()`                                               | Safely duplicate an rbiom object.                |
| [`<rbiom>[`](https://cmmr.github.io/rbiom/reference/subset.md)  | Subset to a specific set of sample names.        |
| [`subset()`](https://cmmr.github.io/rbiom/reference/subset.md)  | Subset samples according to metadata properties. |
| [`slice()`](https://dplyr.tidyverse.org/reference/slice.html)   | Subset to a specific number of samples.          |
| [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) | Create, modify, and delete metadata fields.      |
| [`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md)  | Sub-sample OTU counts to an even sampling depth. |

## Examples

``` r
    library(rbiom)
    
    # Duplicate the HMP50 example dataset.
    biom <- hmp50$clone()
    
    
    # Display an overall summary of the rbiom object.
    biom
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
    
    
    # Markdown syntax for comments is recommended.
    biom$comment %>% cli::cli_text()
#> Oral, nasal, vaginal, and fecal samples from a diverse set of healthy
#> volunteers. Source: [Human Microbiome Project](https://hmpdacc.org).
    
    
    # Demonstrate a few accessors.
    biom$n_samples
#> [1] 50
    biom$fields
#> [1] ".sample"   "Age"       "BMI"       "Body Site" "Sex"      
    biom$metadata
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
    
    
    # Edit the metadata table.
    biom$metadata$rand <- sample(1:50)
    biom %<>% mutate(Obese = BMI >= 30, Sex = NULL)
    biom %<>% rename('Years Old' = "Age")
    biom$metadata
#> # A tibble: 50 × 6
#>    .sample `Years Old`   BMI `Body Site`     rand Obese
#>  * <chr>         <dbl> <dbl> <fct>          <int> <lgl>
#>  1 HMP01            22    20 Buccal mucosa     31 FALSE
#>  2 HMP02            24    23 Buccal mucosa     16 FALSE
#>  3 HMP03            28    26 Saliva            30 FALSE
#>  4 HMP04            25    23 Saliva             6 FALSE
#>  5 HMP05            27    24 Buccal mucosa     43 FALSE
#>  6 HMP06            32    25 Saliva             8 FALSE
#>  7 HMP07            26    22 Buccal mucosa     22 FALSE
#>  8 HMP08            27    26 Saliva            44 FALSE
#>  9 HMP09            33    32 Saliva            39 TRUE 
#> 10 HMP10            22    20 Anterior nares    50 FALSE
#> # ℹ 40 more rows
    
    
    # Subset the rbiom object
    biom %<>% subset(`Body Site` == "Saliva" & !Obese)
    biom$metadata
#> # A tibble: 8 × 6
#>   .sample `Years Old`   BMI `Body Site`  rand Obese
#> * <chr>         <dbl> <dbl> <fct>       <int> <lgl>
#> 1 HMP03            28    26 Saliva         30 FALSE
#> 2 HMP04            25    23 Saliva          6 FALSE
#> 3 HMP06            32    25 Saliva          8 FALSE
#> 4 HMP08            27    26 Saliva         44 FALSE
#> 5 HMP18            28    24 Saliva         35 FALSE
#> 6 HMP28            23    19 Saliva         28 FALSE
#> 7 HMP29            36    25 Saliva         24 FALSE
#> 8 HMP30            24    21 Saliva          7 FALSE
    
    
    # Rarefy to an even sampling depth
    sample_sums(biom)
#> HMP03 HMP04 HMP06 HMP08 HMP18 HMP28 HMP29 HMP30 
#>  1353  1895  4150  1695  2202  1695  2423  3938 
    
    biom %<>% rarefy()
    sample_sums(biom)
#> HMP03 HMP04 HMP06 HMP08 HMP18 HMP28 HMP29 HMP30 
#>  1353  1353  1353  1353  1353  1353  1353  1353 
```
