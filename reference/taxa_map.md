# Map OTUs names to taxa names at a given rank.

Map OTUs names to taxa names at a given rank.

## Usage

``` r
taxa_map(
  biom,
  rank = NULL,
  taxa = Inf,
  unc = "singly",
  lineage = FALSE,
  other = FALSE
)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- rank:

  When `NULL`, the entire biom\$taxonomy data.frame is returned,
  transformed as per `unc`. Alternatively, a single taxonomic rank (rank
  name or integer position in `biom$ranks`) which returns a named
  character vector for mapping OTUs to taxa names.

- taxa:

  Which taxa to display. An integer value will show the top n most
  abundant taxa. A value 0 \<= n \< 1 will show any taxa with that mean
  abundance or greater (e.g. `0.1` implies \>= 10%). A character vector
  of taxa names will show only those named taxa. Default: `6`.

- unc:

  How to handle unclassified, uncultured, and similarly ambiguous taxa
  names. Options are:

  `"singly"` -

  :   Replaces them with the OTU name.

  `"grouped"` -

  :   Replaces them with a higher rank's name.

  `"drop"` -

  :   Excludes them from the result.

  `"asis"` -

  :   To not check/modify any taxa names.

  Abbreviations are allowed. Default: `"singly"`

- lineage:

  Include all ranks in the name of the taxa. For instance, setting to
  `TRUE` will produce
  `Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales`.
  Otherwise the taxa name will simply be `Coriobacteriales`. You want to
  set this to TRUE when `unc = "asis"` and you have taxa names (such as
  *Incertae_Sedis*) that map to multiple higher level ranks. Default:
  `FALSE`

- other:

  Sum all non-itemized taxa into an "Other" taxa. When `FALSE`, only
  returns taxa matched by the `taxa` argument. Specifying `TRUE` adds
  "Other" to the returned set. A string can also be given to imply
  `TRUE`, but with that value as the name to use instead of "Other".
  Default: `FALSE`

## Value

A tibble data.frame when `rank=NULL`, or a character vector named with
the OTU names.

## See also

[`pull.rbiom()`](https://cmmr.github.io/rbiom/reference/pull.rbiom.md)

## Examples

``` r
    library(rbiom)
    library(dplyr, warn.conflicts = FALSE)
    
    # In $taxonomy, .otu is the first column (like a row identifier)  -----
    hmp50$taxonomy %>% head(4)
#> # A tibble: 4 × 7
#>   .otu     Kingdom  Phylum         Class          Order             Family Genus
#>   <chr>    <fct>    <fct>          <fct>          <fct>             <fct>  <fct>
#> 1 Unc01yki Bacteria Firmicutes     Bacilli        Lactobacillales   Lacto… Lact…
#> 2 Unc53100 Bacteria Firmicutes     Bacilli        Lactobacillales   Strep… Stre…
#> 3 LtbAci52 Bacteria Firmicutes     Bacilli        Lactobacillales   Lacto… Lact…
#> 4 CnbTube3 Bacteria Actinobacteria Actinobacteria Corynebacteriales Coryn… Cory…
    
    # In taxa_map, .otu is the last column (most precise rank)  -----------
    taxa_map(hmp50) %>% head(4)
#> # A tibble: 4 × 7
#>   Kingdom  Phylum         Class          Order             Family    Genus .otu 
#>   <fct>    <fct>          <fct>          <fct>             <fct>     <fct> <chr>
#> 1 Bacteria Firmicutes     Bacilli        Lactobacillales   Lactobac… Lact… Unc0…
#> 2 Bacteria Firmicutes     Bacilli        Lactobacillales   Streptoc… Stre… Unc5…
#> 3 Bacteria Firmicutes     Bacilli        Lactobacillales   Lactobac… Lact… LtbA…
#> 4 Bacteria Actinobacteria Actinobacteria Corynebacteriales Coryneba… Cory… CnbT…
    
    # Generate an OTU to Genus mapping  -----------------------------------
    taxa_map(hmp50, "Genus") %>% head(4)
#>          Unc01yki          Unc53100          LtbAci52          CnbTube3 
#>     Lactobacillus     Streptococcus     Lactobacillus Corynebacterium 1 
#> 251 Levels: Abiotrophia Acidaminococcus Acinetobacter ... Veillonella
    
    # Sometimes taxonomic names are incomplete ----------------------------
    otus <- c('GemAsacc', 'GcbBacte', 'Unc58411')
    taxa_map(hmp50, unc = "asis") %>% filter(.otu %in% otus) %>% select(Phylum:.otu)
#> # A tibble: 3 × 6
#>   Phylum          Class      Order          Family    Genus   .otu    
#>   <fct>           <fct>      <fct>          <fct>     <fct>   <chr>   
#> 1 Gracilibacteria c          o              f         g       GcbBacte
#> 2 Tenericutes     Mollicutes Mollicutes RF9 f         g       Unc58411
#> 3 Firmicutes      Bacilli    Bacillales     Family XI Gemella GemAsacc
    
    # rbiom can replace them with unique placeholders ---------------------
    taxa_map(hmp50, unc = "singly") %>% filter(.otu %in% otus) %>% select(Class:.otu)
#> # A tibble: 3 × 5
#>   Class         Order          Family        Genus         .otu    
#>   <fct>         <fct>          <fct>         <fct>         <chr>   
#> 1 Unc. GcbBacte Unc. GcbBacte  Unc. GcbBacte Unc. GcbBacte GcbBacte
#> 2 Mollicutes    Mollicutes RF9 Unc. Unc58411 Unc. Unc58411 Unc58411
#> 3 Bacilli       Bacillales     Bacillales XI Gemella       GemAsacc
    
    # Or collapse them into groups ----------------------------------------
    taxa_map(hmp50, unc = "grouped") %>% filter(.otu %in% otus) %>% select(Class:Genus)
#> # A tibble: 3 × 4
#>   Class                Order                Family               Genus          
#>   <fct>                <fct>                <fct>                <fct>          
#> 1 Unc. Gracilibacteria Unc. Gracilibacteria Unc. Gracilibacteria Unc. Graciliba…
#> 2 Mollicutes           Mollicutes RF9       Unc. Mollicutes RF9  Unc. Mollicute…
#> 3 Bacilli              Bacillales           Bacillales XI        Gemella        
```
