# Save an rbiom object to a file.

Automatically creates directories and adds compression based on file
name.

- `write_biom()` - :

  According to [BIOM format](http://biom-format.org/documentation/)
  specification.

- `write_xlsx()` - :

  Raw data and summary tables in Excel file format. See details.

- `write_fasta()` - :

  Sequences only in fasta format. `biom` may also be a named character
  vector.

- `write_tree()` - :

  Phylogenetic tree only in newick format. `biom` may also be a phylo
  object.

- `write_counts()`, `write_metadata()`, `write_taxonomy()` - :

  Tab-separated values.

## Usage

``` r
write_biom(biom, file, format = "json")

write_metadata(biom, file, quote = FALSE, sep = "\t", ...)

write_counts(biom, file, quote = FALSE, sep = "\t", ...)

write_taxonomy(biom, file, quote = FALSE, sep = "\t", ...)

write_fasta(biom, file = NULL)

write_tree(biom, file = NULL)

write_xlsx(biom, file, depth = NULL, seed = 0, unc = "singly")
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- file:

  Path to the output file. File names ending in `.gz` or `.bz2` will be
  compressed accordingly. Setting `file=NULL` for `write_fasta()`,
  `write_tree()`, and `write_biom(format='json')`, and returns a string
  of the output which would have been written. For
  `write_biom(format='tab')`, `file=NULL` returns the tibble that would
  have been written.

- format:

  Options are `"tab"`, `"json"`, and `"hdf5"`, corresponding to classic
  tabular format, BIOM format version 1.0 and biom version 2.1,
  respectively. NOTE: to write HDF5 formatted BIOM files, the `h5lite` R
  package must be installed. Default: `"json"`

- quote, sep, ...:

  Parameters passed on to
  [`write.table()`](https://rdrr.io/r/utils/write.table.html). Default:
  `quote=FALSE, sep="\t"`

- depth:

  Passed on to
  [`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md). For
  `write_xlsx()` only, `depth=0` disables rarefaction. Default: `NULL`

- seed:

  Random seed to use in rarefying. See
  [`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md)
  function for details. Must be a non-negative integer. Default: `0`

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

## Value

The normalized filepath that was written to (invisibly), unless
`file=NULL` (see `file` argument above).

## Details

For `write_xlsx()`, `attributes(biom)` are saved as additional
worksheets if the attribute is a data frame, matrix, or dist -class
object. An attribute named 'Reads Per Step' is treated specially and
merged with the usual 'Reads Per Sample' tab.

## Examples

``` r
    library(rbiom)
    
    write_tree(hmp50) %>% substr(1, 50)
#> [1] "(((((((((((((EschC738:0.03627,(((Unc92490:0.05748,"
    
    if (FALSE) {
    
      hmp10        <- hmp50$clone()
      hmp10$counts <- hmp10$counts[,1:10] %>% rarefy_cols()
      
      attr(hmp10, "Weighted UniFrac") <- bdiv_distmat(hmp10, 'wunifrac')
      attr(hmp10, "Jaccard")          <- bdiv_distmat(hmp10, 'jaccard')
      
      outfile <- write_xlsx(hmp10, tempfile(fileext = ".xlsx"))
    }
```
