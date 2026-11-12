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
write_biom(biom, file, format = "json", overwrite = FALSE)

write_metadata(biom, file, quote = FALSE, sep = "\t", overwrite = FALSE, ...)

write_counts(biom, file, quote = FALSE, sep = "\t", overwrite = FALSE, ...)

write_taxonomy(biom, file, quote = FALSE, sep = "\t", overwrite = FALSE, ...)

write_fasta(biom, file = NULL, overwrite = FALSE)

write_tree(biom, file = NULL, overwrite = FALSE)

write_xlsx(
  biom,
  file,
  depth = NULL,
  seed = 0,
  unc = "singly",
  overwrite = FALSE
)
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
  `write_tree()`, and `write_biom(format='json')` returns a string of
  the output which would have been written. For
  `write_biom(format='tab')`, `file=NULL` returns the tibble that would
  have been written.

- format:

  Options are `"tab"`, `"json"`, and `"hdf5"`, corresponding to classic
  tabular format, BIOM format version 1.0 and biom version 2.1,
  respectively. NOTE: to write HDF5 formatted BIOM files, the `h5lite` R
  package must be installed. Default: `"json"`

- overwrite:

  Replace `file` if it already exists? All of the writers refuse to
  overwrite an existing file unless `overwrite = TRUE`. Default: `FALSE`

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

Output is first written to a temporary file in the same directory and
then moved into place, so a failed write never leaves a partial file
behind.

In the tabular (`"tab"`) format, OTU IDs, sample names and taxonomy
strings that contain a tab, a double quote, or a line break are enclosed
in double quotes (with embedded quotes doubled), and missing taxonomic
ranks are written as empty strings. In the JSON format, numbers are
written at full precision and missing values are `null`. In the HDF5
format, the order of metadata fields, logical fields, and `Date` fields
(stored as ISO 8601 strings) are preserved through attributes on the
datasets, which other BIOM readers ignore.

For `write_xlsx()`, `attributes(biom)` are saved as additional
worksheets if the attribute is a data frame, matrix, or dist -class
object. An attribute named 'Reads Per Step' is treated specially and
merged with the usual 'Reads Per Sample' tab.

## Examples

``` r
    library(rbiom)

    write_tree(hmp50) %>% substr(1, 50)
#> [1] "(((((((((((((EschC738:0.03627,(((Unc92490:0.05748,"

    outfile <- write_biom(hmp50, tempfile(fileext = ".biom.gz"))
    file.size(outfile)
#> [1] 72975
    unlink(outfile)

    if (requireNamespace("openxlsx", quietly = TRUE)) {

      hmp10 <- rarefy(hmp50[1:10])

      attr(hmp10, "Weighted UniFrac") <- bdiv_distmat(hmp10, 'weighted_unifrac')
      attr(hmp10, "Jaccard")          <- bdiv_distmat(hmp10, 'jaccard')

      outfile <- write_xlsx(hmp10, tempfile(fileext = ".xlsx"))
      unlink(outfile)
    }
```
