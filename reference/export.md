# Export data to QIIME 2 or mothur.

Populates a directory with the following files, formatted according to
QIIME 2 or mothur's specifications.

- `biom_counts.tsv`

- `biom_metadata.tsv`

- `biom_taxonomy.tsv`

- `biom_tree.nwk`

- `biom_seqs.fna`

`biom_counts.tsv` will always be created. The others are dependent on
whether the content is present in the `biom` argument.

## Usage

``` r
write_mothur(biom, dir = tempfile(), prefix = "biom_")

write_qiime2(biom, dir = tempfile(), prefix = "biom_")
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

- dir:

  Where to save the files. If the directory doesn't exist, it will be
  created. Default: [`tempfile()`](https://rdrr.io/r/base/tempfile.html)

- prefix:

  A string to prepend to each file name. Default: `'biom_'`

## Value

The normalized directory path that was written to (invisibly).

## Examples

``` r
    library(rbiom)
    
    tdir <- tempfile()
    
    write_qiime2(hmp50, tdir, 'qiime2_')
    write_mothur(hmp50, tdir, 'mothur_')
    
    list.files(tdir)
#>  [1] "mothur_counts.tsv"   "mothur_metadata.tsv" "mothur_seqs.fna"    
#>  [4] "mothur_taxonomy.tsv" "mothur_tree.nwk"     "qiime2_counts.tsv"  
#>  [7] "qiime2_metadata.tsv" "qiime2_seqs.fna"     "qiime2_taxonomy.tsv"
#> [10] "qiime2_tree.nwk"    
    
    readLines(file.path(tdir, 'qiime2_metadata.tsv'), n = 4)
#> [1] "sample-id\tAge\tBMI\t\"Body Site\"\tSex"              
#> [2] "#q2:types\tnumeric\tnumeric\tcategorical\tcategorical"
#> [3] "HMP01\t22\t20\t\"Buccal mucosa\"\tFemale"             
#> [4] "HMP02\t24\t23\t\"Buccal mucosa\"\tMale"               
    
    readLines(file.path(tdir, 'mothur_taxonomy.tsv'), n = 3)
#> [1] "OTU\tSize\tTaxonomy"                                                                 
#> [2] "1\t24096\tBacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus"
#> [3] "2\t23409\tBacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus"
    
    unlink(tdir, recursive = TRUE)
```
