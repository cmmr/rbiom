# Parse a fasta file into a named character vector.

Parse a fasta file into a named character vector.

## Usage

``` r
read_fasta(file, ids = NULL)
```

## Arguments

- file:

  A file path, URL, or literal string with fasta-formatted sequences.
  Files can optionally be compressed with gzip or bzip2.

- ids:

  Character vector of IDs to retrieve. The default, `NULL`, will
  retrieve everything.

## Value

A named character vector in which names are the fasta headers (the text
after `>` up to the first whitespace) and values are the sequences.
Sequences spanning multiple lines are joined, and whitespace within
sequences is removed.

## Examples

``` r
    library(rbiom)

    fasta <- ">OTU1 first sequence\nATCG\nGGCC\n>OTU2\nTTAA\n"
    read_fasta(fasta)
#>       OTU1       OTU2 
#> "ATCGGGCC"     "TTAA" 

    outfile <- write_fasta(hmp50, tempfile(fileext = ".fna.gz"))
    seqs    <- read_fasta(outfile)
    unlink(outfile)

    identical(seqs, hmp50$sequences)
#> [1] TRUE
```
