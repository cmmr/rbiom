# Parse a fasta file into a named character vector.

Parse a fasta file into a named character vector.

## Usage

``` r
read_fasta(file, ids = NULL)
```

## Arguments

- file:

  A file/URL with fasta-formatted sequences. Can optionally be
  compressed with gzip, bzip2, xz, or lzma.

- ids:

  Character vector of IDs to retrieve. The default, `NULL`, will
  retrieve everything.

## Value

A named character vector in which names are the fasta headers and values
are the sequences.
