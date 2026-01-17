# documentation_clusters

documentation_clusters

## Arguments

- k:

  Number of clusters. Default: `5L`

- rank:

  Which taxa rank to use. E.g. `"Phylum"`, `"Genus"`, `".otu"`, etc. An
  integer can also be given, where `1` is the highest rank, `2` is the
  second highest, `-1` is the lowest rank, `-2` is the second lowest,
  and `0` is the OTU "rank". Run `biom$ranks` to see all options for a
  given rbiom object. Default: `.otu`.

## Value

A numeric factor assigning samples to clusters.
