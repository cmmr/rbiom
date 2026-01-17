# documentation_rank.NULL

documentation_rank.NULL

## Arguments

- rank:

  What rank(s) of taxa to compute biplot coordinates and statistics for,
  or `NULL` to disable. E.g. `"Phylum"`, `"Genus"`, `".otu"`, etc. An
  integer vector can also be given, where `1` is the highest rank, `2`
  is the second highest, `-1` is the lowest rank, `-2` is the second
  lowest, and `0` is the OTU "rank". Run `biom$ranks` to see all options
  for a given rbiom object. Default: `NULL`.
