# documentation_dist_test

documentation_dist_test

## Arguments

- stat.by:

  The categorical or numeric metadata field over which statistics should
  be calculated. Required.

- test:

  Permutational test for accessing significance. Options are:

  `"adonis2"` -

  :   Permutational MANOVA;
      [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html).

  `"mrpp"` -

  :   Multiple response permutation procedure;
      [`vegan::mrpp()`](https://vegandevs.github.io/vegan/reference/mrpp.html).

  `"none"` -

  :   Don't run any statistics.

  Abbreviations are allowed. Default: `"adonis2"`
