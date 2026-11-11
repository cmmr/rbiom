# documentation_dist_test

documentation_dist_test

## Arguments

- stat.by:

  The categorical or numeric metadata field over which statistics should
  be calculated. Required.

- test:

  Permutational test for assessing significance. Options are:

  `"adonis2"` -

  :   Permutational MANOVA;
      [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html).

  `"mrpp"` -

  :   Multiple response permutation procedure;
      [`vegan::mrpp()`](https://vegandevs.github.io/vegan/reference/mrpp.html).

  `"betadisper"` -

  :   Homogeneity of group dispersions;
      [`vegan::betadisper()`](https://vegandevs.github.io/vegan/reference/betadisper.html)
      with
      [`vegan::permutest()`](https://vegandevs.github.io/vegan/reference/anova.cca.html).

  `"mantel"` -

  :   Mantel test of the distance matrix against a numeric variable;
      [`vegan::mantel()`](https://vegandevs.github.io/vegan/reference/mantel.html).

  `"none"` -

  :   Don't run any statistics.

  Abbreviations are allowed. Default: `"adonis2"` for categorical groups
  and `"mantel"` for numeric groups.

- strata:

  Dataset field (or vector) naming the blocking / subject variable for
  repeated-measures designs. Permutations are restricted to within
  strata via
  [`permute::how()`](https://gavinsimpson.github.io/permute/reference/how.html),
  and regression models include a random intercept per stratum. Default:
  `NULL`
