# Changelog

## rbiom 4.0.0

### Reproducibility

- Every rbiom object now records how it was made. `biom$history` lists
  the constructor call and each transformation
  ([`subset()`](https://cmmr.github.io/rbiom/reference/subset.md),
  [`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md),
  [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html), `[`,
  [`biom_merge()`](https://cmmr.github.io/rbiom/reference/biom_merge.md),
  …), and `biom$code` prints them as a runnable script. Direct
  assignments such as `biom$metadata <- df` are recorded as comments,
  since their values exist only in memory.
- `$code` on tables, statistics, plots, matrices, and distance matrices
  is now a complete script: library calls, the object’s history, the
  data step, and the statistics / ggplot2 code. It is covered by tests
  that re-run the code in a fresh environment and compare the result to
  the original, including multi-panel figures, heatmaps, ordinations
  with taxa arrows, and data whose names contain quotes and other R
  syntax.
- Data-derived strings (sample, OTU, taxa, and field names, metadata
  values, file paths) only ever appear in generated code as quoted
  literals, and rbiom never evaluates generated code itself.
- Numbers in generated code parse back to exactly the same values on
  every platform: up to 15 significant digits when that is exact, and
  hexadecimal notation (e.g. `0x1.e999999999999p+3` for `9 * 1.7`)
  otherwise. 17-digit decimals were read back one unit in the last place
  off on macOS (Apple silicon), which changed rank-based statistics.
- `$cmd` no longer re-evaluates the caller’s arguments (so
  `adiv_boxplot(rarefy(biom))` rarefies once) and no longer records
  `cpus`.
- Functions that use random numbers
  ([`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md), rank
  transforms with `ties = "random"`, permutation tests, ordinations,
  clustering) run under `seed` and leave the caller’s RNG state
  untouched.
- [`read_biom()`](https://cmmr.github.io/rbiom/reference/read_biom.md)
  keeps the file’s `generated_by` and creation date.

### Statistics

- Regression statistics with a numeric `regr` no longer scan ~100
  positions along the x-axis for the smallest p-value (which inflated
  false positives to ~17%). `at` defaults to the median of `regr`, and
  p-values are adjusted across the whole table. Without `stat.by`, the
  reported test is the association test (the slope for `lm`/`log`; the
  smooth term for `gam`) rather than a test of “mean = 0”.
- `emtrends` with `stat.by` now evaluates slope differences at each `at`
  position instead of averaging over them.
- Beta diversity statistics are permutation tests on samples.
  [`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md)
  and
  [`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md)
  use PERMANOVA (`adonis2`), MRPP, or a dispersion test (`betadisper`),
  with pairwise PERMANOVA for bracket labels; `regr` and
  [`bdiv_corrplot()`](https://cmmr.github.io/rbiom/reference/bdiv_corrplot.md)
  use a Mantel test. Wilcoxon / t-tests / regressions on pairwise
  distances, which treated dependent distances as independent
  observations, are no longer offered.
- New `strata` argument for repeated-measures designs, on the
  `*_stats()`, `*_boxplot()`, `*_corrplot()`,
  [`bdiv_ord_plot()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_plot.md),
  and
  [`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md)
  functions. Permutation tests permute samples within strata when the
  tested variable varies within strata, and permute whole strata when it
  is constant within each stratum (exact enumeration when there are few
  arrangements). Regression models include a random intercept per
  stratum. Rank tests report an error rather than ignoring `strata`.
- New `covariates` argument for
  [`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md),
  [`adiv_stats()`](https://cmmr.github.io/rbiom/reference/adiv_stats.md),
  and
  [`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md):
  fields to adjust for, added to the regression model as linear terms or
  factors, so that the reported means, slopes, and differences are
  adjusted (emmeans holds numeric covariates at their mean and averages
  over the levels of categorical ones). For `test = "emmeans"` and
  `"emtrends"` only; rank tests report an error. A covariate that is
  constant within a comparison, or confounded with `stat.by` or `regr`,
  is reported as such (an NA row and a warning, or an error when there
  is only one test), rather than as an obscure error or, for `gam` fits,
  a silently dropped coefficient.
- [`bdiv_ord_table()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_table.md)
  /
  [`bdiv_ord_plot()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_plot.md)
  compute distances on the data as given (no hidden relativization), use
  [`stats::cmdscale()`](https://rdrr.io/r/stats/cmdscale.html) for PCoA
  with “% variance explained” axis labels, and test taxa biplot vectors
  with
  [`vegan::envfit()`](https://vegandevs.github.io/vegan/reference/envfit.html).
  `p.top` handles tied p-values and adjusts before filtering.
- [`bdiv_clusters()`](https://cmmr.github.io/rbiom/reference/bdiv_clusters.md)
  uses k-medoids
  ([`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html)) on the
  distance matrix;
  [`taxa_clusters()`](https://cmmr.github.io/rbiom/reference/taxa_clusters.md)
  uses 10 k-means starts. Both accept `seed`.
- Wilcoxon tests use `exact = NULL`; the estimate column is now
  `.loc.shift` (Hodges-Lehmann). `level` is honored. Captions name the
  p-value adjustment method correctly.
- When a Wilcoxon confidence interval can’t reach the requested `level`
  (very small groups), it is reported as unbounded (`-Inf`, `Inf`) on
  every R version, as R \>= 4.6 does, rather than as a narrower interval
  at a lower level. On R \< 4.6, the warnings about this and about
  falling back to the normal approximation for ties are no longer shown.
- [`adiv_boxplot()`](https://cmmr.github.io/rbiom/reference/adiv_boxplot.md),
  [`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md),
  and
  [`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md)
  gain a `test` argument.
- Dropped rows with missing values and comparisons that could not be
  computed are reported in one message / warning instead of silently
  ignored.
- [`bdiv_heatmap()`](https://cmmr.github.io/rbiom/reference/bdiv_heatmap.md)
  clusters samples on the chosen distance metric.
- The log transform in
  [`taxa_matrix()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)
  uses a pseudocount so absent taxa rank lowest.

### Differential abundance

- New [`taxa_da()`](https://cmmr.github.io/rbiom/reference/taxa_da.md)
  finds differentially abundant taxa with methods that account for
  compositionality, between groups (`stat.by`; every pair of groups) or
  along a numeric variable (`regr`), with optional `covariates`,
  `strata` (a random intercept per subject), `split.by`, and several
  ranks. Two methods:
  - `"linda"` (default): LinDA (Zhou et al. 2022), built in. It gives
    the results of `MicrobiomeStat::linda()` (to within 1e-10 on hmp50),
    except that it imputes zeros when sequencing depth is associated
    with the model, as the paper describes (`MicrobiomeStat::linda()`
    1.4 reports imputation but adds a pseudo-count). Comparisons between
    two non-reference groups are corrected by the mode of their own
    contrast, so they don’t depend on the order of the groups.
  - `"ancombc2"`: ANCOM-BC2 (Lin & Peddada 2024) via
    [`ANCOMBC::ancombc2()`](https://rdrr.io/pkg/ANCOMBC/man/ancombc2.html),
    with its pseudo-count sensitivity analysis (`.robust`) and, for
    three or more groups, its mixed directional FDR pairwise tests.

  Both report log2 fold changes with standard errors, confidence
  intervals, t statistics, and p-values adjusted across taxa within each
  comparison. `$code` names the method, its model, and the exact
  `ancombc2()` call.
- New
  [`taxa_forest()`](https://cmmr.github.io/rbiom/reference/taxa_forest.md)
  plots [`taxa_da()`](https://cmmr.github.io/rbiom/reference/taxa_da.md)
  results: each taxon’s log2 fold change and confidence interval,
  colored by direction when significant, with one panel per comparison
  (and rank and `split.by` subset) and the `p.top` most significant taxa
  ordered by fold change.

### Data handling

- [`subset()`](https://cmmr.github.io/rbiom/reference/subset.md) and
  [`subset_taxa()`](https://cmmr.github.io/rbiom/reference/subset.md)
  treat `NA` conditions as `FALSE`, like
  [`base::subset()`](https://rdrr.io/r/base/subset.html).
- Single-sample and single-OTU objects work (`biom[1]`).
- `biom$relative` is `TRUE` for relative-abundance data; `biom$depth` is
  no longer reported as `1` for such data.
- Clear errors for `NA`, negative, or non-numeric counts and unknown
  sample/OTU names.
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md)
  accepts data frames of counts.
- `write_biom(format = "json")` no longer rounds values to four
  decimals.
- Newick output quotes labels containing metacharacters, spaces, or
  underscores, keeps node labels, and handles large trees. With ecodive
  2.3.1’s standard-conforming parser, labels containing apostrophes and
  square brackets (e.g. `o__[Saprospirales]`) now round-trip exactly,
  and malformed Newick input is an error rather than a garbled tree.
- [`read_fasta()`](https://cmmr.github.io/rbiom/reference/read_fasta.md)
  returns a named character vector, accepts URLs and compressed files,
  and is linear-time.
- [`tree_subset()`](https://cmmr.github.io/rbiom/reference/tree_subset.md)
  handles unsorted numeric tips and keeps node labels.
- [`biom_merge()`](https://cmmr.github.io/rbiom/reference/biom_merge.md)
  accepts anything
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md)
  accepts, keeps a tree that covers all OTUs, and reports conflicts in
  one warning.
- Tabular imports no longer depend on readr; sample IDs like `"1.10"`
  stay character, and `Body Site` is never renamed to `Body.Site`.
- Many fixes to TSV/HDF5 round trips (dates, column order, logical
  columns, blank cells),
  [`write_xlsx()`](https://cmmr.github.io/rbiom/reference/write_biom.md),
  `convert_to_*()`, and the `exec/` scripts.
- Caching works when disabled; keys include package versions; only
  rbiom’s own files are ever evicted from the cache directory.

### Plots

- `ggbeeswarm`, `ggrepel`, and `fillpattern` are optional. Without them,
  dot and strip layers use jittered points and taxa labels use
  [`geom_label()`](https://ggplot2.tidyverse.org/reference/geom_text.html).
- Arguments in `...` that don’t match any layer are an error that lists
  the valid layer prefixes (previously they produced cryptic ggplot2
  errors). `b.` now refers to the bar layer and `box.` to the box layer.
- Fixed: flipped boxplots without `x`, confidence intervals clamped at
  zero, heatmap dendrogram placement and `tree_height`, single-taxon
  heatmaps, `rare_stacked(labels = TRUE)`, `xlab.angle`, empty
  statistics tables, and ordination ellipses for small groups.
- Plot titles use readable metric names.

### Figure style

- New
  [`theme_rbiom()`](https://cmmr.github.io/rbiom/reference/theme_rbiom.md):
  one publication-ready theme for every rbiom plot, following common
  journal figure guidelines (Nature, Science, Cell). Black sans-serif
  text sized relative to `base_size` (11 pt; previously 14 pt), open x/y
  axes with outward ticks, no gridlines or background shading, facet
  labels without boxes, and titles and captions aligned with the
  figure’s left edge. A plot’s `$code` now ends with
  [`theme_rbiom()`](https://cmmr.github.io/rbiom/reference/theme_rbiom.md)
  instead of
  [`theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
  plus a block of
  [`theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
  settings. Options: `frame` (`"axes"`, `"box"`, `"none"`), `grid`,
  `base_size`, `base_family`; pass them to any plot function,
  e.g. `adiv_boxplot(biom, base_size = 8)`. Use
  [`theme_rbiom()`](https://cmmr.github.io/rbiom/reference/theme_rbiom.md)
  on your own ggplots to match rbiom’s figures.
- Default colors are now Okabe-Ito for up to 8 groups. The previous
  default (`"classic"`, the D3/Tableau 10 palette) has red/green and
  orange/green pairs that are hard to tell apart with color vision
  deficiencies. The `"okabe"` palette is reordered so that its first
  colors contrast most, with yellow (hard to see on white) last. The
  “Other” group in
  [`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md)
  is light gray.
- Default shapes are solid (circle, triangle, square, diamond) before
  outlines.
- p-values in figures have two significant digits, are never in
  scientific notation (`p = 2e-04`), and are shown as `p < 0.001` below
  0.001.
- Ordination facet labels no longer overflow: statistics are shortened
  to `F = 10.7, R² = 0.63, p = 0.001` (the caption names the test). Taxa
  labels have a fixed size; mean points are scaled by abundance.
- Flipped boxplots draw p-values alongside their brackets instead of
  across them. Heatmap and
  [`rare_stacked()`](https://cmmr.github.io/rbiom/reference/rare_stacked.md)
  sample labels are rotated 90° when there are many of them.

### Breaking changes

- `write_*()` functions refuse to overwrite an existing file unless
  `overwrite = TRUE`. Text output uses LF line endings on all platforms.
- The plyr functions `.`, `ddply`, `ldply`, and `llply` are no longer
  re-exported.
- Removed `transform = "percent"` (never worked), and the
  never-functional `apcoa` and tree-plotting code.
- `bdiv = "unifrac"` is now an error because it is ambiguous; name the
  variant (e.g. `"unweighted_unifrac"`).
- `ord = "tSNE"` uses
  [`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html) instead
  of [`tsne::tsne()`](https://rdrr.io/pkg/tsne/man/tsne.html). It is
  about 8-15x faster from a few hundred samples up and keeps more of
  each sample’s nearest neighbors together. Coordinates differ from
  earlier versions for the same `seed`. Rtsne arguments such as `theta`,
  `max_iter`, `Y_init`, and `num_threads` can be passed through `...`;
  the old tsne arguments (`initial_config`, `whiten`, `epoch`, …) are
  ignored. `perplexity` now defaults to `(n - 1) / 3` for fewer than 91
  samples, since the previous default of 30 cannot be reached with so
  few.
- Statistics table columns changed as described above.

### Dependencies

- Removed: ape, crayon, ggtext, labeling, pillar, plyr, readr. Moved to
  Suggests: fillpattern, ggbeeswarm, ggrepel, readxl. Packages that were
  always required but undeclared (cli, glue, lifecycle, Matrix, R6,
  rlang, scales, tibble) are now in Imports, along with cluster and
  permute (both already installed with vegan).
- Suggests Rtsne instead of tsne.
- Suggests ANCOMBC (for `taxa_da(method = "ancombc2")`), and lmerTest
  and lme4 (for `taxa_da(method = "linda")` with `strata`).
- Requires ecodive (\>= 2.3.1).

## rbiom 3.1.1

- [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md) and
  [`read_biom()`](https://cmmr.github.io/rbiom/reference/read_biom.md)
  now accept count matrices that are triangular, diagonal, or 1x1.
- [`taxa_matrix()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)
  and
  [`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)
  no longer fail with “No taxa match the criteria” when `taxa` is a
  number and a sample has no abundance at the requested rank, e.g. with
  `unc = "drop"`.

## rbiom 3.1.0

CRAN release: 2026-05-08

- Can now read and write `biom` objects from the `biomformat` package
  (microbiome/mia#818).
- Fixed handling of objects from the `TreeSummarizedExperiment` package
  (microbiome/mia#822).

## rbiom 3.0.1

- Re-exports
  [`ecodive::list_metrics`](https://cmmr.github.io/ecodive/reference/metrics.html),
  [`ecodive::match_metric`](https://cmmr.github.io/ecodive/reference/metrics.html),
  and
  [`ggplot2::ggsave`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## rbiom 3.0.0

### Announcements

The code for `rbiom`’s very fast implementations of diversity metrics,
rarefaction, and tree importing has been spun off into its own
dependency-free package, [`ecodive`](https://github.com/cmmr/ecodive),
to encourage adoption by other R packages. This also means `rbiom` will
be easier to install going forward as it will no longer need to be
compiled.

### Breaking Changes

- `weighted` and `normalized` arguments are no longer used by any
  functions.
- Diversity metric names will now be all lowercase in returned matrices,
  tables, etc.

The following functions have been deprecated since `rbiom` v2.0 and are
now removed completely. They are generally functions which have been
renamed for improved organization.

`alpha.div`, `as.percent`, `beta.div`, `comments`, `comments<-`,
`counts`, `counts<-`, `depth`, `depths_barplot`, `has.phylogeny`,
`has.sequences`, `id`, `id<-`, `info`, `is.rarefied`, `metadata`,
`metadata<-`, `nsamples`, `ntaxa`, `phylogeny`, `phylogeny<-`, `rarefy`,
`read.biom`, `read_biom`, `read.fasta`, `read.tree`, `repair`,
`sample.names`, `sample.names<-`, `sample_subset`, `sample.sums`,
`sequences`, `sequences<-`, `subtree`, `taxa.means`, `taxa.names`,
`taxa.names<-`, `taxa.ranks`, `taxa.ranks<-`, `taxa.rollup`,
`taxa.sums`, `taxonomy`, `taxonomy<-`, `tips`, `top.taxa`, `top_taxa`,
`unifrac`, `write.biom`, `write.fasta`, `write.tree`, `write.xlsx`

### Additions

- [`adiv_vector()`](https://cmmr.github.io/rbiom/reference/adiv_table.md)
  returns a simple named vector of alpha diversity values.
- Replaced `rhdf5` dependency with `h5lite`.
- Updated list of alpha diversity metrics:
  `c("ace", "berger", "brillouin", "chao1", "faith", "fisher", "simpson", "inv_simpson", "margalef", "mcintosh", "menhinick", "observed", "shannon", "squares")`
- Updated list of beta diversity metrics:
  `c("aitchison", "bhattacharyya", "bray", "canberra", "chebyshev", "chord", "clark", "sorensen", "divergence", "euclidean", "generalized_unifrac", "gower", "hamming", "hellinger", "horn", "jaccard", "jensen", "jsd", "lorentzian", "manhattan", "matusita", "minkowski", "morisita", "motyka", "normalized_unifrac", "ochiai", "psym_chisq", "soergel", "squared_chisq", "squared_chord", "squared_euclidean", "topsoe", "unweighted_unifrac", "variance_adjusted_unifrac", "wave_hedges", "weighted_unifrac")`

## rbiom 2.2.1

CRAN release: 2025-06-27

- Compatibility fixes for ggplot2 4.0.0
  [\#30](https://github.com/cmmr/rbiom/issues/30)

## rbiom 2.2.0

CRAN release: 2025-04-04

- Satisfy new R CMD check: no setting attributes on built-in functions.
  [\#28](https://github.com/cmmr/rbiom/issues/28)

## rbiom 2.1.2

CRAN release: 2025-03-18

- Fix for `write_biom(format='json')` where `{...}` should be `[...]`.
- Restore `rline` functionality in
  [`rare_corrplot()`](https://cmmr.github.io/rbiom/reference/rare_corrplot.md).

## rbiom 2.1.1

CRAN release: 2025-03-10

- `unifrac()` no longer normalizes weighted values, for
  back-compatibility.
- Option `underscores` allows keeping underscores in tree IDs.

## rbiom 2.0.13

CRAN release: 2025-01-24

- Major release with significant new features.
- Plotting added.
- Statistics added.
- Caching added.
- Clearer naming scheme for functions. Previous names still work but are
  deprecated.
- Generation of human-readable ggplot2 code for modifying plots outside
  of rbiom.
- Provenance tracking for BIOM objects and derivatives.

## rbiom 1.0.3

CRAN release: 2021-11-05

- Corrects for breaking changes in ‘rhdf5’ package.

## rbiom 1.0.2

CRAN release: 2020-05-29

- Improved compatibility with Debian and Solaris.
- The ‘rhdf5’ package is now an optional dependency.
- The select() and subset() functions now subset sequences too.

## rbiom 1.0.0

CRAN release: 2020-05-26

- Initial Release
