# rbiom 4.0.0

## Reproducibility

* `$methods` on any table, statistics table, plot, matrix, distance matrix, 
  or rbiom object returns a paragraph for a manuscript's methods section: 
  data preparation, metrics, ordination, tests, permutations, p-value 
  adjustment, and software versions, followed by the references it cites. 
  `fn_schema()` lists `"methods"` in `returns$has`.
* Every rbiom object now records how it was made. `biom$history` lists the
  constructor call and each transformation (`subset()`, `rarefy()`, 
  `mutate()`, `[`, `biom_merge()`, ...), and `biom$code` prints them as a 
  runnable script. Direct assignments such as `biom$metadata <- df` are 
  recorded as comments, since their values exist only in memory.
* `$code` on tables, statistics, plots, matrices, and distance matrices is now
  a complete script: library calls, the object's history, the data step, and 
  the statistics / ggplot2 code. It is covered by tests that re-run the code 
  in a fresh environment and compare the result to the original, including 
  multi-panel figures, heatmaps, ordinations with taxa arrows, and data whose 
  names contain quotes and other R syntax.
* Data-derived strings (sample, OTU, taxa, and field names, metadata values, 
  file paths) only ever appear in generated code as quoted literals, and rbiom
  never evaluates generated code itself.
* Numbers in generated code are decimals with the fewest digits (up to 17)
  that read back as the same value, e.g. `15.299999999999999` for `9 * 1.7`;
  never hexadecimal. Plot positions that rbiom computes (p-value labels, axis
  breaks, dendrogram segments) are rounded so they read well. R 4.3 on macOS
  (Apple silicon) can read 16-17 digit literals back one unit in the last 
  place off; newer R versions on the same machines replay them exactly.
* `$cmd` no longer re-evaluates the caller's arguments (so 
  `adiv_boxplot(rarefy(biom))` rarefies once) and no longer records `cpus`.
* Functions that use random numbers (`rarefy()`, rank transforms with 
  `ties = "random"`, permutation tests, ordinations, clustering) run under 
  `seed` and leave the caller's RNG state untouched.
* `read_biom()` keeps the file's `generated_by` and creation date.

## Statistics

* Regression statistics with a numeric `regr` no longer scan ~100 positions 
  along the x-axis for the smallest p-value (which inflated false positives to
  ~17%). `at` defaults to the median of `regr`, and p-values are adjusted 
  across the whole table. Without `stat.by`, the reported test is the 
  association test (the slope for `lm`/`log`; the smooth term for `gam`) 
  rather than a test of "mean = 0".
* `emtrends` with `stat.by` now evaluates slope differences at each `at` 
  position instead of averaging over them.
* Beta diversity statistics are permutation tests on samples. `bdiv_stats()`
  and `bdiv_boxplot()` use PERMANOVA (`adonis2`), MRPP, or a dispersion test 
  (`betadisper`), with pairwise PERMANOVA for bracket labels; `regr` and 
  `bdiv_corrplot()` use a Mantel test. Wilcoxon / t-tests / regressions on 
  pairwise distances, which treated dependent distances as independent 
  observations, are no longer offered.
* `within` and `between` are deprecated in `bdiv_stats()`. They select 
  sample pairs, and the tests now run on samples, so they have no effect; 
  giving either is a warning. Use `split.by` to test within each group of a
  field, or `test = "betadisper"` to compare how spread out the groups are.
  They still select the pairs shown by `bdiv_table()` and the `bdiv_*` 
  plots.
* New `strata` argument for repeated-measures designs, on the `*_stats()`, 
  `*_boxplot()`, `*_corrplot()`, `bdiv_ord_plot()`, and `distmat_stats()` 
  functions. Permutation tests permute samples within strata when the tested 
  variable varies within strata, and permute whole strata when it is constant 
  within each stratum (exact enumeration when there are few arrangements). 
  Regression models include a random intercept per stratum. Rank tests report
  an error rather than ignoring `strata`.
* New `covariates` argument for `stats_table()`, `adiv_stats()`, and 
  `taxa_stats()`: fields to adjust for, added to the regression model as 
  linear terms or factors, so that the reported means, slopes, and 
  differences are adjusted (emmeans holds numeric covariates at their mean 
  and averages over the levels of categorical ones). For `test = "emmeans"` 
  and `"emtrends"` only; rank tests report an error. A covariate that is 
  constant within a comparison, or confounded with `stat.by` or `regr`, is 
  reported as such (an NA row and a warning, or an error when there is only 
  one test), rather than as an obscure error or, for `gam` fits, a silently
  dropped coefficient.
* `bdiv_ord_table()` / `bdiv_ord_plot()` compute distances on the data as 
  given (no hidden relativization), use `stats::cmdscale()` for PCoA with 
  "% variance explained" axis labels, and test taxa biplot vectors with 
  `vegan::envfit()`. `p.top` handles tied p-values and adjusts before 
  filtering.
* `bdiv_clusters()` uses k-medoids (`cluster::pam()`) on the distance matrix;
  `taxa_clusters()` uses 10 k-means starts. Both accept `seed`.
* Wilcoxon tests use `exact = NULL`; the estimate column is now `.loc.shift` 
  (Hodges-Lehmann). `level` is honored. Captions name the p-value adjustment 
  method correctly.
* When a Wilcoxon confidence interval can't reach the requested `level` (very
  small groups), it is reported as unbounded (`-Inf`, `Inf`) on every R 
  version, as R >= 4.6 does, rather than as a narrower interval at a lower 
  level. On R < 4.6, the warnings about this and about falling back to the 
  normal approximation for ties are no longer shown.
* `adiv_boxplot()`, `bdiv_boxplot()`, and `taxa_boxplot()` gain a `test` 
  argument.
* `adiv_boxplot()`, `adiv_corrplot()`, `adiv_stats()`, and `write_xlsx()` 
  gain a `cpus` argument, like their beta diversity counterparts.
* `cpus` also bounds R's linear algebra library (BLAS). OpenBLAS, MKL, and 
  BLIS start up to one thread per core for operations such as the 
  eigendecomposition in a PCoA, whatever `cpus` says: with `cpus = 1` and 
  OpenBLAS on six cores, a PCoA of 4,020 samples used 55 s of processor 
  time in 15.5 s on an idle machine, and up to 148 s in 31 s on a busy 
  one. When `cpus` is less than the library's thread count, rbiom now runs
  it on one thread until the call returns (that PCoA: 17.5 s, all of it on
  one core) and then restores the setting, also after an error or an 
  interrupt. More BLAS threads were never more than about 20% faster in 
  rbiom's computations, and two were five times slower than one. The 
  count is read and set with `ecodive::blas_threads()`, new in ecodive 
  2.4.0; where it can't be read (R's reference BLAS, Apple's Accelerate), 
  nothing changes.
* `taxa_stats()`, `taxa_boxplot()`, `taxa_corrplot()`, `taxa_da()`, 
  `taxa_forest()`, `stats_table()`, `stats_boxplot()`, `stats_corrplot()`, 
  `distmat_stats()`, and `distmat_ord_table()` gain a `cpus` argument, 
  which bounds the linear algebra threads of their model fits, permutation 
  tests, and ordinations.
* Dropped rows with missing values and comparisons that could not be computed
  are reported in one message / warning instead of silently ignored.
* `bdiv_heatmap()` clusters samples on the chosen distance metric.
* The log transform in `taxa_matrix()` uses a pseudocount so absent taxa rank 
  lowest.

## Differential abundance

* New `taxa_da()` finds differentially abundant taxa with LinDA (Zhou et al.
  2022), which accounts for compositionality, between groups (`stat.by`; 
  every pair of groups) or along a numeric variable (`regr`), with optional
  `covariates`, `strata` (a random intercept per subject), `split.by`, and 
  several ranks. LinDA is built in. It gives the results of 
  `MicrobiomeStat::linda()` (to within 1e-10 on hmp50), except that it 
  imputes zeros when sequencing depth is associated with the model, as the 
  paper describes (`MicrobiomeStat::linda()` 1.4 reports imputation but adds
  a pseudo-count). Comparisons between two non-reference groups are 
  corrected by the mode of their own contrast, so they don't depend on the 
  order of the groups.
  
  It reports log2 fold changes with standard errors, confidence intervals,
  t statistics, and p-values adjusted across taxa within each comparison. 
  `$code` names the method and its model.
* New `taxa_forest()` plots `taxa_da()` results: each taxon's log2 fold 
  change and confidence interval, colored by direction when significant, 
  with one panel per comparison (and rank and `split.by` subset) and the 
  `p.top` most significant taxa ordered by fold change.

## Data handling

* `subset()` and `subset_taxa()` treat `NA` conditions as `FALSE`, like 
  `base::subset()`.
* Single-sample and single-OTU objects work (`biom[1]`).
* `biom$relative` is `TRUE` for relative-abundance data; `biom$depth` is no 
  longer reported as `1` for such data.
* Clear errors for `NA`, negative, or non-numeric counts and unknown 
  sample/OTU names. `as_rbiom()` accepts data frames of counts.
* `write_biom(format = "json")` no longer rounds values to four decimals.
* `biom_inflate()` and `suggest_inflate_depths()` work on values of any 
  scale: they previously assumed each sample summed to 1, so counts divided 
  by 4 were inflated to depths of about 3. A sample whose values are all 
  whole multiples of its smallest value now gets its original depth exactly 
  (counts divided by a constant, or proportions of counts with a singleton); 
  the singleton-peak heuristic, and its documented accuracy, covers the rest.
* Newick output quotes labels containing metacharacters, spaces, or 
  underscores, keeps node labels, and handles large trees. With ecodive 
  2.3.1's standard-conforming parser, labels containing apostrophes and 
  square brackets (e.g. `o__[Saprospirales]`) now round-trip exactly, and 
  malformed Newick input is an error rather than a garbled tree.
* `read_fasta()` returns a named character vector, accepts URLs and 
  compressed files, and is linear-time.
* `tree_subset()` handles unsorted numeric tips and keeps node labels.
* `biom_merge()` accepts anything `as_rbiom()` accepts, keeps a tree that 
  covers all OTUs, and reports conflicts in one warning.
* Tabular imports no longer depend on readr; sample IDs like `"1.10"` stay 
  character, and `Body Site` is never renamed to `Body.Site`.
* Many fixes to TSV/HDF5 round trips (dates, column order, logical columns, 
  blank cells), `write_xlsx()`, `convert_to_*()`, and the `exec/` scripts.
* Caching works when disabled; keys include package versions; only rbiom's 
  own files are ever evicted from the cache directory.
* Several R processes can share one cache directory (e.g. a pool of workers
  with the same `rbiom.cache_dir`). Entries are written to a temporary file 
  and renamed into place, so no process ever reads a partly written entry. 
  The directory is listed now and then rather than on every write, and the 
  least recently used entries are evicted down to 80% of `rbiom.cache_size`.
  Temporary files left by a process that died are removed after an hour. A 
  cache that can't be read or written never causes an error or a warning.
* The default cache size is now 50 MB (previously 200 MB, although the 
  documentation said 20 MB).
* Distance matrices that differ only in provenance share cache entries in 
  `distmat_ord_table()` and `distmat_stats()`, as do `tree_subset()` calls 
  that differ only in `underscores`.

## Importing files

* QIIME 2 artifacts are read directly: `read_biom()` and `as_rbiom()` take a 
  feature table `.qza`, `taxonomy =` a taxonomy `.qza`, `tree =` and 
  `read_tree()` a phylogeny `.qza`, and `read_fasta()` and `sequences =` a 
  sequence `.qza`. Visualizations (`.qzv`), other artifact types, and an 
  artifact given to the wrong reader are errors that name what the file is.
* `biom$sequences <- ` and `as_rbiom(sequences = )` accept a fasta file, URL, 
  or fasta string as well as a named character vector.
* `read_biom()` reads mothur `.shared` files. mothur taxonomy files have 
  their `(100)` confidence scores and `Size` column removed, and lineages 
  separated by commas are split into ranks.
* `write_biom(format = "hdf5")` records the order of factor levels, and 
  `read_biom()` restores it. A factor whose labels look like numbers is now 
  read back as a factor. The JSON and tabular formats and `write_metadata()` 
  cannot record the order.
* `write_biom(format = "hdf5")` writes the tree in the form that the Python 
  biom-format package expects; it previously rejected rbiom's HDF5 files that 
  had a tree.
* Import errors name the file and the problem: decimal commas in counts (an 
  error) or metadata (a warning; the field stays categorical), Excel or 
  binary files given as counts, UTF-16 text, NEXUS and other non-Newick tree 
  files, and trees whose tips don't match the OTUs (with counts and examples 
  from both sides). Sample or OTU names that differ only by case, 
  whitespace, or punctuation are pointed out.
* Text files in Latin-1 encoding are read instead of failing.

## Shared and dominant taxa

* `taxa_top()` lists the most abundant taxa overall, per metadata group, or 
  per sample.
* `taxa_shared()` reports which taxa each group has, at a chosen detection 
  level and prevalence, and the size of every combination of groups. With no 
  grouping it gives the core taxa of the whole dataset.
* `taxa_upset()` draws those combinations as an UpSet plot.

## Describing and checking arguments

* New `fn_schema()`: a machine-readable description of the arguments of the
  plot functions, `adiv_stats()`, `bdiv_stats()`, `taxa_stats()`, 
  `distmat_stats()`, `taxa_da()`, `taxa_clusters()`, `bdiv_clusters()`, 
  `rarefy()`, `adiv_matrix()`, `bdiv_distmat()`, `taxa_matrix()`, 
  `adiv_table()`, `bdiv_table()`, `taxa_table()`, `bdiv_ord_table()`, 
  `distmat_ord_table()`, `write_biom()`, and `write_xlsx()`. For each 
  argument: its kind, default, 
  fixed choices (with names for layers and metrics), numeric limits, how 
  many values it takes, what `NULL` means, where dataset-dependent choices 
  come from (a categorical or numeric field, a rank, taxa, a metric that 
  may need a tree), the forms it accepts, when it has any effect, and the 
  first sentence of its help text; and the constraints between arguments 
  (e.g. `test = "kruskal"` needs `stat.by`). Plain lists, or JSON with 
  `json = TRUE`, for programs that build input controls or check requests.
  Entries also say which counts a function wants (`counts`: `"rarefied"` 
  or `"unrarefied"`) and which optional package it needs (`needs`, e.g. 
  `write_xlsx()`); choices that need an optional package are listed in 
  `choices_needs` (`ord`'s `"tSNE"` and `"UMAP"`, `bdiv_ord_plot()`'s 
  `"density"` layer, `write_biom()`'s `"hdf5"` format).
* New `fn_check()` and `fn_check_data()`: check arguments against a 
  dataset, or a summary of one, without running the function. Errors are 
  the function's own. Deprecated arguments (`weighted`, `normalized`, and 
  `bdiv_stats()`'s `within` and `between`) are refused, with what to use 
  instead.
* New `list_palettes()`: palette names, kinds (discrete, sequential, 
  diverging), and colors.
* The covered functions check their arguments against the same description 
  on entry, and the validators take their option lists and limits from it,
  so the description and the behavior can't disagree.
* Statistical tests that don't apply to the plotted values are refused 
  instead of being ignored: `adonis2`, `mrpp`, and `betadisper` on alpha 
  diversity and taxa boxplots (and `stats_boxplot()`), which ran 
  Kruskal-Wallis tests or failed later; `wilcox`, `emmeans`, and `kruskal` 
  on `bdiv_boxplot()`, which failed later; and `mantel` on alpha diversity,
  taxa, and rarefaction corrplots (and `stats_corrplot()`), which skipped 
  the statistics with a warning.
* `clust` accepts `"ward.D"` and `"ward.D2"`, as documented (both were 
  refused). `"ward"` still means `"ward.D"`, without `hclust()`'s message 
  about the renamed method.
* `taxa_stacked()`'s `xlab.angle` default is `"90"`, a string like the 
  other choices (the number 90 was, and is, accepted too).
* Heatmap `label`, `trees`, `rescale`, `legend`, `xlab.angle`, 
  `label_size`, `tree_height`, and `track_height` errors name the argument
  and list the choices, instead of quoting `stopifnot()` code. Like other 
  choices, they may be abbreviated.
* A phylogenetic metric (`faith`, the UniFrac family) without a tree is an 
  error naming the metric, before any computing (previously 
  "inherits(tree, "phylo") is not TRUE").
* `adiv_stats()` and `taxa_stats()` report a test's missing or unsupported 
  arguments (`stat.by` for `"kruskal"`, `regr` for `"emtrends"`, `strata` 
  or `covariates` with a rank test) before computing anything.
* A misspelled field name is reported under the argument it was given in 
  (e.g. `x`) rather than `md`. `k` in `taxa_clusters()` and 
  `bdiv_clusters()` must be a whole number of at least 1.
* A metric name that doesn't match is reported under the caller's argument
  (`adiv = "nosuch"`, `bdiv = "nosuch"`), not ecodive's `metric`.
* Beta diversity heat maps, ordinations, clusterings, box plots, scatter 
  plots, and statistics stop with a clear message when the metric gives 
  infinite or missing distances (e.g. `"bhattacharyya"` between samples with
  no taxa in common), rather than an internal error.

## Plots

* `ggbeeswarm`, `ggrepel`, and `fillpattern` are optional. Without them, dot 
  and strip layers use jittered points and taxa labels use `geom_label()`.
* Arguments in `...` that don't match any layer are an error that lists the 
  valid layer prefixes (previously they produced cryptic ggplot2 errors). 
  `b.` now refers to the bar layer and `box.` to the box layer.
* Fixed: flipped boxplots without `x`, confidence intervals clamped at zero, 
  heatmap dendrogram placement and `tree_height`, single-taxon heatmaps, 
  `rare_stacked(labels = TRUE)`, `xlab.angle`, empty statistics tables, and 
  ordination ellipses for small groups.
* Plot titles use readable metric names.
* Every plot's caption states the input table (rarefied to a depth, 
  unrarefied, or relative abundances) and the method: the metric and its 
  parameters, the test, the p-value adjustment, and the number of 
  permutations. Heatmaps, stacked bar plots, and rarefaction plots have 
  captions and a `caption` argument, and `bdiv_ord_plot()` has a caption 
  without `stat.by`.
* NMDS stress is reported: in the `.stress` column of the `axes` attribute of 
  `distmat_ord_table()` and `bdiv_ord_table()`, and in the caption of 
  `bdiv_ord_plot()`. Faceted ordinations show each panel's percent variance 
  or stress in its strip.
* Fixed: `bdiv_boxplot(x = "==Sex")` captioned "within: Sex, Sex", and 
  deprecation warnings named `eval()` rather than the function called.
* Named `colors`, `shapes`, and `patterns` no longer need a value for every 
  group, and names for groups that aren't in the plot are ignored, so one 
  map can serve every subset of a dataset. Groups without a name get the 
  default palette, skipping values already in use. A warning is given only 
  when an unrecognised name meets a group left without a value. Between-group
  pairs such as "Female vs Male" can be named in either order. Heatmap track 
  colors work the same way, and the documented shorthand 
  `tracks = list('Body Site' = c(Saliva = "red", ...))` works again.
* `taxa_forest(colors = )` accepts a named vector, keyed by `"Increased"`, 
  `"Decreased"`, and `"Not significant"`, or for two groups by the `stat.by` 
  groups. A named vector of two or three colors was previously applied by 
  position.

## Figure style

* New `theme_rbiom()`: one publication-ready theme for every rbiom plot, 
  following common journal figure guidelines (Nature, Science, Cell). Black
  sans-serif text sized relative to `base_size` (11 pt; previously 14 pt), 
  open x/y axes with outward ticks, no gridlines or background shading, 
  facet labels without boxes, and titles and captions aligned with the 
  figure's left edge. A plot's `$code` now ends with `theme_rbiom()` instead
  of `theme_bw()` plus a block of `theme()` settings. Options: `frame` 
  (`"axes"`, `"box"`, `"none"`), `grid`, `base_size`, `base_family`; pass 
  them to any plot function, e.g. `adiv_boxplot(biom, base_size = 8)`. Use 
  `theme_rbiom()` on your own ggplots to match rbiom's figures.
* Default colors are now Okabe-Ito for up to 8 groups. The previous default
  (`"classic"`, the D3/Tableau 10 palette) has red/green and orange/green 
  pairs that are hard to tell apart with color vision deficiencies. The 
  `"okabe"` palette is reordered so that its first colors contrast most, 
  with yellow (hard to see on white) last. The "Other" group in 
  `taxa_stacked()` is light gray.
* Default shapes are solid (circle, triangle, square, diamond) before 
  outlines.
* Every facet panel has its own x and y axis lines and ticks (e.g. with 
  `facet.ncol = 1`, panels above the bottom one previously had no x-axis 
  line). Tick labels stay on the outer panels unless the scale is free. Pass
  `facet.axes = "margins"` for ggplot2's default. Boxplots without `x` no 
  longer have a lone, unlabeled x-axis tick.
* p-values in figures have two significant digits, are never in scientific
  notation (`p = 2e-04`), and are shown as `p < 0.001` below 0.001.
* Ordination facet labels no longer overflow: statistics are shortened to 
  `F = 10.7, R² = 0.63, p = 0.001` (the caption names the test). Taxa labels 
  have a fixed size; mean points are scaled by abundance.
* Flipped boxplots draw p-values alongside their brackets instead of across 
  them. Heatmap and `rare_stacked()` sample labels are rotated 90° when 
  there are many of them.

## Breaking changes

* `write_*()` functions refuse to overwrite an existing file unless 
  `overwrite = TRUE`. Text output uses LF line endings on all platforms.
* The plyr functions `.`, `ddply`, `ldply`, and `llply` are no longer 
  re-exported.
* Removed `transform = "percent"` (never worked), and the never-functional 
  `apcoa` and tree-plotting code.
* Removed the `rbiom.cache_hash` option; it is ignored if set. Cache keys 
  are always `rlang::hash()`, as `biom$hash` already was. A custom hash 
  function could only make keys collide or, if it didn't return 32 hex 
  digits, stop the cache from ever being trimmed.
* `bdiv = "unifrac"` is now an error because it is ambiguous; name the 
  variant (e.g. `"unweighted_unifrac"`).
* `ord = "tSNE"` uses `Rtsne::Rtsne()` instead of `tsne::tsne()`. It is 
  about 8-15x faster from a few hundred samples up and keeps more of each 
  sample's nearest neighbors together. Coordinates differ from earlier 
  versions for the same `seed`. Rtsne arguments such as `theta`, 
  `max_iter`, `Y_init`, and `num_threads` can be passed through `...`; the 
  old tsne arguments (`initial_config`, `whiten`, `epoch`, ...) are 
  ignored. `perplexity` now defaults to `(n - 1) / 3` for fewer than 91 
  samples, since the previous default of 30 cannot be reached with so few.
* Statistics table columns changed as described above.

## Dependencies

* Removed: ape, crayon, ggtext, labeling, pillar, plyr, readr. Moved to 
  Suggests: fillpattern, ggbeeswarm, ggrepel, readxl. Packages that were 
  always required but undeclared (cli, glue, lifecycle, Matrix, R6, rlang, 
  scales, tibble) are now in Imports, along with cluster and permute (both 
  already installed with vegan).
* Suggests Rtsne instead of tsne.
* Suggests lmerTest and lme4 (for `taxa_da()` with `strata`).
* Requires ecodive (>= 2.4.0).


# rbiom 3.1.1

* `as_rbiom()` and `read_biom()` now accept count matrices that are triangular, diagonal, or 1x1.
* `taxa_matrix()` and `taxa_table()` no longer fail with "No taxa match the criteria" when `taxa` is a number and a sample has no abundance at the requested rank, e.g. with `unc = "drop"`.


# rbiom 3.1.0

* Can now read and write `biom` objects from the `biomformat` package (microbiome/mia#818).
* Fixed handling of objects from the `TreeSummarizedExperiment` package (microbiome/mia#822).


# rbiom 3.0.1

* Re-exports `ecodive::list_metrics`, `ecodive::match_metric`, and `ggplot2::ggsave`.


# rbiom 3.0.0

## Announcements

The code for `rbiom`'s very fast implementations of diversity metrics,
rarefaction, and tree importing has been spun off into its own dependency-free
package, [`ecodive`](https://github.com/cmmr/ecodive), to encourage adoption by
other R packages. This also means `rbiom` will be easier to install going
forward as it will no longer need to be compiled.


## Breaking Changes

* `weighted` and `normalized` arguments are no longer used by any functions.
* Diversity metric names will now be all lowercase in returned matrices, tables, etc.

The following functions have been deprecated since `rbiom` v2.0 and are now
removed completely. They are generally functions which have been renamed for
improved organization.

`alpha.div`, `as.percent`, `beta.div`, `comments`, `comments<-`, `counts`,
`counts<-`, `depth`, `depths_barplot`, `has.phylogeny`, `has.sequences`, `id`,
`id<-`, `info`, `is.rarefied`, `metadata`, `metadata<-`, `nsamples`, `ntaxa`,
`phylogeny`, `phylogeny<-`, `rarefy`, `read.biom`, `read_biom`, `read.fasta`,
`read.tree`, `repair`, `sample.names`, `sample.names<-`, `sample_subset`,
`sample.sums`, `sequences`, `sequences<-`, `subtree`, `taxa.means`,
`taxa.names`, `taxa.names<-`, `taxa.ranks`, `taxa.ranks<-`, `taxa.rollup`,
`taxa.sums`, `taxonomy`, `taxonomy<-`, `tips`, `top.taxa`, `top_taxa`,
`unifrac`, `write.biom`, `write.fasta`, `write.tree`, `write.xlsx`


## Additions

* `adiv_vector()` returns a simple named vector of alpha diversity values.
* Replaced `rhdf5` dependency with `h5lite`.
* Updated list of alpha diversity metrics: `c("ace", "berger", "brillouin", "chao1", "faith", "fisher", "simpson", "inv_simpson", "margalef", "mcintosh", "menhinick", "observed", "shannon", "squares")`
* Updated list of beta diversity metrics: `c("aitchison", "bhattacharyya", "bray", "canberra", "chebyshev", "chord", "clark", "sorensen", "divergence", "euclidean", "generalized_unifrac", "gower", "hamming", "hellinger", "horn", "jaccard", "jensen", "jsd", "lorentzian", "manhattan", "matusita", "minkowski", "morisita", "motyka", "normalized_unifrac", "ochiai", "psym_chisq", "soergel", "squared_chisq", "squared_chord", "squared_euclidean", "topsoe", "unweighted_unifrac", "variance_adjusted_unifrac", "wave_hedges", "weighted_unifrac")`



# rbiom 2.2.1

* Compatibility fixes for ggplot2 4.0.0 #30


# rbiom 2.2.0

* Satisfy new R CMD check: no setting attributes on built-in functions. #28


# rbiom 2.1.2

* Fix for `write_biom(format='json')` where `{...}` should be `[...]`.
* Restore `rline` functionality in `rare_corrplot()`.


# rbiom 2.1.1

* `unifrac()` no longer normalizes weighted values, for back-compatibility.
* Option `underscores` allows keeping underscores in tree IDs.


# rbiom 2.0.13

* Major release with significant new features.
* Plotting added.
* Statistics added.
* Caching added.
* Clearer naming scheme for functions. Previous names still work but are deprecated.
* Generation of human-readable ggplot2 code for modifying plots outside of rbiom.
* Provenance tracking for BIOM objects and derivatives.


# rbiom 1.0.3

* Corrects for breaking changes in 'rhdf5' package.


# rbiom 1.0.2

* Improved compatibility with Debian and Solaris.
* The 'rhdf5' package is now an optional dependency.
* The select() and subset() functions now subset sequences too.


# rbiom 1.0.0

* Initial Release
