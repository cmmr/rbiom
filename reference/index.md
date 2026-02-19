# Package index

## Input / Output

Read and write BIOM, FASTA, Newick, and tabular data. Includes tools to
convert to/from other popular formats like phyloseq,
TreeSummarizedExperiment, and QIIME 2.

- [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md) :
  Convert a variety of data types to an rbiom object.
- [`read_biom()`](https://cmmr.github.io/rbiom/reference/read_biom.md) :
  Parse counts, metadata, taxonomy, and phylogeny from a BIOM file.
- [`read_fasta()`](https://cmmr.github.io/rbiom/reference/read_fasta.md)
  : Parse a fasta file into a named character vector.
- [`read_tree()`](https://cmmr.github.io/rbiom/reference/read_tree.md) :
  Read a newick formatted phylogenetic tree.
- [`write_mothur()`](https://cmmr.github.io/rbiom/reference/export.md)
  [`write_qiime2()`](https://cmmr.github.io/rbiom/reference/export.md) :
  Export data to QIIME 2 or mothur.
- [`write_biom()`](https://cmmr.github.io/rbiom/reference/write_biom.md)
  [`write_metadata()`](https://cmmr.github.io/rbiom/reference/write_biom.md)
  [`write_counts()`](https://cmmr.github.io/rbiom/reference/write_biom.md)
  [`write_taxonomy()`](https://cmmr.github.io/rbiom/reference/write_biom.md)
  [`write_fasta()`](https://cmmr.github.io/rbiom/reference/write_biom.md)
  [`write_tree()`](https://cmmr.github.io/rbiom/reference/write_biom.md)
  [`write_xlsx()`](https://cmmr.github.io/rbiom/reference/write_biom.md)
  : Save an rbiom object to a file.
- [`convert_to_SE()`](https://cmmr.github.io/rbiom/reference/convert_to.md)
  [`convert_to_TSE()`](https://cmmr.github.io/rbiom/reference/convert_to.md)
  [`convert_to_phyloseq()`](https://cmmr.github.io/rbiom/reference/convert_to.md)
  [`convert_to_animalcules()`](https://cmmr.github.io/rbiom/reference/convert_to.md)
  : Convert biom data to an external package class.

## The rbiom Object

The [rbiom
object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md) serves
as the core container for data, providing user-friendly accessors for
`$counts`, `$metadata`, `$taxonomy`, `$tree`, and `$sequences`.

- [`rbiom_objects`](https://cmmr.github.io/rbiom/reference/rbiom_objects.md)
  : Working with rbiom Objects.
- [`as.list(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/as.list.rbiom.md)
  : Convert an rbiom object to a base R list.
- [`as.matrix(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/as.matrix.rbiom.md)
  : Convert an rbiom object to a simple count matrix.

## Sample Metadata

Access and manipulate sample metadata directly from the rbiom object
using familiar dplyr-like verbs.

- [`pull(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/pull.rbiom.md)
  : Map sample names to metadata field values.
- [`with(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/with.md)
  [`within(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/with.md)
  : Evaluate expressions on metadata.
- [`mutate(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/modify_metadata.md)
  [`rename(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/modify_metadata.md)
  : Create, modify, and delete metadata fields.
- [`glimpse(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/glimpse.rbiom.md)
  : Get a glimpse of your metadata.

## Subsetting

Filter samples based on metadata conditions or apply functions to
subsets of the data (split-apply-combine).

- [`subset(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/subset.md)
  [`` `[`( ``*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/subset.md)
  [`na.omit(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/subset.md)
  [`subset_taxa()`](https://cmmr.github.io/rbiom/reference/subset.md) :
  Subset an rbiom object by sample names, OTU names, metadata, or
  taxonomy.
- [`slice(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/slice_metadata.md)
  [`slice_head(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/slice_metadata.md)
  [`slice_tail(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/slice_metadata.md)
  [`slice_min(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/slice_metadata.md)
  [`slice_max(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/slice_metadata.md)
  [`slice_sample(`*`<rbiom>`*`)`](https://cmmr.github.io/rbiom/reference/slice_metadata.md)
  : Subset to a specific number of samples.
- [`bdply()`](https://cmmr.github.io/rbiom/reference/bdply.md)
  [`blply()`](https://cmmr.github.io/rbiom/reference/bdply.md) : Apply a
  function to each subset of an rbiom object.

## Taxa Abundance

Aggregates counts to specific taxonomic ranks (e.g. Genus, Phylum) to
test for differential abundance and generate visualizations like stacked
bar charts, heatmaps, and boxplots.

- [`taxa_boxplot()`](https://cmmr.github.io/rbiom/reference/taxa_boxplot.md)
  : Visualize BIOM data with boxplots.
- [`taxa_clusters()`](https://cmmr.github.io/rbiom/reference/taxa_clusters.md)
  : Cluster samples by taxa abundances k-means.
- [`taxa_corrplot()`](https://cmmr.github.io/rbiom/reference/taxa_corrplot.md)
  : Visualize taxa abundance with scatterplots and trendlines.
- [`taxa_heatmap()`](https://cmmr.github.io/rbiom/reference/taxa_heatmap.md)
  : Display taxa abundances as a heatmap.
- [`taxa_map()`](https://cmmr.github.io/rbiom/reference/taxa_map.md) :
  Map OTUs names to taxa names at a given rank.
- [`taxa_table()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)
  [`taxa_matrix()`](https://cmmr.github.io/rbiom/reference/taxa_matrix.md)
  : Taxa abundances per sample.
- [`taxa_stacked()`](https://cmmr.github.io/rbiom/reference/taxa_stacked.md)
  : Display taxa abundances as a stacked bar graph.
- [`taxa_stats()`](https://cmmr.github.io/rbiom/reference/taxa_stats.md)
  : Test taxa abundances for associations with metadata.
- [`taxa_sums()`](https://cmmr.github.io/rbiom/reference/taxa_sums.md)
  [`taxa_means()`](https://cmmr.github.io/rbiom/reference/taxa_sums.md)
  [`taxa_apply()`](https://cmmr.github.io/rbiom/reference/taxa_sums.md)
  : Get summary taxa abundances.

## Alpha Diversity

Calculate within-sample diversity metrics (Shannon, Simpson, Richness,
etc.) and correlate them with metadata using boxplots and statistical
tests.

- [`adiv_boxplot()`](https://cmmr.github.io/rbiom/reference/adiv_boxplot.md)
  : Visualize alpha diversity with boxplots.
- [`adiv_corrplot()`](https://cmmr.github.io/rbiom/reference/adiv_corrplot.md)
  : Visualize alpha diversity with scatterplots and trendlines.
- [`adiv_stats()`](https://cmmr.github.io/rbiom/reference/adiv_stats.md)
  : Test alpha diversity for associations with metadata.
- [`adiv_table()`](https://cmmr.github.io/rbiom/reference/adiv_table.md)
  [`adiv_matrix()`](https://cmmr.github.io/rbiom/reference/adiv_table.md)
  [`adiv_vector()`](https://cmmr.github.io/rbiom/reference/adiv_table.md)
  : Calculate the alpha diversity of each sample.
- [`sample_sums()`](https://cmmr.github.io/rbiom/reference/sample_sums.md)
  [`sample_apply()`](https://cmmr.github.io/rbiom/reference/sample_sums.md)
  : Summarize the taxa observations in each sample.

## Beta Diversity

Compare samples to one another using distance metrics (Bray-Curtis,
UniFrac, etc.) and visualize patterns with Ordinations (PCoA, NMDS) and
clustering.

- [`bdiv_boxplot()`](https://cmmr.github.io/rbiom/reference/bdiv_boxplot.md)
  : Visualize BIOM data with boxplots.
- [`bdiv_clusters()`](https://cmmr.github.io/rbiom/reference/bdiv_clusters.md)
  : Cluster samples by beta diversity k-means.
- [`bdiv_corrplot()`](https://cmmr.github.io/rbiom/reference/bdiv_corrplot.md)
  : Visualize beta diversity with scatterplots and trendlines.
- [`bdiv_heatmap()`](https://cmmr.github.io/rbiom/reference/bdiv_heatmap.md)
  : Display beta diversities in an all vs all grid.
- [`bdiv_ord_plot()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_plot.md)
  : Ordinate samples and taxa on a 2D plane based on beta diversity
  distances.
- [`bdiv_ord_table()`](https://cmmr.github.io/rbiom/reference/bdiv_ord_table.md)
  : Calculate PCoA and other ordinations, including taxa biplots and
  statistics.
- [`bdiv_stats()`](https://cmmr.github.io/rbiom/reference/bdiv_stats.md)
  : Test beta diversity for associations with metadata.
- [`bdiv_table()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md)
  [`bdiv_matrix()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md)
  [`bdiv_distmat()`](https://cmmr.github.io/rbiom/reference/bdiv_table.md)
  : Distance / dissimilarity between samples.

## Rarefaction Plots

Visualize how sequencing depth affects diversity metrics and taxa
discovery (rarefaction curves).

- [`rare_corrplot()`](https://cmmr.github.io/rbiom/reference/rare_corrplot.md)
  : Visualize rarefaction curves with scatterplots and trendlines.
- [`rare_multiplot()`](https://cmmr.github.io/rbiom/reference/rare_multiplot.md)
  : Combines rare_corrplot and rare_stacked into a single figure.
- [`rare_stacked()`](https://cmmr.github.io/rbiom/reference/rare_stacked.md)
  : Visualize the number of observations per sample.

## Transforming Counts

Normalize and transform count data, including rarefaction to a fixed
depth, converting to relative abundance, and inflation/rescaling.

- [`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md) :
  Rarefy Counts to a Constant Depth
- [`biom_inflate()`](https://cmmr.github.io/rbiom/reference/biom_inflate.md)
  : Inflate Relative Abundances to Integer Counts
- [`biom_relativize()`](https://cmmr.github.io/rbiom/reference/biom_relativize.md)
  : Relativize Counts to Proportions
- [`biom_rescale()`](https://cmmr.github.io/rbiom/reference/biom_rescale.md)
  : Rescale Counts to a Specific Range
- [`suggest_rarefy_depth()`](https://cmmr.github.io/rbiom/reference/suggest_rarefy_depth.md)
  : Suggest Rarefaction Depth
- [`suggest_inflate_depths()`](https://cmmr.github.io/rbiom/reference/suggest_inflate_depths.md)
  : Suggest Inflation Depths
- [`biom_merge()`](https://cmmr.github.io/rbiom/reference/biom_merge.md)
  : Combine several rbiom objects into one.

## Low Level Functions

Generic statistical and plotting functions that work on standard R data
structures (matrices, data.frames, trees) rather than rbiom objects.

- [`distmat_ord_table()`](https://cmmr.github.io/rbiom/reference/distmat_ord_table.md)
  : Run ordinations on a distance matrix.
- [`distmat_stats()`](https://cmmr.github.io/rbiom/reference/distmat_stats.md)
  : Run statistics on a distance matrix vs a categorical or numeric
  variable.
- [`stats_boxplot()`](https://cmmr.github.io/rbiom/reference/stats_boxplot.md)
  : Visualize categorical metadata effects on numeric values.
- [`stats_corrplot()`](https://cmmr.github.io/rbiom/reference/stats_corrplot.md)
  : Visualize regression with scatterplots and trendlines.
- [`stats_table()`](https://cmmr.github.io/rbiom/reference/stats_table.md)
  : Run non-parametric statistics on a data.frame.
- [`tree_subset()`](https://cmmr.github.io/rbiom/reference/tree_subset.md)
  : Create a subtree by specifying tips to keep.
- [`plot_heatmap()`](https://cmmr.github.io/rbiom/reference/plot_heatmap.md)
  : Create a heatmap with tracks and dendrograms from any matrix.

## Datasets

Built-in microbiome datasets for testing and examples.

- [`hmp50`](https://cmmr.github.io/rbiom/reference/hmp50.md) : Human
  Microbiome Project - demo dataset (n = 50)
- [`gems`](https://cmmr.github.io/rbiom/reference/gems.md) : Global
  Enteric Multicenter Study (n = 1,006)
- [`babies`](https://cmmr.github.io/rbiom/reference/babies.md) :
  Longitudinal Stool Samples from Infants (n = 2,684)

## Deprecated Functions

Functions that have been superseded by newer alternatives.

- [`mtx_rarefy()`](https://cmmr.github.io/rbiom/reference/matrix_ops.md)
  [`mtx_percent()`](https://cmmr.github.io/rbiom/reference/matrix_ops.md)
  [`mtx_rescale()`](https://cmmr.github.io/rbiom/reference/matrix_ops.md)
  [`rarefy_cols()`](https://cmmr.github.io/rbiom/reference/matrix_ops.md)
  [`rescale_rows()`](https://cmmr.github.io/rbiom/reference/matrix_ops.md)
  [`rescale_cols()`](https://cmmr.github.io/rbiom/reference/matrix_ops.md)
  : Deprecated matrix transformations
