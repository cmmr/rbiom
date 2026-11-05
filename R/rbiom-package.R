#' rbiom: Integrated Analysis and Visualization of Microbiome Data
#'
#' A toolkit for working with Biological Observation Matrix (BIOM) files.
#' Features include reading/writing all BIOM formats, rarefaction, alpha
#' diversity, beta diversity (including UniFrac), summarizing counts by
#' taxonomic level, sample subsetting, statistics, and visualization.
#' Standalone functions for reading, writing, and subsetting phylogenetic
#' trees are also provided.
#'
#' @section Reproducibility:
#'
#' Every rbiom object records the steps that produced it, and every table,
#' statistic, and plot carries the R code needed to regenerate it. See
#' `biom$code`, `table$code`, `stats$code`, and `plot$code`. The
#' [rbiom objects][rbiom_objects] page has details.
#'
#' @section Multithreading:
#'
#' Diversity calculations are performed by the 'ecodive' package, which runs
#' on all available cores by default. Use the `cpus` argument of individual
#' functions to limit the number of threads, e.g. `adiv_table(biom, cpus = 1)`.
#'
#' @keywords internal
#' @aliases rbiom-package
"_PACKAGE"

NULL
