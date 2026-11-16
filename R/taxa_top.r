#' The most abundant taxa overall, per group, or per sample.
#'
#' [taxa_table()] with `taxa = 5` keeps the five taxa that are most abundant
#' across all samples. `taxa_top()` ranks the taxa separately within each
#' group of samples (or each sample), so the taxa listed can differ from one
#' group to the next.
#'
#' @inherit documentation_default
#'
#' @family taxa_abundance
#'
#' @param n   How many taxa to list for each group. Default: `5`
#'
#' @param by   Dataset field(s) that define the groups of samples to rank
#'        taxa within. Must be categorical. Several fields give one group per
#'        combination of their values. Use `".sample"` to rank taxa within
#'        each sample. `NULL` ranks taxa across all samples together.
#'        Samples with `NA` in a `by` field are left out. Default: `NULL`
#'
#' @param transform   Units for the abundances in the result. Options are:
#'        `c("none", "rank", "log", "log1p", "sqrt", "percent")`. `"percent"`
#'        gives relative abundances as fractions (`0.25` is 25%), and
#'        `"none"` gives counts. See [taxa_table()] for the others. The
#'        choice does not change which taxa are listed. Default: `"percent"`
#'
#' @details
#' Taxa are ranked by their mean relative abundance in the group's samples,
#' as [taxa_table()] ranks them for its `taxa` argument. Ties are broken
#' alphabetically by taxon name, so the result never depends on a random
#' seed. Taxa that are absent from a group are not listed, so a group can
#' have fewer than `n` rows.
#'
#' @return A tibble data frame with one row per listed taxon and columns:
#' \describe{
#'   \item{`.rank`}{ The taxonomic rank. }
#'   \item{`by` field(s)}{ The group. Not present when `by = NULL`. }
#'   \item{`.top`}{ The taxon's position in its group: `1` for the most
#'     abundant, up to `n`. }
#'   \item{`.taxa`}{ The taxon's name. }
#'   \item{`.mean`}{ Mean abundance over the group's samples, in `transform`
#'     units. }
#'   \item{`.prevalence`}{ Fraction of the group's samples in which the taxon
#'     was observed. }
#' }
#' When `by` includes `".sample"`, each group is one sample, and the table
#' has `.abundance` (in `transform` units) in place of `.mean` and
#' `.prevalence`.
#'
#' The result carries `$code`: the R code that regenerates it, including
#' the steps that produced `biom`.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(hmp50)
#'
#'     # The five most abundant genera overall.
#'     taxa_top(biom)
#'
#'     # The three most abundant phyla at each body site.
#'     taxa_top(biom, rank = "Phylum", n = 3, by = "Body Site")
#'
#'     # The dominant genus in each sample.
#'     taxa_top(biom, n = 1, by = ".sample")
#'

taxa_top <- function (
    biom, rank = -1, n = 5, by = NULL, unc = "singly",
    lineage = FALSE, transform = "percent" ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  cmd    <- current_cmd('taxa_top', params)
  schema_check('taxa_top', params)


  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  validate_rank(max = Inf)
  validate_biom_field('by', max = Inf, null_ok = TRUE, col_type = "cat")
  validate_var_choices('transform', CHOICES$transform)
  validate_var_choices('unc', CHOICES$unc)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('taxa_top', params)
  tbl        <- get_cache_value(cache_file)

  if (is.null(tbl)) {

    by <- unique(by)
    md <- as.data.frame(biom$metadata, stringsAsFactors = FALSE)
    rownames(md) <- md[['.sample']]
    md <- md[biom$samples, , drop = FALSE]

    per_sample <- '.sample' %in% by


    #________________________________________________________
    # The samples in each group, and each group's `by` values.
    #________________________________________________________
    if (is.null(by)) {

      groups <- list(biom$samples)
      keys   <- NULL

    } else {

      md <- md[stats::complete.cases(md[, by, drop = FALSE]), , drop = FALSE]

      if (nrow(md) == 0)
        cli_abort("All samples have NA values in {.field {by}}.")

      if (per_sample) {
        md[['.sample']] <- factor(md[['.sample']], levels = md[['.sample']])
        groups <- as.list(rownames(md))

      } else {
        grp    <- interaction(lapply(md[, by, drop = FALSE], as.factor), drop = TRUE, lex.order = TRUE)
        groups <- unname(split(rownames(md), grp))
      }

      keys <- md[vapply(groups, `[[`, character(1), 1), by, drop = FALSE]
      for (i in by) if (is.factor(keys[[i]])) keys[[i]] <- droplevels(keys[[i]])
    }


    #________________________________________________________
    # Return multiple ranks in a single table.
    #________________________________________________________
    tbl <- NULL

    for (r in rank) {

      raw <- taxa_matrix(biom = biom, rank = r, lineage = lineage, sparse = TRUE, unc = unc)

      val <- if (eq(transform, 'none')) raw else taxa_matrix(
        biom      = biom,
        rank      = r,
        lineage   = lineage,
        sparse    = TRUE,
        unc       = unc,
        transform = transform,
        ties      = "average" )

      val <- val[rownames(raw), colnames(raw), drop = FALSE]

      # Relative abundances, as top_taxa() computes them. Samples with no
      # abundance at this rank (e.g. from unc = "drop") can't be ranked.
      totals <- colSums(raw)
      taxa   <- rownames(raw)

      for (g in seq_along(groups)) {

        ids <- groups[[g]]
        ids <- ids[totals[ids] > 0]
        if (length(ids) == 0) next

        rel <- as.matrix(raw[, ids, drop = FALSE])
        rel <- rowMeans(t(t(rel) / totals[ids]))

        # Most abundant first; ties alphabetically.
        top <- order(-rel, tolower(taxa), taxa)
        top <- head(top[rel[top] > 0], n)
        if (length(top) == 0) next

        row <- tibble('.rank' = rep(r, length(top)))

        for (i in by)
          row[[i]] <- rep(keys[[i]][[g]], length(top))

        row[['.top']]  <- seq_along(top)
        row[['.taxa']] <- taxa[top]

        if (per_sample) {
          row[['.abundance']] <- as.numeric(val[top, ids])

        } else {
          ids <- groups[[g]]
          row[['.mean']]       <- unname(rowMeans(val[top, ids, drop = FALSE]))
          row[['.prevalence']] <- unname(rowSums(raw[top, ids, drop = FALSE] > 0)) / length(ids)
        }

        tbl %<>% dplyr::bind_rows(row)
      }
    }

    if (is.null(tbl))
      cli_abort("No taxa with any abundance at rank {.val {rank}}.")

    tbl[['.rank']] %<>%  factor(., levels = rank)
    tbl[['.taxa']] %<>% {factor(., levels = unique(.))}

    tbl %<>% as_rbiom_tbl()
    attr(tbl, 'response') <- if (per_sample) ".abundance" else ".mean"

    set_cache_value(cache_file, tbl)
  }


  return (prov_attach(tbl, cmd, biom))
}
