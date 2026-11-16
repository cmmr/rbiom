#' Taxa shared between groups, and the core microbiome.
#'
#' Lists which taxa are found in which groups of samples. With a prevalence
#' cut-off (`prev`) the taxa listed for a group are its core microbiome.
#' [taxa_upset()] plots the result.
#'
#' @inherit documentation_default
#'
#' @family taxa_abundance
#'
#' @param stat.by   Dataset field with the groups to compare. Must be
#'        categorical. Samples with `NA` in this field are left out. `NULL`
#'        puts every sample in a single group named `"All"`, which gives a
#'        plain core microbiome list. Default: `NULL`
#'
#' @param rank   Which taxa rank to use. E.g. `"Phylum"`, `"Genus"`,
#'        `".otu"`, etc. An integer can also be given, where `1` is the
#'        highest rank, `2` is the second highest, `-1` is the lowest rank,
#'        `-2` is the second lowest, and `0` is the OTU "rank". Run
#'        `biom$ranks` to see all options for a given rbiom object.
#'        Default: `-1`.
#'
#' @param detect   Detection threshold. A taxon is present in a sample when
#'        its relative abundance there is greater than `detect`, a fraction
#'        from 0 to 1 (`0.001` is 0.1%). The default counts any observation.
#'        Default: `0`
#'
#' @param prev   Prevalence threshold. A taxon belongs to a group when it is
#'        present in at least this fraction of the group's samples, and in
#'        at least one of them. So `0` needs a single sample, `0.9` needs
#'        90% of the group's samples, and `1` needs all of them.
#'        Default: `0`
#'
#' @details
#' Whether a rare taxon is observed in a sample depends on how deeply the
#' sample was sequenced, so give this function rarefied counts (see
#' [rarefy()]).
#'
#' @return A tibble data frame with one row per taxon that belongs to at
#' least one group, and columns:
#' \describe{
#'   \item{`.taxa`}{ The taxon's name. }
#'   \item{One column per group}{ Named after the levels of `stat.by`, in
#'     level order: `TRUE` when the taxon belongs to that group. }
#'   \item{`.n_groups`}{ How many groups the taxon belongs to. }
#'   \item{`.set`}{ The combination of groups the taxon belongs to, e.g.
#'     `"Saliva & Stool"`. }
#'   \item{`.prev.<group>`}{ One column per group: the fraction of the
#'     group's samples in which the taxon is present. }
#' }
#' Rows are sorted by `.set` (the largest set first), then by name.
#'
#' The number of taxa in each set is attached as `attr(, "sets")`: a tibble
#' with `.set`, `.n` (taxa belonging to exactly that combination of
#' groups), `.n_groups`, and one logical column per group, sorted from the
#' largest set to the smallest.
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
#'     # Genera found at each body site.
#'     shared <- taxa_shared(biom, stat.by = "Body Site")
#'     shared
#'
#'     attr(shared, "sets")
#'
#'     # Core phyla: above 0.1% in at least 90% of a body site's samples.
#'     taxa_shared(biom, "Body Site", rank = "Phylum", detect = 0.001, prev = 0.9)
#'
#'     # Core genera across all samples.
#'     taxa_shared(biom, prev = 0.5)
#'

taxa_shared <- function (
    biom, stat.by = NULL, rank = -1, detect = 0, prev = 0,
    unc = "singly", lineage = FALSE ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  cmd    <- current_cmd('taxa_shared', params)
  schema_check('taxa_shared', params)


  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  validate_rank()
  validate_biom_field('stat.by', null_ok = TRUE, col_type = "cat")
  validate_var_choices('unc', CHOICES$unc)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('taxa_shared', params)
  tbl        <- get_cache_value(cache_file)

  if (is.null(tbl)) {

    #________________________________________________________
    # Assign samples to groups.
    #________________________________________________________
    if (is.null(stat.by)) {
      groups <- stats::setNames(factor(rep("All", biom$n_samples)), biom$samples)

    } else {
      groups <- pull(biom, stat.by)
      groups <- groups[!is.na(groups)]
      groups <- droplevels(as.factor(groups))

      if (length(groups) == 0)
        cli_abort("All samples have NA values in {.field {stat.by}}.")
    }

    lvls     <- levels(groups)
    reserved <- lvls[lvls %in% c('.taxa', '.n_groups', '.set', '.n') | startsWith(lvls, '.prev.')]

    if (length(reserved) > 0)
      cli_abort("{.field {stat.by}} has {qty(length(reserved))} level{?s} named {.val {reserved}}, which {?is a/are} column name{?s} of the result. Rename {?it/them} first.")


    #________________________________________________________
    # Presence in each sample, then prevalence in each group.
    #________________________________________________________
    mtx <- taxa_matrix(
      biom      = biom,
      rank      = rank,
      lineage   = lineage,
      sparse    = TRUE,
      unc       = unc,
      transform = "percent" )

    present <- mtx > detect
    taxa    <- rownames(mtx)

    prevalence <- matrix(0,     nrow = length(taxa), ncol = length(lvls), dimnames = list(taxa, lvls))
    member     <- matrix(FALSE, nrow = length(taxa), ncol = length(lvls), dimnames = list(taxa, lvls))

    for (lvl in lvls) {
      ids <- names(groups)[groups == lvl]
      n   <- rowSums(present[, ids, drop = FALSE])
      prevalence[, lvl] <- n / length(ids)
      member[, lvl]     <- n > 0 & n >= prev * length(ids) - 1e-9
    }

    keep       <- rowSums(member) > 0
    member     <- member[keep, , drop = FALSE]
    prevalence <- prevalence[keep, , drop = FALSE]
    taxa       <- taxa[keep]


    #________________________________________________________
    # The set (combination of groups) each taxon falls in.
    #________________________________________________________
    key   <- apply(member, 1L, function (x) paste(as.integer(x), collapse = ""))
    label <- apply(member, 1L, function (x) paste(lvls[x], collapse = " & "))
    key   <- as.character(key)   # apply() gives list() for zero rows
    label <- as.character(label)

    first <- which(!duplicated(key))
    sets  <- tibble(
      '.set'      = label[first],
      '.n'        = as.integer(table(key)[key[first]]),
      '.n_groups' = as.integer(rowSums(member[first, , drop = FALSE])) )

    for (lvl in lvls)
      sets[[lvl]] <- unname(member[first, lvl])

    # Largest first; then sets of more groups; then by group order.
    sets <- sets[order(-sets[['.n']], -sets[['.n_groups']], -xtfrm(key[first])), ]
    sets[['.set']] %<>% {factor(., levels = .)}


    #________________________________________________________
    # One row per taxon.
    #________________________________________________________
    tbl <- tibble('.taxa' = taxa)

    for (lvl in lvls)
      tbl[[lvl]] <- unname(member[, lvl])

    tbl[['.n_groups']] <- as.integer(rowSums(member))
    tbl[['.set']]      <- factor(label, levels = levels(sets[['.set']]))

    for (lvl in lvls)
      tbl[[paste0('.prev.', lvl)]] <- unname(prevalence[, lvl])

    tbl <- tbl[order(tbl[['.set']], tolower(taxa), taxa), ]
    tbl[['.taxa']] %<>% {factor(., levels = .)}

    tbl %<>% as_rbiom_tbl()
    attr(tbl, 'sets') <- sets

    set_cache_value(cache_file, tbl)
  }


  return (prov_attach(tbl, cmd, biom))
}
