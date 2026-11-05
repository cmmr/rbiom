
#' Combine several rbiom objects into one.
#'
#' WARNING: It is generally ill-advised to merge BIOM datasets, as OTUs
#' mappings are dependent on upstream clustering and are not equivalent
#' between BIOM files.
#'
#' @inherit documentation_return.biom return
#'
#' @family biom
#'
#' @param ...  Any number of rbiom objects, or anything accepted by the
#'        `biom` parameter of [as_rbiom()] (matrices, file names, lists with a
#'        `counts` element, phyloseq objects, ...). Lists of such inputs are
#'        flattened.
#'
#' @param metadata,taxonomy,tree,sequences,id,comment   Replace the
#'        corresponding data in the merged rbiom object with these values.
#'        The default, `NA`, derives the component from the inputs (see
#'        details). Set to `NULL` to omit the component.
#'
#' @details
#'
#' Samples present in more than one input must have identical counts;
#' otherwise an error is raised. OTUs are matched by name.
#'
#' With the default `NA` arguments:
#'
#' * `metadata` and `taxonomy` are the union of the inputs' tables. A field
#'   whose type differs between inputs (for example numeric in one and
#'   categorical in another) is converted to categorical, with a single
#'   warning listing the affected fields. An OTU with different lineages in
#'   different inputs keeps the most detailed one (one summarizing warning).
#' * `sequences` are the union of the inputs' sequences. If only some inputs
#'   have sequences, OTUs without one get `NA` and a warning is issued once.
#' * `tree` is taken from the first input whose tree contains every OTU of
#'   the merged object (pruned to those OTUs). This is the case when the
#'   inputs are subsets of one dataset. If no input tree covers all of the
#'   OTUs, the merged object has no tree and a message says so; a merged
#'   tree can't be inferred, so supply `tree=` in that case.
#' * `id` and `comment` are kept when all inputs agree, and are otherwise
#'   `"Merged BIOM"` and `""`.
#'
#' The merged object's `$history` nests the histories of its inputs, so its
#' `$code` reproduces the entire pipeline.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     b1 <- as_rbiom(hmp50$counts[,1:4])
#'     b2 <- as_rbiom(hmp50$counts[,5:8])
#'
#'     biom <- biom_merge(b1, b2)
#'     print(biom)
#'
#'     biom$tree     <- hmp50$tree
#'     biom$metadata <- subset(hmp50$metadata, .sample %in% biom$samples)
#'     print(biom)
#'
#'     # Subsets of one dataset keep its tree, and record their history.
#'     biom <- biom_merge(hmp50[1:10], hmp50[11:20])
#'     biom$code

biom_merge <- function (..., metadata = NA, taxonomy = NA, tree = NA, sequences = NA, id = NA, comment = NA) {

  #________________________________________________________
  # Flatten and validate the inputs.
  #________________________________________________________
  inputs <- flatten_merge_inputs(list(...))

  if (length(inputs) == 0)
    cli_abort("No datasets provided to biom_merge().")

  bioms <- lapply(inputs, as_rbiom)
  n     <- length(bioms)


  #________________________________________________________
  # Non-default arguments, for the provenance record.
  #________________________________________________________
  extra <- list()
  for (arg in c('metadata', 'taxonomy', 'tree', 'sequences', 'id', 'comment'))
    if (!is_na_arg(get(arg)))
      extra[arg] <- list(get(arg))


  #________________________________________________________
  # Global OTU / sample names. Samples shared by several inputs
  # must have identical counts; they are taken from the first.
  #________________________________________________________
  all_otus    <- unique(unlist(lapply(bioms, function (b) b$otus),    use.names = FALSE))
  all_samples <- unique(unlist(lapply(bioms, function (b) b$samples), use.names = FALSE))

  owner    <- integer(length(all_samples))
  triplets <- list()

  for (k in seq_len(n)) {

    mtx  <- bioms[[k]]$counts
    sids <- colnames(mtx)
    pos  <- match(sids, all_samples)

    for (j in which(owner[pos] > 0)) {

      sid  <- sids[[j]]
      prev <- bioms[[owner[pos[[j]]]]]$counts

      a <- mtx[, sid];  a <- a[a != 0]
      b <- prev[, sid]; b <- b[b != 0]

      if (!isTRUE(all.equal(a[order(names(a))], b[order(names(b))])))
        cli_abort(c(
          'x' = "Sample {.val {sid}} has different counts in different inputs.",
          'i' = "Samples present in more than one input must be identical." ))
    }

    new <- which(owner[pos] == 0)
    if (length(new) == 0) next
    owner[pos[new]] <- k

    sm <- Matrix::summary(mtx[, new, drop = FALSE])

    triplets[[k]] <- list(
      i = match(rownames(mtx), all_otus)[sm$i],
      j = pos[new][sm$j],
      x = sm$x )
  }

  counts <- sparseMatrix(
    i        = unlist(lapply(triplets, `[[`, "i"), use.names = FALSE),
    j        = unlist(lapply(triplets, `[[`, "j"), use.names = FALSE),
    x        = unlist(lapply(triplets, `[[`, "x"), use.names = FALSE),
    dims     = c(length(all_otus), length(all_samples)),
    dimnames = list(all_otus, all_samples) )


  # Warn if no OTU overlap
  if (n > 1 && length(all_otus) == sum(vapply(bioms, function (b) b$n_otus, 0L)))
    cli_warn("No overlapping OTU names. Likely incompatible datasets.")


  #________________________________________________________
  # Metadata and Taxonomy
  #________________________________________________________
  if (is_na_arg(metadata))
    metadata <- merge_metadata(lapply(bioms, function (b) b$metadata))

  if (is_na_arg(taxonomy))
    taxonomy <- merge_taxonomy(lapply(bioms, function (b) b$taxonomy))


  #________________________________________________________
  # Sequences
  #________________________________________________________
  if (is_na_arg(sequences))
    sequences <- merge_sequences(lapply(bioms, function (b) b$sequences), all_otus)


  #________________________________________________________
  # Tree: the first input tree that covers every merged OTU.
  #________________________________________________________
  if (is_na_arg(tree)) {

    # The `$tree` getter prunes lazily, so an input's stored tree may still
    # cover OTUs beyond its own (e.g. subsets of one dataset). Use the first
    # stored tree that has every OTU of the merged object, pruned to them.
    stored <- lapply(bioms, function (b) b$.__enclos_env__$private$.tree)
    covers <- vapply(stored, function (t) !is.null(t) && all(all_otus %in% t$tip.label), TRUE)

    if (any(covers)) {
      tree <- stored[[which(covers)[[1]]]]
      if (length(tree$tip.label) > length(all_otus) && length(all_otus) >= 2)
        tree <- tree_subset(tree, all_otus)

    } else {
      if (!all(vapply(stored, is.null, TRUE)))
        cli::cli_inform(c('i' = "No input tree covers every OTU of the merged object, so it has no tree. Use `tree=` to set one."))
      tree <- NULL
    }
  }


  #________________________________________________________
  # ID and comment
  #________________________________________________________
  if (is_na_arg(id)) {
    ids   <- unique(vapply(bioms, function (b) b$id, ""))
    named <- setdiff(ids, c('Untitled Dataset', ''))
    id    <- if (length(named) == 1) { named
    } else if (length(ids) == 1)     { ids
    } else                           { "Merged BIOM" }
  }

  if (is_na_arg(comment)) {
    cm      <- setdiff(unique(vapply(bioms, function (b) b$comment, "")), '')
    comment <- if (length(cm) == 1) cm else ""
  }


  #________________________________________________________
  # Assemble. The date is the most recent input date, so that
  # repeated merges of the same inputs are identical.
  #________________________________________________________
  merged <- rbiom$new(
    counts    = counts,
    metadata  = metadata,
    taxonomy  = taxonomy,
    sequences = sequences,
    tree      = tree,
    id        = id,
    comment   = comment,
    date      = max(vapply(bioms, function (b) b$date, "")) )


  #________________________________________________________
  # Provenance: biom <- biom_merge(biom1, biom2, ...), with the
  # inputs' histories nested.
  #________________________________________________________
  pc   <- prov_call('biom_merge', extra)
  call <- as.call(c(
    as.name('biom_merge'),
    lapply(seq_len(n), function (i) as.name(paste0('biom', i))),
    as.list(pc$call)[-1] ))

  prov_set(merged, list())
  prov_log(merged, call = call, reproducible = pc$reproducible, inputs = lapply(bioms, prov_get))


  return (merged)
}



#________________________________________________________
# Is `x` the default `NA` (as opposed to NULL or a value)?
#________________________________________________________
is_na_arg <- function (x) {
  is.atomic(x) && length(x) == 1 && is.na(x)
}



#________________________________________________________
# Turn the `...` of biom_merge() into a flat list of inputs for
# as_rbiom(): nested lists (without a 'counts' element) and character
# vectors of file names are expanded.
#________________________________________________________
flatten_merge_inputs <- function (x) {

  out <- list()

  for (el in x) {

    if (is.list(el) && !is.data.frame(el) && !inherits(el, 'rbiom') && !hasName(el, 'counts')) {
      out <- c(out, flatten_merge_inputs(el))

    } else if (is.character(el) && length(el) > 1) {
      out <- c(out, as.list(el))

    } else {
      out <- c(out, list(el))
    }
  }

  return (out)
}



#________________________________________________________
# Bind tables by row after harmonizing column types, with one warning
# that lists every column whose type differed between inputs.
#________________________________________________________
bind_harmonized <- function (tbls, key, what) {

  type_of <- function (x) {
    if (is.numeric(x) || is.logical(x)) return ("numeric")
    if (is.factor(x) || is.character(x)) return ("categorical")
    class(x)[[1]]
  }

  conflicts <- character(0)

  for (col in setdiff(unique(unlist(lapply(tbls, names))), key)) {

    types <- unique(unlist(lapply(tbls, function (t) if (hasName(t, col)) type_of(t[[col]]))))
    if (length(types) < 2) next

    conflicts <- c(conflicts, sprintf("%s (%s)", col, paste(types, collapse = ", ")))

    for (k in seq_along(tbls))
      if (hasName(tbls[[k]], col))
        tbls[[k]][[col]] <- as.character(tbls[[k]][[col]])
  }

  if (length(conflicts) > 0)
    cli_warn(c(
      '!' = "Converted {length(conflicts)} {what}{?s} to categorical because the type differs between inputs:",
      setNames(conflicts, rep('*', length(conflicts))) ))

  dplyr::bind_rows(tbls)
}



#________________________________________________________
# Union of the inputs' metadata. A sample in several inputs keeps the
# first row; differing rows are reported in one warning.
#________________________________________________________
merge_metadata <- function (tbls) {

  md <- bind_harmonized(tbls, key = '.sample', what = "metadata field")

  if (as.logical(anyDuplicated(md[['.sample']]))) {

    rows <- do.call(paste, c(lapply(md, as.character), sep = "\x01"))
    dups <- unique(md[['.sample']][duplicated(md[['.sample']])])

    differ <- dups[vapply(dups, function (s) length(unique(rows[md[['.sample']] == s])) > 1, TRUE)]

    if (length(differ) > 0)
      cli_warn(c(
        '!' = "{length(differ)} sample{?s} {?has/have} different metadata in different inputs; keeping the first.",
        'i' = "Affected sample{?s}: {.val {differ}}." ))

    md <- md[!duplicated(md[['.sample']]),,drop=FALSE]
  }

  return (md)
}



#________________________________________________________
# Union of the inputs' taxonomy. An OTU in several inputs keeps its
# most detailed lineage (fewest NA ranks); conflicts are reported once.
#________________________________________________________
merge_taxonomy <- function (tbls) {

  rank_sets <- unique(lapply(tbls, function (t) names(t)[-1]))
  if (length(rank_sets) > 1)
    cli_warn(c(
      '!' = "Taxonomic ranks differ between inputs; missing ranks are NA.",
      '*' = "Ranks: {.val {unique(unlist(rank_sets))}}" ))

  tax <- bind_harmonized(tbls, key = '.otu', what = "taxonomic rank")

  if (as.logical(anyDuplicated(tax[['.otu']]))) {

    lineage <- rep("", nrow(tax))
    n_na    <- integer(nrow(tax))

    if (ncol(tax) > 1) {
      lineage <- do.call(paste, c(lapply(tax[-1], as.character), sep = "\x01"))
      n_na    <- rowSums(is.na(as.matrix(tax[,-1,drop=FALSE])))
    }

    n_conflict <- sum(tapply(lineage, tax[['.otu']], function (l) length(unique(l)) > 1))

    if (n_conflict > 0)
      cli_warn(c(
        '!' = "{n_conflict} OTU{?s} {?has/have} different taxonomy in different inputs.",
        'i' = "Keeping the most detailed lineage for each." ))

    tax <- tax[order(n_na, seq_len(nrow(tax))),,drop=FALSE]
    tax <- tax[!duplicated(tax[['.otu']]),,drop=FALSE]
    tax <- tax[order(match(tax[['.otu']], unique(unlist(lapply(tbls, `[[`, '.otu'))))),,drop=FALSE]
  }

  for (i in seq_len(ncol(tax))[-1])
    if (is.character(tax[[i]])) tax[[i]] <- as.factor(tax[[i]])

  return (tax)
}



#________________________________________________________
# Union of the inputs' sequences; NA for OTUs without one.
#________________________________________________________
merge_sequences <- function (seq_list, all_otus) {

  have <- !vapply(seq_list, is.null, TRUE)
  if (!any(have)) return (NULL)

  seqs <- unlist(unname(seq_list[have]))

  dups <- unique(names(seqs)[duplicated(names(seqs))])
  n_conflict <- sum(vapply(dups, function (otu) length(unique(seqs[names(seqs) == otu])) > 1, TRUE))

  seqs    <- seqs[!duplicated(names(seqs))]
  missing <- setdiff(all_otus, names(seqs[!is.na(seqs)]))

  if (n_conflict > 0)
    cli_warn("{n_conflict} OTU{?s} {?has/have} different sequences in different inputs; keeping the first.")

  if (length(missing) > 0)
    cli_warn(c(
      '!' = "{length(missing)} OTU{?s} {?has/have} no sequence; set to NA.",
      'i' = if (!all(have)) "Not all inputs include sequences." ))

  res <- setNames(rep(NA_character_, length(all_otus)), all_otus)
  res[names(seqs)] <- unname(seqs)

  return (res)
}
