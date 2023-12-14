## rbiom-deprecated.r
#' @name rbiom-deprecated
#' @title Deprecated functions in package \pkg{rbiom}.
#' @description The functions listed below are deprecated and will be defunct in
#'   the near future. When possible, alternative functions with similar
#'   functionality are also mentioned. Help pages for deprecated functions are
#'   available at \code{help("<function>-deprecated")}.
#' @keywords internal
NULL



# Most functions have been renamed since last CRAN version. Reasons:
# 
#  * Use underscores instead of periods.
#  * Better description of their purpose.
#  * Combine near-duplicates.



#________________________________________________________
# Changes impacting users of CRAN's rbiom 1.0.3
#________________________________________________________

#' @name alpha.div-deprecated
#' @rdname rbiom-deprecated
#' @keywords internal
#' @section \code{alpha.div}:
#' Use [adiv_matrix()] or [adiv_table()] instead.
#' @export
alpha.div <- function (biom, rarefy = FALSE) {
  
  deprecate_warn(
    when    = "2.0.0",
    what    = "alpha.div()",
    details = "Please use `adiv_matrix()` or `adiv_table()` instead." )
  
  biom <- as_rbiom(biom)
  
  
  # Log intervals until rLvl/2, then even intervals until rLvl*2
  if (eq(rarefy, "multi"))
    rarefy <- local({
      rLvl  <- rare_suggest(biom)
      rLvls <- 10 ** (c(1,3,5,7,8,9) * (log10(rLvl / 2) / 10))
      rLvls <- c(rLvls, seq(from = rLvl / 2, to = rLvl * 2, length.out = 5))
      rLvls <- floor(rLvls)
      return (rLvls)
    })
  
  
  stopifnot(is_scalar_logical(rarefy) || is_integerish(rarefy))
  
  res <- NULL
  for (i in rarefy) {
    
    df <- adiv_matrix(biom = biom, rarefy = i) %>%
      as_tibble(x, rownames = 'Sample') %>%
      as.data.frame()
    
    if (length(rarefy == 1))
      rownames(df) <- df[['Sample']]
    
    res <- dplyr::bind_rows(res, df)
  }
  
  return (res)
}


# Depending on its arguments, beta.div() returned either a dist object or a 
# data.frame. Now it's split into two functions - bdiv_distmat() which always 
# returns a dist, and bdiv_table() which always returns a data.frame.

#' @name beta.div-deprecated
#' @rdname rbiom-deprecated
#' @section \code{beta.div}:
#' Use [bdiv_table()] or [bdiv_distmat()] instead.
#' @export
beta.div <- function (biom, method, weighted = TRUE, tree = NULL, long = FALSE, md = FALSE) {
  
  if (isTRUE(long) || !isFALSE(md)) {
    
    deprecate_warn(
      when    = "2.0.0",
      what    = "beta.div()",
      details = "Please use `bdiv_table()` instead for generating a data.frame." )
    
    if (isTRUE(md)) md <- ".all"
    
    bdiv_table(biom = biom, bdiv = method, weighted = weighted, tree = tree, md = md)
    
  } else {
    
    deprecate_warn(
      when    = "2.0.0",
      what    = "beta.div()",
      details = "Please use `bdiv_distmat()` instead for generating a distance matrix." )
    
    bdiv_distmat(biom = biom, bdiv = method, weighted = weighted, tree = tree)
  }
}


#' @name counts-deprecated
#' @rdname rbiom-deprecated
#' @section \code{counts}:
#' Use `$counts` instead.
#' @export
counts <- function (biom) {
  deprecate_warn("2.0.0", "counts()", "otu_matrix()")
  as.matrix(as_rbiom(biom)$counts)
}


#' @name info-deprecated
#' @rdname rbiom-deprecated
#' @section \code{info}:
#' Use `biom$id`, `biom$comment`, etc instead.
#' @export
info <- function (biom) {
  deprecate_warn("2.0.0", "info()", "Use `biom$id`, `biom$comment`, etc instead.")
  list(
    id                  = biom$id, 
    comment             = biom$comment, 
    date                = biom$date, 
    generated_by        = biom$generated_by,
    format              = "1.0.0", 
    type                = "OTU table",
    format_url          = "http://biom-format.org",
    matrix_type         = "sparse",
    matrix_element_type = ifelse(all(biom$counts[['v']] %% 1 == 0), "int", "float"),
    shape               = dim(biom$counts) )
}


#' @name metadata-deprecated
#' @rdname rbiom-deprecated
#' @section \code{metadata}:
#' Use `biom$metadata` or `pull(biom, field)` instead.
#' @export
metadata <- function (biom, field = NULL, cleanup = FALSE) {
  deprecate_warn("2.0.0", "metadata()", "Use `biom$metadata` or `pull(biom, field)` instead.")
  if (!missing(cleanup)) warning("`cleanup` is defunct")
  
  biom <- as_rbiom(biom)
  if (is.null(field)) biom$metadata else pull(biom = biom, field = field)
}


#' @name nsamples-deprecated
#' @rdname rbiom-deprecated
#' @section \code{nsamples}:
#' Use `biom$n_samples` instead.
#' @export
nsamples <- function (biom) {
  deprecate_warn("2.0.0", "nsamples()", "biom$n_samples")
  as_rbiom(biom)$n_samples
}


#' @name ntaxa-deprecated
#' @rdname rbiom-deprecated
#' @section \code{ntaxa}:
#' Use `biom$n_otus` instead.
#' @export
ntaxa <- function (biom) {
  deprecate_warn("2.0.0", "ntaxa()", "biom$n_otus")
  as_rbiom(biom)$n_otus
}


#' @name phylogeny-deprecated
#' @rdname rbiom-deprecated
#' @section \code{phylogeny}:
#' Use `biom$tree` instead.
#' @export
phylogeny <- function (biom) {
  deprecate_warn("2.0.0", "phylogeny()", "biom$tree")
  as_rbiom(biom)$tree
}


#' @name read.biom-deprecated
#' @rdname rbiom-deprecated
#' @section \code{read.biom}:
#' Use [read_biom()] instead.
#' @export
read.biom <- function (src, tree = "auto", prune = FALSE) {
  
  deprecate_warn("2.0.0", "read.biom()", "as_rbiom()")
  if (!missing(prune)) warning("`prune` argument is defunct")
  
  biom <- as_rbiom(biom = src)
  
  if (isTRUE(tree) && is.null(biom$tree))
    cli_abort("Dataset does not include a phylogenetic tree.")
  
  if (isTRUE(tree))  tree <- "auto"
  if (isFALSE(tree)) tree <- NULL
  
  if (!identical(tree, "auto"))
    biom$tree <- tree
  
  return (biom)
}


#' @name read.fasta-deprecated
#' @rdname rbiom-deprecated
#' @section \code{read.fasta}:
#' Use [read_fasta()] instead.
#' @export
read.fasta <- function (file, ids = NULL) {
  deprecate_warn("2.0.0", "read.fasta()", "read_fasta()")
  read_fasta(file = file, ids = ids)
}


#' @name read.tree-deprecated
#' @rdname rbiom-deprecated
#' @section \code{read.tree}:
#' Use [read_tree()] instead.
#' @export
read.tree <- function (src) {
  deprecate_warn("2.0.0", "read.tree()", "read_tree()")
  read_tree(src = src)
}


#' @name sample.names-deprecated
#' @rdname rbiom-deprecated
#' @section \code{sample.names}:
#' Use `biom$samples` instead.
#' @export
sample.names <- function (biom) {
  deprecate_warn("2.0.0", "sample.names()", "biom$samples")
  as_rbiom(biom)$samples
}


#' @name select-deprecated
#' @rdname rbiom-deprecated
#' @section \code{select}:
#' Use [slice()] instead.
#' @export
select.rbiom <- function (biom, samples = NULL, nTop = NULL, nRandom = NULL, seed = 0) {
  
  deprecate_warn("2.0.0", "select()", "slice()")
  biom <- biom$clone()
  
  if (!is.null(samples)) biom$counts <- biom$counts[,samples]
  if (!is.null(nTop))    biom$counts <- biom$counts[,names(head(sample_sums(biom), nTop))]
  if (!is.null(nRandom)) {
    set.seed(seed)
    stopifnot(nRandom <= biom$n_samples)
    biom$counts <- biom$counts[,sample(seq_len(ncol(biom$counts)), nRandom)]
  }
  
  return (biom)
}


#' @name sequences-deprecated
#' @rdname rbiom-deprecated
#' @section \code{sequences}:
#' Use `biom$sequences` instead.
#' @export
sequences <- function (biom) {
  deprecate_warn("2.0.0", "sequences()", "biom$sequences")
  as_rbiom(biom)$tree
}


#' @name subtree-deprecated
#' @rdname rbiom-deprecated
#' @section \code{subtree}:
#' Use [tree_subset()] instead.
#' @export
subtree <- function (tree, tips) {
  deprecate_warn("2.0.0", "subtree()", "tree_subset()")
  tree_subset(tree = tree, tips = tips)
}


#' @name taxa.names-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxa.names}:
#' Use `biom$otus` instead.
#' @export
taxa.names <- function (biom) {
  deprecate_warn("2.0.0", "taxa.names()", "biom$otus")
  as_rbiom(biom)$otus
}


#' @name taxa.ranks-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxa.ranks}:
#' Use `biom$ranks` instead.
#' @export
taxa.ranks <- function (biom) {
  deprecate_warn("2.0.0", "taxa.ranks()", "biom$ranks")
  as_rbiom(biom)$ranks
}


# Depending on its arguments, taxa.rollup() returned either a matrix or a 
# data.frame. Now it's split into two functions - taxa_matrix() which always 
# returns a matrix, and taxa_table() which always returns a data.frame.

#' @name taxa.rollup-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxa.rollup}:
#' Use [taxa_table()] [taxa_matrix()] instead.
#' @export
taxa.rollup <- function (
    biom, rank = 'OTU', map = NULL, lineage = FALSE, 
    sparse = FALSE, taxa = NULL, long = FALSE, md = FALSE) {
  
  
  if (!is.null(map)) {
    biom <- as_rbiom(biom)$clone()
    biom$taxonomy <- map
  }
  
  
  if (isTRUE(long) || !isFALSE(md)) {
    
    deprecate_warn(
      when    = "2.0.0",
      what    = "taxa.rollup()",
      details = "Please use `taxa_table()` instead for generating a data.frame." )
    
    if (isTRUE(md)) md <- ".all"
    
    taxa_table(
      biom    = biom, 
      rank    = rank, 
      taxa    = taxa, 
      lineage = lineage, 
      md      = md )
    
  } else {
    
    deprecate_warn(
      when    = "2.0.0",
      what    = "beta.div()",
      details = "Please use `taxa_matrix()` instead for generating a matrix." )
    
    taxa_matrix(
      biom    = biom, 
      rank    = rank, 
      taxa    = taxa, 
      lineage = lineage, 
      sparse  = sparse )
  }
}


#' @name taxonomy-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxonomy}:
#' Use `$taxonomy` instead.
#' @export
taxonomy <- function (biom, ranks = NULL, unc = "asis") {
  deprecate_warn("2.0.0", "taxonomy()", "Use `biom$taxonomy` or `taxa_map(biom, rank)` instead.")
  biom <- as_rbiom(biom)
  if (is.null(ranks)) biom$taxonomy else taxa_map(biom, ranks)
}


#' @name tips-deprecated
#' @rdname rbiom-deprecated
#' @section \code{tips}:
#' Use `tree$tip.label` instead.
#' @export
tips <- function (x) {
  deprecate_warn("2.0.0", "tips()", "tree$tip.label")
  validate_tree("x")
  x$tip.label
}


#' @name unifrac-deprecated
#' @rdname rbiom-deprecated
#' @section \code{unifrac}:
#' Use [bdiv_distmat()] or [bdiv_table()] instead.
#' @export
unifrac <- function (biom, weighted=TRUE, tree=NULL) {
  deprecate_soft("2.0.0", "unifrac()", "bdiv_distmat()")
  bdiv_distmat(biom = biom, bdiv = "unifrac", weighted = weighted, tree = tree)
}


#' @name write.biom-deprecated
#' @rdname rbiom-deprecated
#' @section \code{write.biom}:
#' Use [write_biom()] instead.
#' @export
write.biom <- function (biom, file, format="json") {
  deprecate_warn("2.0.0", "write.biom()", "write_biom()")
  write_biom(biom = biom, file = file, format = format)
}


#' @name write.fasta-deprecated
#' @rdname rbiom-deprecated
#' @section \code{write.fasta}:
#' Use [write_fasta()] instead.
#' @export
write.fasta <- function (seqs, outfile = NULL) {
  deprecate_warn("2.0.0", "write.fasta()", "write_fasta()")
  write_fasta(biom = seqs, file = outfile)
}


#' @name write.tree-deprecated
#' @rdname rbiom-deprecated
#' @section \code{write.tree}:
#' Use [write_tree()] instead.
#' @export
write.tree <- function (tree, file = NULL) {
  deprecate_warn("2.0.0", "write.tree()", "write_tree()")
  write_tree(biom = tree, file = file)
}


#' @name write.xlsx-deprecated
#' @rdname rbiom-deprecated
#' @section \code{write.xlsx}:
#' Use [write_xlsx()] instead.
#' @export
write.xlsx <- function (biom, outfile, depth = NULL, seed = 0) {
  deprecate_warn("2.0.0", "write.xlsx()", "write_xlsx()")
  write_xlsx(biom = biom, file = outfile, depth = depth, seed = seed)
}





#________________________________________________________
# Changes impacting users of github's development version
#________________________________________________________


#' @name as.percent-deprecated
#' @rdname rbiom-deprecated
#' @section \code{as.percent}:
#' Use `biom$counts %<>% rescale_cols()` instead.
#' @export
as.percent <- function (biom) {
  
  deprecate_warn(
    when    = "2.0.0",
    what    = "as.percent()",
    details = "Please use `biom$counts %<>% rescale_cols()` instead." )
  biom <- biom$clone()
  
  biom$counts %<>% rescale_cols()
  return (biom)
}


#' @name comments-deprecated
#' @rdname rbiom-deprecated
#' @section \code{comments}:
#' Use `biom$comment` instead.
#' @export
comments <- function (biom) {
  deprecate_warn("2.0.0", "comments()", "biom$comment")
  as_rbiom(biom)$comment
}


#' @name depth-deprecated
#' @rdname rbiom-deprecated
#' @section \code{depth}:
#' Use [sample_sums()] instead.
#' @export
depth <- function (biom) {
  deprecate_warn("2.0.0", "depth()", "sample_sums()")
  sample_sums(biom = biom)
}


#' @name depths_barplot-deprecated
#' @rdname rbiom-deprecated
#' @section \code{depths_barplot}:
#' Use [rare_barplot()] instead.
#' @export
depths_barplot <- function (biom, rline = TRUE, counts = TRUE, labels = TRUE, trans = "log10", ...) {
  deprecate_warn("2.0.0", "depths_barplot()", "rare_barplot()")
  rare_barplot(biom = biom, rline = rline, counts = counts, labels = labels, trans = trans, ...)
}


#' @name has.phylogeny-deprecated
#' @rdname rbiom-deprecated
#' @section \code{has.phylogeny}:
#' Use `!is.null(biom$tree)` instead.
#' @export
has.phylogeny <- function (biom) {
  deprecate_warn("2.0.0", "has.phylogeny()", "!is.null(biom$tree)")
  !is.null(as_rbiom(biom)$tree)
}


#' @name has.sequences-deprecated
#' @rdname rbiom-deprecated
#' @section \code{has.sequences}:
#' Use `!is.null(biom$sequences)` instead.
#' @export
has.sequences <- function (biom) {
  deprecate_warn("2.0.0", "has.sequences()", "!is.null(biom$sequences)")
  !is.null(as_rbiom(biom)$sequences)
}



#' @name id-deprecated
#' @rdname rbiom-deprecated
#' @section \code{id}:
#' Use `biom$id` instead.
#' @export
id <- function (biom) {
  deprecate_warn("2.0.0", "id()", "biom$id")
  as_rbiom(biom)$id
}



#' @name is.rarefied-deprecated
#' @rdname rbiom-deprecated
#' @section \code{is.rarefied}:
#' Use `!is.null(biom$depth)` instead.
#' @export
is.rarefied <- function (biom) {
  deprecate_warn("2.0.0", "is.rarefied()", "!is.null(biom$depth)")
  !is.null(as_rbiom(biom)$depth)
}


#' @name repair-deprecated
#' @rdname rbiom-deprecated
#' @section \code{repair}:
#' Use [as_rbiom(as.list(biom))] instead.
#' @export
repair <- function (biom) {
  deprecate_warn("2.0.0", "repair()", "as_rbiom(as.list(biom))")
  as_rbiom(as.list(biom))
}


#' @name sample_subset-deprecated
#' @rdname rbiom-deprecated
#' @section \code{sample_subset}:
#' Use `biom$metadata %<>% base::subset()` instead.
#' @export
sample_subset <- function (x, ...) {
  deprecate_warn("2.0.0", "sample_subset()", "subset()")
  x <- x$clone()
  x$metadata %<>% base::subset(...)
  return (x)
}


#' @name sample.sums-deprecated
#' @rdname rbiom-deprecated
#' @section \code{sample.sums}:
#' Use [sample_sums()] or [adiv_table()] instead.
#' @export
sample.sums <- function (biom, long = FALSE, md = FALSE) {
  
  if (isFALSE(long) && isFALSE(md)) {
    
    deprecate_warn("2.0.0", "sample.sums()", "sample_sums()")
    sample_sums(biom = biom)
    
  } else {
    
    deprecate_warn("2.0.0", "sample.sums()", "adiv_table()")
    
    if (isTRUE(md))  md <- ".all"
    if (isFALSE(md)) md <- NULL
    
    adiv_table(biom = biom, md = md) %>% 
      dplyr::mutate(.keep = "none", Sample = .sample, Reads = .depth) %>%
      as.data.frame()
  }
}


#' @name stats.table-deprecated
#' @rdname rbiom-deprecated
#' @section \code{stats.table}:
#' Use [stats_table()] instead.
#' @export
stats.table <- function (...) {
  deprecate_warn("2.0.0", "stats.table()", "stats_table()")
  stats_table(...)
}


#' @name taxa.means-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxa.means}:
#' Use [taxa_means()] instead.
#' @export
taxa.means <- function (biom, rank = NULL) {
  deprecate_warn("2.0.0", "taxa.means()", "taxa_means()")
  taxa_means(biom = biom, rank = rank)
}


#' @name taxa.sums-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxa.sums}:
#' Use [taxa_sums()] instead.
#' @export
taxa.sums <- function (biom, rank = NULL) {
  deprecate_warn("2.0.0", "taxa.sums()", "taxa_sums()")
  taxa_sums(biom = biom, rank = rank)
}


#' @name top.taxa-deprecated
#' @rdname rbiom-deprecated
#' @section \code{top.taxa}:
#' Use [taxa_sums()] instead.
#' @export
top.taxa <- function (biom, rank = 'OTU', n = Inf) {
  deprecate_warn("2.0.0", "top.taxa()", "taxa_sums()")
  names(head(taxa_sums(biom, rank), n))
}


#' @name top_taxa-deprecated
#' @rdname rbiom-deprecated
#' @section \code{top_taxa}:
#' Use [taxa_sums()] instead.
#' @export
top_taxa <- function (biom, rank = 'OTU', n = Inf) {
  deprecate_warn("2.0.0", "top_taxa()", "taxa_sums()")
  names(head(taxa_sums(biom, rank), n))
}


#' @name comments-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{comments-set}:
#' Use [$comment] instead.
#' @export
`comments<-` <- function (x, value) {
  deprecate_warn("2.0.0", "comments()", "$comment")
  x$comment <- value
}


#' @name counts-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{counts-set}:
#' Use [$counts] instead.
#' @export
`counts<-` <- function (x, value) {
  deprecate_warn("2.0.0", "counts()", "$counts")
  x$counts <- value
}


#' @name id-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{id-set}:
#' Use [$id] instead.
#' @export
`id<-` <- function (x, value) {
  deprecate_warn("2.0.0", "id()", "$id")
  x$id <- value
}


#' @name metadata-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{metadata-set}:
#' Use [$metadata] instead.
#' @export
`metadata<-` <- function (x, value) {
  deprecate_warn("2.0.0", "metadata()", "$metadata")
  x$metadata <- value
}


#' @name phylogeny-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{phylogeny-set}:
#' Use [$tree] instead.
#' @export
`phylogeny<-` <- function (x, value) {
  deprecate_warn("2.0.0", "phylogeny()", "$tree")
  x$tree <- value
}


#' @name sample.names-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{sample.names-set}:
#' @export
`sample.names<-` <- function (x, value) {
  deprecate_warn("2.0.0", "sample.names()")
  
  stopifnot(isa(x, "rbiom"))
  stopifnot(is.character(value))
  stopifnot(length(value) == ncol(x$counts))
  
  b <- as.list(x)
  colnames(b$counts)      <- value
  b$metadata[['.sample']] <- value
  if (!is.null(b$sequences))
    names(b$sequences) <- value
  
  as_rbiom(b)
}


#' @name sequences-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{sequences-set}:
#' Use [otu_sequences()] instead.
#' @export
`sequences<-` <- function (x, value) {
  deprecate_warn("2.0.0", "sequences()", "$sequences")
  x$sequences <- value
}


#' @name taxa.names-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxa.names-set}:
#' @export
`taxa.names<-` <- function (x, value) {
  deprecate_warn("2.0.0", "taxa.names()")
  
  stopifnot(isa(x, "rbiom"))
  stopifnot(is.character(value))
  stopifnot(length(value) == nrow(x$counts))
  
  b <- as.list(x)
  if (!is.null(b$tree))
    b$tree$tip.label <- value[match(rownames(b$counts), b$tree$tip.label)]
  rownames(b$counts)   <- value
  b$taxonomy[['.otu']] <- value
  
  as_rbiom(b)
}


#' @name taxa.ranks-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxa.ranks-set}:
#' @export
`taxa.ranks<-` <- function (x, value) {
  deprecate_warn("2.0.0", "taxa.ranks()")
  
  df <- x$taxonomy
  colnames(df) <- value
  x$taxonomy <- df
  
  return (x)
}


#' @name taxonomy-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxonomy-set}:
#' Use [$taxonomy] instead.
#' @export
`taxonomy<-` <- function (x, value) {
  deprecate_warn("2.0.0", "taxonomy()", "$taxonomy")
  x$taxonomy <- value
}

