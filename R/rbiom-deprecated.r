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
  
  validate_biom()
  
  
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
#' Use [otu_matrix()] instead.
#' @export
counts <- function (biom) {
  deprecate_warn("2.0.0", "counts()", "otu_matrix()")
  otu_matrix(biom = biom)
}


#' @name info-deprecated
#' @rdname rbiom-deprecated
#' @section \code{info}:
#' Use [biom_info()] instead.
#' @export
info <- function (biom) {
  deprecate_warn("2.0.0", "info()", "biom_info()")
  biom_info(biom = biom)
}


#' @name metadata-deprecated
#' @rdname rbiom-deprecated
#' @section \code{metadata}:
#' Use [sample_metadata()] instead.
#' @export
metadata <- function (biom, field = NULL, cleanup = FALSE) {
  deprecate_warn("2.0.0", "metadata()", "sample_metadata()")
  if (!missing(cleanup)) warning("`cleanup` is defunct")
  sample_metadata(biom = biom, field = field)
}


#' @name nsamples-deprecated
#' @rdname rbiom-deprecated
#' @section \code{nsamples}:
#' Use [n_samples()] instead.
#' @export
nsamples <- function (biom) {
  deprecate_warn("2.0.0", "nsamples()", "n_samples()")
  n_samples(biom = biom)
}


#' @name ntaxa-deprecated
#' @rdname rbiom-deprecated
#' @section \code{ntaxa}:
#' Use [n_otus()] instead.
#' @export
ntaxa <- function (biom) {
  deprecate_warn("2.0.0", "ntaxa()", "n_otus()")
  n_otus(biom = biom)
}


#' @name phylogeny-deprecated
#' @rdname rbiom-deprecated
#' @section \code{phylogeny}:
#' Use [otu_tree()] instead.
#' @export
phylogeny <- function (biom) {
  deprecate_warn("2.0.0", "phylogeny()", "otu_tree()")
  otu_tree(biom = biom)
}


#' @name rarefy-deprecated
#' @rdname rbiom-deprecated
#' @section \code{sample_rarefy}:
#' Use [otu_tree()] instead.
#' @export
rarefy <- function (biom, depth = NULL, seed = 0) {
  deprecate_warn("2.0.0", "rarefy()", "sample_rarefy()")
  sample_rarefy(biom = biom, depth = depth, seed = seed)
}


#' @name read.biom-deprecated
#' @rdname rbiom-deprecated
#' @section \code{read.biom}:
#' Use [read_biom()] instead.
#' @export
read.biom <- function (src, tree = "auto", prune = FALSE) {
  deprecate_warn("2.0.0", "read.biom()", "read_biom()")
  if (!missing(prune)) warning("`prune` argument is defunct")
  read_biom(src = src, tree = tree)
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
#' Use [sample_names()] instead.
#' @export
sample.names <- function (biom) {
  deprecate_warn("2.0.0", "sample.names()", "sample_names()")
  sample_names(biom = biom)
}


#' @name select-deprecated
#' @rdname rbiom-deprecated
#' @section \code{select}:
#' Use [subset()] instead.
#' @export
select <- function (biom, samples = NULL, nTop = NULL, nRandom = NULL, seed = 0) {
  deprecate_warn("2.0.0", "select()", "sample_select()")
  sample_select(biom = biom, samples = samples, nTop = nTop, nRandom = nRandom, seed = seed)
}


#' @name sequences-deprecated
#' @rdname rbiom-deprecated
#' @section \code{sequences}:
#' Use [otu_sequences()] instead.
#' @export
sequences <- function (biom) {
  deprecate_warn("2.0.0", "sequences()", "otu_sequences()")
  otu_sequences(biom = biom)
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
#' Use [otu_names()] instead.
#' @export
taxa.names <- function (biom) {
  deprecate_warn("2.0.0", "taxa.names()", "otu_names()")
  otu_names(biom = biom)
}


#' @name taxa.ranks-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxa.ranks}:
#' Use [taxa_ranks()] instead.
#' @export
taxa.ranks <- function (biom) {
  deprecate_warn("2.0.0", "taxa.ranks()", "taxa_ranks()")
  taxa_ranks(biom = biom)
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
    validate_biom(clone = FALSE)
    otu_taxonomy(biom) <- map
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
#' Use [otu_taxonomy()] instead.
#' @export
taxonomy <- function (biom, ranks = NULL, unc = "asis") {
  deprecate_warn("2.0.0", "taxonomy()", "otu_taxonomy()")
  otu_taxonomy(biom = biom, rank = NULL, unc = unc)
}


#' @name tips-deprecated
#' @rdname rbiom-deprecated
#' @section \code{tips}:
#' Use [tree_tips()] instead.
#' @export
tips <- function (x) {
  deprecate_warn("2.0.0", "tips()", "tree_tips()")
  tree_tips(x = x)
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
  write_fasta(seqs = seqs, file = outfile)
}


#' @name write.tree-deprecated
#' @rdname rbiom-deprecated
#' @section \code{write.tree}:
#' Use [write_tree()] instead.
#' @export
write.tree <- function (tree, file = NULL) {
  deprecate_warn("2.0.0", "write.tree()", "write_tree()")
  write_tree(tree = tree, file = file)
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
#' Use [as_percent()] instead.
#' @export
as.percent <- function (biom) {
  .Deprecated("as_percent")
  deprecate_warn("2.0.0", "as.percent()", "as_percent()")
  as_percent(biom = biom)
}


#' @name comments-deprecated
#' @rdname rbiom-deprecated
#' @section \code{comments}:
#' Use [biom_comment()] instead.
#' @export
comments <- function (biom) {
  deprecate_warn("2.0.0", "comments()", "biom_comment()")
  biom_comment(biom)
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
#' Use [has_tree()] instead.
#' @export
has.phylogeny <- function (biom) {
  deprecate_warn("2.0.0", "has.phylogeny()", "has_tree()")
  has_tree(biom = biom)
}


#' @name has.sequences-deprecated
#' @rdname rbiom-deprecated
#' @section \code{has.sequences}:
#' Use [has_sequences()] instead.
#' @export
has.sequences <- function (biom) {
  deprecate_warn("2.0.0", "has.sequences()", "has_sequences()")
  has_sequences(biom = biom)
}



#' @name id-deprecated
#' @rdname rbiom-deprecated
#' @section \code{id}:
#' Use [biom_id()] instead.
#' @export
id <- function (biom) {
  deprecate_warn("2.0.0", "id()", "biom_id()")
  biom_id(biom)
}



#' @name is.rarefied-deprecated
#' @rdname rbiom-deprecated
#' @section \code{is.rarefied}:
#' Use [is_rarefied()] instead.
#' @export
is.rarefied <- function (biom) {
  deprecate_warn("2.0.0", "is.rarefied()", "is_rarefied()")
  is_rarefied(biom = biom)
}


#' @name repair-deprecated
#' @rdname rbiom-deprecated
#' @section \code{repair}:
#' Use [biom_repair()] instead.
#' @export
repair <- function (biom) {
  deprecate_warn("2.0.0", "repair()", "biom_repair()")
  biom_repair(biom = biom)
}


#' @name sample_subset-deprecated
#' @rdname rbiom-deprecated
#' @section \code{sample_subset}:
#' Use [subset()] instead.
#' @export
sample_subset <- function (x, ...) {
  deprecate_warn("2.0.0", "sample_subset()", "subset()")
  subset(x, ...)
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


# #' @name stats.table-deprecated
# #' @rdname rbiom-deprecated
# #' @section \code{stats.table}:
# #' Use [biom_stats()] instead.
# #' @export
# stats.table <- function (...) {
#   .Deprecated("biom_stats")
#   biom_stats(...)
# }


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
#' Use [biom_comment()] instead.
#' @export
`comments<-` <- function (x, value) {
  deprecate_warn("2.0.0", "comments()", "biom_comment()")
  biom_comment(x) <- value
}


#' @name counts-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{counts-set}:
#' Use [otu_matrix()] instead.
#' @export
`counts<-` <- function (x, value) {
  deprecate_warn("2.0.0", "counts()", "otu_matrix()")
  otu_matrix(x) <- value
}


#' @name id-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{id-set}:
#' Use [biom_id()] instead.
#' @export
`id<-` <- function (x, value) {
  deprecate_warn("2.0.0", "id()", "biom_id()")
  biom_id(x) <- value
}


#' @name metadata-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{metadata-set}:
#' Use [sample_metadata()] instead.
#' @export
`metadata<-` <- function (x, value) {
  deprecate_warn("2.0.0", "metadata()", "sample_metadata()")
  sample_metadata(x) <- value
}


#' @name phylogeny-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{phylogeny-set}:
#' Use [otu_tree()] instead.
#' @export
`phylogeny<-` <- function (x, value) {
  deprecate_warn("2.0.0", "phylogeny()", "otu_tree()")
  otu_tree(x) <- value
}


#' @name sample.names-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{sample.names-set}:
#' Use [sample_names()] instead.
#' @export
`sample.names<-` <- function (x, value) {
  deprecate_warn("2.0.0", "sample.names()", "sample_names()")
  `sample_names<-`(x, value)
}


#' @name sequences-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{sequences-set}:
#' Use [otu_sequences()] instead.
#' @export
`sequences<-` <- function (x, value) {
  deprecate_warn("2.0.0", "sequences()", "otu_sequences()")
  otu_sequences(x) <- value
}


#' @name taxa.names-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxa.names-set}:
#' Use [otu_names()] instead.
#' @export
`taxa.names<-` <- function (x, value) {
  deprecate_warn("2.0.0", "taxa.names()", "otu_names()")
  `otu_names<-`(x, value)
}


#' @name taxa.ranks-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxa.ranks-set}:
#' Use [taxa_ranks()] instead.
#' @export
`taxa.ranks<-` <- function (x, value) {
  deprecate_warn("2.0.0", "taxa.ranks()", "taxa_ranks()")
  `taxa_ranks<-`(x, value)
}


#' @name taxonomy-set-deprecated
#' @rdname rbiom-deprecated
#' @section \code{taxonomy-set}:
#' Use [otu_taxonomy()] instead.
#' @export
`taxonomy<-` <- function (x, value) {
  deprecate_warn("2.0.0", "taxonomy()", "otu_taxonomy()")
  otu_taxonomy(x) <- value
}

