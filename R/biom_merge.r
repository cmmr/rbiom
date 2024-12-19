
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
#' @param ...  Any number of rbiom objects (e.g. from [as_rbiom()]), lists of
#'        rbiom objects, or valid arguments to the \code{biom} parameter of 
#'        [as_rbiom()] (for instance file names).
#'        
#' @param metadata,taxonomy,tree,sequences,id,comment Replace the corresponding 
#'        data in the merged rbiom object with these values. Set to `NULL` to
#'        not inherit a particular component. The default, `NA`, will attempt
#'        to create the component based on `...` values. The merged 
#'        phylogenetic tree cannot be inferred.
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
#'     biom$metadata <- hmp50$metadata
#'     print(biom)

biom_merge <- function (..., metadata = NA, taxonomy = NA, tree = NULL, sequences = NA, id = NA, comment = NA) { 
  
  dots      <- list(...)
  biom_list <- lapply(dots, function (dot) {
      
    if (inherits(dot, 'rbiom'))         return (dot)
    if (length(dot) > 1)          return (do.call(biom_merge, dot))
    if (length(dot) < 1)          return (NULL)
    if (inherits(dot[[1]], 'rbiom'))    return (dot[[1]])
    if (is_scalar_character(dot)) return (read_biom(src = dot))
    cli::cli_abort("Unknown argument to biom_merge(): {.typeof {dot}}")
    
  })
  biom_list <- biom_list[sapply(biom_list, inherits, 'rbiom')]
  
  if (length(biom_list) == 0) cli::cli_abort("No BIOM datasets provided to biom_merge().")
  if (length(biom_list) == 1) return (biom_list[[1]]) 
  
  
  samples <- do.call(c, lapply(biom_list, function (x) { colnames(x[['counts']]) }))
  otus    <- do.call(c, lapply(biom_list, function (x) { rownames(x[['counts']]) }))
  
  if (length(dups <- unique(samples[duplicated(samples)])) > 0)
    cli::cli_abort("Sample names are not unique among BIOM datasets: {.val {dups}}")
  
  if (!any(duplicated(otus)))
    cli::cli_warn("No overlapping OTU names. Likely incompatible datasets.")
  
  otus <- unique(otus)
  
  
  counts <- simple_triplet_matrix(
    i        = match(do.call(c, lapply(biom_list, function (b) { rownames(b$counts)[b$counts$i] })), otus),
    j        = match(do.call(c, lapply(biom_list, function (b) { colnames(b$counts)[b$counts$j] })), samples),
    v        = do.call(c, lapply(biom_list, function (b) { b$counts$v })),
    nrow     = length(otus),
    ncol     = length(samples),
    dimnames = list(otus, samples) )
  
  
  if (is.na(metadata))  metadata  <- dplyr::bind_rows(lapply(biom_list, `[[`, 'metadata'))
  if (is.na(taxonomy))  taxonomy  <- dplyr::bind_rows(lapply(biom_list, `[[`, 'taxonomy'))
  if (is.na(sequences)) sequences <- do.call(c,       lapply(biom_list, `[[`, 'sequences'))
  
  if (is.na(id)) {
    id <- unique(do.call(c, lapply(biom_list, `[[`, 'id')))
    id <- setdiff(id, 'Untitled Dataset')
    if (length(id) != 1) id <- "Merged BIOM"
  }
  
  if (is.na(comment)) {
    comment <- unique(do.call(c, lapply(biom_list, `[[`, 'comment')))
    comment <- setdiff(comment, '')
    if (length(comment) != 1) comment <- NULL
  }
  
  
  rbiom$new(
    id        = id,
    comment   = comment,
    counts    = counts, 
    metadata  = metadata, 
    taxonomy  = taxonomy,
    sequences = sequences,
    tree      = tree )
}



