
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
#' @param ...  Any number of rbiom objects (e.g. from [read_biom()]), lists of
#'        rbiom objects, or valid arguments to the \code{src} parameter of 
#'        [read_biom()] (for instance file names).
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

biom_merge <- function (...) { 
  
  dots      <- list(...)
  biom_list <- lapply(dots, function (dot) {
      
    if (is(dot, 'rbiom'))      return (dot)
    if (length(dot) > 1)       return (do.call(biom_merge, dot))
    if (length(dot) < 1)       return (NULL)
    if (is(dot[[1]], 'rbiom')) return (dot[[1]])
    if (is.character(dot))     return (read_biom(src = dot))
    stop("Unknown argument to biom_merge(): ", dot)
    
  })
  biom_list <- biom_list[sapply(biom_list, is, 'rbiom')]
  
  if (length(biom_list) == 0) stop("No BIOM datasets provided to biom_merge().")
  if (length(biom_list) == 1) return (biom_list[[1]]) 
  
  
  samples <- do.call(c, lapply(biom_list, function (x) { colnames(x[['counts']]) }))
  otus    <- do.call(c, lapply(biom_list, function (x) { rownames(x[['counts']]) }))
  
  if (any(duplicated(samples))) stop("Sample names are not unique among BIOM datasets.")
  if (!any(duplicated(otus)))   warning("No overlapping OTU names. Likely incompatible datasets.")
  
  otus <- unique(otus)
  
  
  counts <- slam::simple_triplet_matrix(
    i        = match(do.call(c, lapply(biom_list, function (b) { rownames(b$counts)[b$counts$i] })), otus),
    j        = match(do.call(c, lapply(biom_list, function (b) { colnames(b$counts)[b$counts$j] })), samples),
    v        = do.call(c, lapply(biom_list, function (b) { b$counts$v })),
    nrow     = length(otus),
    ncol     = length(samples),
    dimnames = list(otus, samples) )
  
  
  metadata <- dplyr::bind_rows(lapply(biom_list, `[[`, 'metadata'))
  taxonomy <- dplyr::bind_rows(lapply(biom_list, `[[`, 'taxonomy'))
  
  if (ncol(taxonomy) > 1) {
    taxstrs <- apply(taxonomy, 1L, paste, collapse = "; ")
    for (otu in otus)
      if (length(strs <- unique(taxstrs[which(names(taxstrs) == otu)])) > 1)
        warning("OTU '", otu, "' has multiple taxonomic mappings:", paste("\n  ", strs))
  }
  
  taxonomy <- taxonomy[!duplicated(taxonomy[['.otu']]),]
  
  
  sequences <- do.call(c, lapply(biom_list, `[[`, 'sequences'))
  
  
  rbiom$new(
    id        = "Merged BIOM",
    counts    = counts, 
    metadata  = metadata, 
    taxonomy  = taxonomy,
    sequences = sequences )
}



