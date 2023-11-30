

#' Create an rbiom object.
#' 
#' @inherit read_biom return
#' 
#' @family biom
#' 
#' @param counts  The count data as a numeric matrix, where column names are
#'        the sample names and row names are the OTU names. An rbiom object or
#'        a filename/URL compatible with [read_biom()] is also acceptable.
#' 
#' @param metadata  A data.frame with sample names as the row names and 
#'        metadata names as the column names. Or a filename/URL with comma- or 
#'        tab-separated data.
#' 
#' @param taxonomy  A character matrix with OTU names as the row names. Or a 
#'        filename/URL with comma- or tab-separated data.
#' 
#' @param tree  A \code{phylo} object with tip labels matching the OTU names. 
#'        (E.g. from [read_tree()]). Or a filename/URL with newick formatted 
#'        data.
#' 
#' @param sequences  A named character vector of DNA sequences. Or a 
#'        filename/URL with fasta format data.
#'        (Currently not used by \code{rbiom}.)
#' 
#' @param id,comment  A character vector of length one with text of the user's
#'        choosing.
#'        
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- sample_select(hmp50, 1:5)
#'     
#'     # Save components to separate files.
#'     ct  <- write_counts(biom, tempfile())
#'     md  <- write_metadata(biom, tempfile())
#'     tax <- write_taxonomy(biom, tempfile())
#'     tre <- write_tree(biom, tempfile())
#'     fas <- write_fasta(biom, tempfile())
#'     
#'     # Peek at the file structures.
#'     cat(readLines(ct, n = 2L),  '', sep="\n")
#'     cat(readLines(md, n = 2L),  '', sep="\n")
#'     cat(readLines(tax, n = 2L), '', sep="\n")
#'     cat(readChar(tre, nchars = 50L), "\n\n")
#'     cat(readChar(fas, nchars = 50L), "\n\n")
#'     
#'     # Re-assemble the rbiom object.
#'     biom <- biom_build(ct, md, tax, tre, id = "New")
#'     print(biom)
#'     
#'     # Remove temporary files.
#'     invisible(file.remove(c(ct, md, tax, fas, tre)))

biom_build <- function (
    counts, metadata = NULL, taxonomy = NULL, tree = NULL, 
    sequences = NULL, id = NULL, comment = NULL ) {
  
  
  
  #________________________________________________________
  # Initialize the biom object based on counts.
  #________________________________________________________
  biom <- local({
    
    if (is(counts, 'rbiom'))
      return (counts)
    
    if (is_scalar_character(counts))
      return (read_biom(src = counts, tree = FALSE))
    
    
    #________________________________________________________
    # Coerce to matrix.
    #________________________________________________________
    if (!is.matrix(counts))
      counts <- tryCatch(
        expr  = as.matrix(counts),
        error = function (e) ({
          msg <- "Can't convert counts object of class '%s' to matrix.\n"
          stop(sprintf(msg, paste(collapse = " ", class(counts))), e) }) )
    
    stopifnot(is_string(typeof(counts), c("integer", "double")))
    
    
    
    #________________________________________________________
    # Construct the minimal rbiom object.
    #________________________________________________________
    return (structure(
      class = c('rbiom', "list"),
      .Data = list(
        counts   = slam::as.simple_triplet_matrix(counts),
        metadata = tibble(.sample = colnames(counts)),
        taxonomy = tibble(.otu    = rownames(counts)),
        info   = list(
          id                  = "",
          type                = "OTU table",
          format              = "1.0.0",
          format_url          = "http://biom-format.org",
          generated_by        = paste("rbiom", packageVersion("rbiom")),
          date                = strftime(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC"),
          matrix_type         = "sparse",
          matrix_element_type = ifelse(typeof(counts) == "integer", "int", "float"),
          shape               = dim(counts),
          comment             = "" ))))
    
  })
  
  
  #________________________________________________________
  # Integrate remaining rbiom components.
  #________________________________________________________
  if (!is.null(metadata))  sample_metadata(biom) <- metadata
  if (!is.null(taxonomy))  otu_taxonomy(biom)    <- taxonomy
  if (!is.null(tree))      otu_tree(biom)        <- tree
  if (!is.null(sequences)) otu_sequences(biom)   <- sequences
  if (!is.null(id))        biom_id(biom)         <- id
  if (!is.null(comment))   biom_comment(biom)    <- comment
  
  
  biom <- biom_repair(biom)
  
  return (biom)
}






#' Combine several rbiom objects into one.
#' 
#' WARNING: It is generally ill-advised to merge BIOM datasets, as OTUs
#' mappings are dependent on upstream clustering and are not equivalent
#' between BIOM files.
#' 
#' @inherit read_biom return
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
#'     b1 <- sample_select(hmp50, 1:4)
#'     b2 <- sample_select(hmp50, 5:8)
#'     
#'     biom <- biom_merge(b1, b2)
#'     print(biom)
#'     
#'     otu_tree(biom) <- otu_tree(hmp50)
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
  
  
  
  # return (structure(
  #   class = c('rbiom', "list"),
  #   .Data = list(
  #     counts    = counts,
  #     metadata  = metadata,
  #     taxonomy  = taxonomy,
  #     sequences = sequences,
  #     info      = list(id = "Merged BIOM", comment = "") )))
  
  
  biom <- biom_build(
    id        = "Merged BIOM",
    counts    = counts, 
    metadata  = metadata, 
    taxonomy  = taxonomy,
    sequences = sequences )
  
  return (biom)
}



