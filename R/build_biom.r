

#' Create a BIOM object.
#' 
#' @param counts  The count data as a numeric matrix, where column names are
#'        the sample names and row names are the taxa IDs. A BIOM object or
#'        a filename/URL compatible with \code{read_biom()} is also acceptable.
#' 
#' @param metadata  A data.frame with sample IDs as the row names and metadata 
#'        names as the column names. Or a filename/URL with comma- or 
#'        tab-separated data.
#' 
#' @param taxonomy  A character matrix with taxa IDs as the row names. Or a 
#'        filename/URL with comma- or tab-separated data.
#' 
#' @param phylogeny  A \code{phylo} object with tip labels matching the taxa 
#'        IDs. (E.g. from \code{rbiom::read_tree()}). Or a filename/URL with 
#'        newick formatted data.
#' 
#' @param sequences  A named character vector of DNA sequences. Or a 
#'        filename/URL with fasta format data.
#'        (Currently not used by \code{rbiom}.)
#' 
#' @param id,comments  A character vector of length one with text of the user's
#'        choosing.
#'        
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- select(hmp50, 1:5)
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
#'     # Re-assemble the BIOM object.
#'     biom <- build_biom(ct, md, tax, tre, id = "New BIOM")
#'     print(biom)
#'     
#'     # Remove temporary files.
#'     invisible(file.remove(c(ct, md, tax, fas, tre)))
#'     
#'

build_biom <- function (
    counts, metadata = NULL, taxonomy = NULL, phylogeny = NULL, 
    sequences = NULL, id = NULL, comments = NULL ) {
  
  
  
  #________________________________________________________
  # Initialize the biom object based on counts.
  #________________________________________________________
  biom <- local({
    
    if (is(counts, 'BIOM'))
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
    # Construct the minimal BIOM object.
    #________________________________________________________
    return (structure(
      class = c("BIOM", "list"),
      .Data = list(
        counts = slam::as.simple_triplet_matrix(counts),
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
  # Integrate remaining BIOM components.
  #________________________________________________________
  if (!missing(metadata))  rbiom::metadata(biom)  <- metadata
  if (!missing(taxonomy))  rbiom::taxonomy(biom)  <- taxonomy
  if (!missing(phylogeny)) rbiom::phylogeny(biom) <- phylogeny
  if (!missing(sequences)) rbiom::sequences(biom) <- sequences
  if (!missing(id))        rbiom::id(biom)        <- id
  if (!missing(comments))  rbiom::comments(biom)  <- comments
  
  
  return (biom)
}
