#' Generate a data.frame of taxa abundance for each sample at the specified taxonomic rank.
#' 
#' @name taxa_table
#'
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'        object, as returned from \link{read_biom}. For matrices, the rows and 
#'        columns are assumed to be the taxa and samples, respectively.
#'        
#' @param rank  The taxonomic rank. E.g. \code{rank = "OTU"}, 
#'        \code{rank = "Phylum"}, etc. May also be given numerically: 0 for
#'        OTU, 1 for the highest level (i.e. Kingdom), and extending to the number
#'        of taxonomic ranks encoded in the original biom file. Run 
#'        \code{taxa_ranks(biom)} to fetch the names of all available ranks.
#'        
#' @param taxa   Limit the number of taxa returned in the matrix. Depending on 
#'        the value type, a different filter is applied.
#'        \describe{
#'          \item{\code{NULL}}{ Retain all taxa. (The default.) }
#'          \item{Integer >= 1}{ The N most abundant taxa are retained. Example:
#'            \code{taxa_table(biom, 'Genus', 5)} returned the five most
#'            abundant genera. }
#'          \item{Numeric between 0 and 1}{ Taxa are retained if their abundance 
#'            is greater than or equal to N. Example: 
#'            \code{taxa_table(biom, 'Phylum', 0.1)} returns all phyla with a
#'            relative abundance of at least 10 \%. }
#'          \item{Character vector}{ Only these taxa names are retained. Example: 
#'            \code{taxa_table(biom, 'Phylum', c("Firmicutes", "Bacteroidetes"))}
#'            will only retain those two phyla. }
#'        }
#'        
#' @param map  A character matrix defining the value that each taxa IDs is
#'        assigned for each taxonomic rank. If \code{map=NULL} and \code{biom} is a
#'        \code{BIOM} class object, the map will be automatically loaded from 
#'        \code{biom$taxonomy}. \code{map} must not be \code{null} when \code{biom}
#'        is a \code{matrix} or \code{simple_triplet_matrix}. See the example below
#'        for an example of \code{map}'s structure.
#'        
#' @param lineage  Include all ranks in the name of the taxa. For instance,
#'        setting to \code{TRUE} will produce 
#'        \code{Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales}. 
#'        Whereas setting to \code{FALSE} (the default) will return simply
#'        \code{Coriobacteriales}. You want to set this to TRUE if you have
#'        genus names (such as \emph{Incertae_Sedis}) that map to multiple higher
#'        level ranks.
#'        
#' @param md  Include metadata in the output data frame? Options are: 
#'        \describe{
#'          \item{\code{FALSE}}{ Don't include metadata. (Default) }
#'          \item{\code{TRUE}}{ Include all metadata. }
#'          \item{\emph{character vector}}{ Include only the specified metadata columns. }
#'        }
#'        
#' @param safe  Should autogenerated columns be prefixed with a "." to avoid conflicting
#'        with metadata column names? (Default: FALSE) \code{safe=NA} can be used with
#'        \code{long=FALSE} to omit the 'Sample' column entirely, leaving the sample
#'        names only as the row names.
#'        
#' @return A numeric matrix with samples as column names, and taxonomic
#'         identifiers as row names. Or a data.frame if metadata is requested 
#'         with \code{md}.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     taxa_ranks(hmp50)
#'     
#'     phyla <- taxa_table(hmp50, 'Phylum')
#'     phyla[1:4,1:6]
#'     
#'     # Custom matrices should be formatted like so:
#'     counts <- counts(hmp50)
#'     map    <- taxonomy(hmp50)
#'     
#'     counts[1:3,1:6]
#'     map[1:3,1:4]
#'     
#'     phyla <- taxa_table(counts, 'Phylum', map=map)
#'     phyla[1:3,1:6]
#'

taxa_table <- function (
    biom, rank = 'OTU', taxa = NULL, map = NULL, lineage = FALSE, md = FALSE, safe = FALSE) {
  
  df <- taxa_rollup(
    biom    = biom, 
    rank    = rank, 
    taxa    = taxa, 
    map     = map, 
    lineage = lineage, 
    md      = md, 
    safe    = safe,
    long    = TRUE )
  
  return (df)
}



#' Generate a matrix of samples by taxa, at the specified taxonomic rank.
#' 
#' @name taxa_matrix
#'
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'        object, as returned from \link{read_biom}. For matrices, the rows and 
#'        columns are assumed to be the taxa and samples, respectively.
#'        
#' @param rank  The taxonomic rank. E.g. \code{rank = "OTU"}, 
#'        \code{rank = "Phylum"}, etc. May also be given numerically: 0 for
#'        OTU, 1 for the highest level (i.e. Kingdom), and extending to the number
#'        of taxonomic ranks encoded in the original biom file. Run 
#'        \code{taxa_ranks(biom)} to fetch the names of all available ranks.
#'        
#' @param taxa   Limit the number of taxa returned in the matrix. Depending on 
#'        the value type, a different filter is applied.
#'        \describe{
#'          \item{\code{NULL}}{ Retain all taxa. (The default.) }
#'          \item{Integer >= 1}{ The N most abundant taxa are retained. Example:
#'            \code{taxa_matrix(biom, 'Genus', 5)} returned the five most
#'            abundant genera. }
#'          \item{Numeric between 0 and 1}{ Taxa are retained if their abundance 
#'            is greater than or equal to N. Example: 
#'            \code{taxa_matrix(biom, 'Phylum', 0.1)} returns all phyla with a
#'            relative abundance of at least 10 \%. }
#'          \item{Character vector}{ Only these taxa names are retained. Example: 
#'            \code{taxa_matrix(biom, 'Phylum', c("Firmicutes", "Bacteroidetes"))}
#'            will only retain those two phyla. }
#'        }
#'        
#' @param map  A character matrix defining the value that each taxa IDs is
#'        assigned for each taxonomic rank. If \code{map=NULL} and \code{biom} is a
#'        \code{BIOM} class object, the map will be automatically loaded from 
#'        \code{biom$taxonomy}. \code{map} must not be \code{null} when \code{biom}
#'        is a \code{matrix} or \code{simple_triplet_matrix}. See the example below
#'        for an example of \code{map}'s structure.
#'        
#' @param lineage  Include all ranks in the name of the taxa. For instance,
#'        setting to \code{TRUE} will produce 
#'        \code{Bacteria; Actinobacteria; Coriobacteriia; Coriobacteriales}. 
#'        Whereas setting to \code{FALSE} (the default) will return simply
#'        \code{Coriobacteriales}. You want to set this to TRUE if you have
#'        genus names (such as \emph{Incertae_Sedis}) that map to multiple higher
#'        level ranks.
#'        
#' @param sparse  If true, returns a sparse matrix as described by 
#'        \code{slam::simple_triplet_matrix}, otherwise returns a normal R
#'        matrix object. Sparse matrices will likely be considerably more
#'        memory efficient in this scenario.
#'        
#' @return A numeric matrix with samples as column names, and taxonomic
#'         identifiers as row names.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     taxa_ranks(hmp50)
#'     
#'     phyla <- taxa_matrix(hmp50, 'Phylum')
#'     phyla[1:4,1:6]
#'     
#'     # Custom matrices should be formatted like so:
#'     counts <- counts(hmp50)
#'     map    <- taxonomy(hmp50)
#'     
#'     counts[1:3,1:6]
#'     map[1:3,1:4]
#'     
#'     phyla <- taxa_matrix(counts, 'Phylum', map=map)
#'     phyla[1:3,1:6]
#'

taxa_matrix <- function (
    biom, rank = 'OTU', taxa = NULL, map = NULL, lineage = FALSE, sparse = FALSE) {
  
  mtx <- taxa_rollup(
    biom    = biom, 
    rank    = rank, 
    taxa    = taxa, 
    map     = map, 
    lineage = lineage, 
    sparse  = sparse )
  
  return (mtx)
}




taxa_rollup <- function (
    biom, rank = 'OTU', taxa = NULL, 
    map = NULL, lineage = FALSE, sparse = FALSE, 
    long = FALSE, md = FALSE, safe = FALSE ) {
  
  #________________________________________________________
  # Get the input into a simple_triplet_matrix
  #________________________________________________________
  
  if (is(biom, "simple_triplet_matrix")) { counts <- biom
  } else if (is(biom, "BIOM"))           { counts <- biom$counts
  } else if (is(biom, "matrix"))         { counts <- slam::as.simple_triplet_matrix(biom)
  } else {
    stop(simpleError("biom must be a matrix, simple_triplet_matrix, or BIOM object."))
  }
  
  if (is_null(map) && is(biom, "BIOM")) map <- biom$taxonomy
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  
  if (!is(biom, "BIOM") && is_null(map))
    stop(simpleError("Please provide a taxonomy map."))
  
  if (!is(map, "matrix") || !identical(typeof(map), "character"))
    stop(simpleError("Taxonomy map must be a character matrix."))
  
  if (!all(rownames(counts) %in% rownames(map)))
    stop(simpleError("Taxonomy map does not include all taxa IDs."))
  
  if (is(rank, "character"))
    rank <- pmatch(tolower(rank), tolower(c('OTU', colnames(map)))) - 1
  
  if (length(rank) != 1 || !is(rank, "numeric") || is.na(rank))
    stop(simpleError("Invalid taxonomic rank specified."))
  
  if (rank > ncol(map))
    stop(simpleError("Invalid taxonomic rank specified."))
  
  
  if (rank == 0) {
    
    #________________________________________________________
    # The OTU table is easy
    #________________________________________________________
    
    res <- biom[['counts']]
    
  } else {
  
    
    #________________________________________________________
    # compute abundance matrix
    #________________________________________________________
    
    if (identical(lineage, TRUE)) rank <- 1:rank
    
    res <- try(silent=TRUE, {
      otuNames <- counts$dimnames[[1]]
      otu2taxa <- map[otuNames, rank, drop=FALSE]
      otu2taxa <- apply(otu2taxa, 1L, paste, collapse="; ")
      slam::rollup(counts, 1L, otu2taxa, sum)
    })
    if (!is(res, "simple_triplet_matrix"))
      stop(simpleError("Unable to group by taxonomic level."))
    
    res <- res[order(tolower(rownames(res))), colnames(counts), drop=FALSE]
  }
  
  
  #________________________________________________________
  # Only retain taxa of interest (by abundance or name)
  #________________________________________________________
  if (!is_null(taxa)) {
    
    if (is.numeric(taxa) && length(taxa) == 1) {
      rel <- sort(slam::row_means(t(t(res) / slam::col_sums(res))), decreasing = TRUE)
      if (taxa >= 1) { taxa <- head(names(rel), taxa) 
      } else         { taxa <- names(rel)[rel >= taxa] }
      
    } else if (is.character(taxa)) {
      taxa <- intersect(as.character(taxa), rownames(res))
      
    } else {
      stop ("Invalid argument for taxa parameter in taxa_rollup().")
    }
    
    if (length(taxa) == 0)
      stop("No taxa match the criteria: ", capture.output(str(taxa)))
    
    res <- res[taxa,,drop=FALSE]
  }
  
  
  #________________________________________________________
  # Depending on params, taxa might be in rows or cols. Set attr.
  #________________________________________________________
  taxa_in <- "rows"
  
  
  #________________________________________________________
  # Return the slam matrix if sparse=TRUE
  #________________________________________________________
  if (isTRUE(sparse)) {
    attr(res, 'taxa_in') <- taxa_in
    return (res)
  }
  
  
  #________________________________________________________
  # Convert to R matrix
  #________________________________________________________
  res <- t(as.matrix(res))
  taxa_in <- "cols"
  
  
  #________________________________________________________
  # Pivot Longer
  #________________________________________________________
  if (isTRUE(long)) {
    res <- data.frame(
      check.names      = FALSE, 
      stringsAsFactors = FALSE,
      '.sample'        = rownames(res)[row(res)],
      '.taxa'          = colnames(res)[col(res)],
      '.value'         = as.numeric(res)
    )
    if (is_null(taxa)) { res[['.taxa']] %<>% factor()
    } else             { res[['.taxa']] %<>% factor(levels = taxa) }
      
    
    taxa_in <- "rows"
    
  } else if (isFALSE(md)) {
    
    # Just a simple numeric matrix. rownames = samples, colnames = taxa
    attr(res, 'taxa_in') <- taxa_in
    return (res)
    
  } else {
    
    # Don't pivot, but convert matrix to data.frame with 'Sample' column
    res <- data.frame(
      check.names      = FALSE, 
      stringsAsFactors = FALSE, 
      '.sample'        = rownames(res), 
      res
    )
  }
  
  
  #________________________________________________________
  # Add Metadata
  #________________________________________________________
  if (!isFALSE(md)) {
    
    df <- metadata(biom)[,unique(md),drop=F]
    
    if (isTRUE(long)) {
      res <- merge(res, df, by.x = ".sample", by.y = "row.names")
    } else {
      res <- merge(df, res, by.x = "row.names", by.y = ".sample")
      names(res)[1] <- ".sample"
      rownames(res) <- res[['.sample']]
    }
    taxa_in <- "rows"
  }
  
  if (isTRUE(long) && length(unique(res[[".taxa"]])) > 1)
    attr(res, 'facet') <- ifelse(isTRUE(safe), ".taxa", "Taxa")
  
  
  #________________________________________________________
  # Convert auto-generated column names to less-safe values
  #________________________________________________________
  if (isFALSE(safe)) {
    f <- c("Sample", "Taxa", "Abundance")
    if (isFALSE(long)) f <- "Sample"
    colnames(res)[seq_along(f)] <- f
    
  } else if (is.na(safe) && isFALSE(long)) {
    res <- res[,-1,drop=F]
  }
  
  attr(res, 'taxa_in') <- taxa_in
  return (res)
}
