

#' Generate a matrix of samples by taxa, at the specified taxonomic rank.
#' 
#' @name taxa_rollup
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
#' @param long  Format of returned data. The default, \code{FALSE}, returns a
#'        matrix of samples by taxa. If set to \code{TRUE}, returns a data
#'        frame with column names 'Samples', 'Taxa', and 'Abundance'.
#'        
#' @param md  Include metadata in the output data frame? Options are: 
#'        \describe{
#'          \item{\code{FALSE}}{ Don't include metadata. (Default) }
#'          \item{\code{TRUE}}{ Include all metadata. }
#'          \item{\emph{character vector}}{ Include only the specified metadata columns. }
#'        }
#'        
#' @param unc  How to handle unclassified, uncultured, and similarly ambiguous
#'        taxa names. The default, \code{"singly"}/\code{"s"}, replaces them 
#'        with the OTU name. Other options are \code{"grouped"}/\code{"g"}, 
#'        which replaces them with a higher rank's name, \code{"drop"}/\code{"d"} 
#'        which excludes them from the result, and \code{"asis"}/\code{"a"} to
#'        not check/modify any taxa names.
#'        
#' @param other  Sum all non-itemized taxa into an "Other" taxa. The default,
#'        \code{FALSE}, only returns taxa matched by the \code{`taxa`} 
#'        argument. Specifying \code{TRUE} adds "Other" to the returned set.
#'        A string can also be given to imply \code{TRUE}, but with that
#'        value as the name to use in place of "Other".
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

taxa_rollup <- function (
    biom, rank = 'OTU', taxa = NULL, map = NULL, 
    lineage = FALSE, sparse = FALSE, long = FALSE, 
    md = FALSE, unc = "singly", other = FALSE ) {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("taxa_rollup", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))

  
  #________________________________________________________
  # Get the input into a simple_triplet_matrix
  #________________________________________________________
  if (is(biom, "simple_triplet_matrix")) { counts <- biom
  } else if (is(biom, "BIOM"))           { counts <- biom$counts
  } else if (is(biom, "matrix"))         { counts <- slam::as.simple_triplet_matrix(biom)
  } else {
    stop(simpleError("biom must be a matrix, simple_triplet_matrix, or BIOM object."))
  }
  
  
  #________________________________________________________
  # Ensure we have a valid taxonomy map.
  #________________________________________________________
  if (is.matrix(map)) {
    stopifnot(is(map, "matrix"))
    stopifnot(identical(typeof(map), "character"))
    stopifnot(!is.null(rownames(biom)))
    stopifnot(!is.null(colnames(biom)))
    stopifnot(all(rownames(counts) %in% rownames(map)))
    
  } else if (is(biom, "BIOM")) {
    map <- biom[['taxonomy']]
    map <- cbind(map, 'OTU' = rownames(map))
    
  } else {
    stop("Invalid taxonomy map.")
  }
  
  
  
  #________________________________________________________
  # Sanity check arguments for rank and unc.
  #________________________________________________________
  if (is_null(rank))      rank <- ncol(map)
  if (is_character(rank)) rank <- pmatch(tolower(rank), tolower(colnames(map)))
  stopifnot(is_integerish(rank) && all(rank > 0) && all(rank <= ncol(map)))
  
  unc <- match.arg(tolower(unc), c("singly", "grouped", "drop", "asis"))
  
  
  
  #________________________________________________________
  # Handle ambiguous taxa names.
  #________________________________________________________
  map <- tryCatch(
    
    expr = local({
    
      if (isTRUE(lineage)) rank <- range(c(1,rank))
      
      map <- taxonomy(biom = map, ranks = rank, unc = unc)
      map <- apply(map, 1L, paste, collapse="; ")
      
      return (map)
    }), 
    error = function (e) {
      stop("Can't resolve taxa names from map and unc. ", e)
    })
  
  
  
  
  mtx <- local({
    
    
    #________________________________________________________
    # The OTU table is easy
    #________________________________________________________
    if (rank == 0) {
      mtx <- biom[['counts']]
      
      return (mtx)
    }
    
    
    
    #________________________________________________________
    # Compute abundance matrix for this rank.
    #________________________________________________________
    mtx <- tryCatch(
      
      expr = local({
        
        mtx <- slam::rollup(counts[names(map),], 1L, map, sum)
        mtx <- mtx[order(tolower(rownames(mtx))), colnames(counts), drop=FALSE]
        
        stopifnot(is(mtx, "simple_triplet_matrix"))
        
        return (mtx)
      }),
      
      error = function (e) {
        stop("Unable to group by taxonomic level: ", e)
      })
    
    
    return (mtx)
  })
  
  
  #________________________________________________________
  # Only retain taxa of interest (by abundance or name).
  #________________________________________________________
  if (!is_null(taxa)) {
    
    if (is.numeric(taxa) && length(taxa) == 1) {
      rel <- sort(slam::row_means(t(t(mtx) / slam::col_sums(mtx))), decreasing = TRUE)
      if (taxa >= 1) { taxa <- head(names(rel), taxa) 
      } else         { taxa <- names(rel)[rel >= taxa] }
      
    } else if (is.character(taxa)) {
      taxa <- intersect(as.character(taxa), rownames(mtx))
      
    } else {
      stop ("Invalid argument for taxa parameter in taxa_rollup().")
    }
    
    if (length(taxa) == 0)
      stop("No taxa match the criteria: ", capture.output(str(taxa)))
    
    mtx <- mtx[taxa,,drop=FALSE]
  }
  
  
  #________________________________________________________
  # Add an 'Other' taxa.
  #________________________________________________________
  if (isTRUE(other)) other <- "Other"
  if (is_scalar_character(other)) {
    
  }
  
  
  #________________________________________________________
  # Depending on params, taxa might be in rows or cols.
  #________________________________________________________
  taxa_in <- "rows"
  
  
  #________________________________________________________
  # Return the slam matrix if sparse=TRUE
  #________________________________________________________
  if (isTRUE(sparse)) {
    attr(mtx, 'taxa_in') <- taxa_in
    return (mtx)
  }
  
  
  #________________________________________________________
  # Convert to R matrix
  #________________________________________________________
  mtx <- t(as.matrix(mtx))
  taxa_in <- "cols"
  
  
  #________________________________________________________
  # Pivot Longer
  #________________________________________________________
  if (isTRUE(long)) {
    df <- data.frame(
      check.names      = FALSE, 
      stringsAsFactors = FALSE,
      '.sample'        = rownames(mtx)[row(mtx)],
      '.taxa'          = colnames(mtx)[col(mtx)],
      '.abundance'     = as.numeric(mtx)
    )
    if (is_null(taxa)) { df[['.taxa']] %<>% factor()
    } else             { df[['.taxa']] %<>% factor(levels = taxa) }
    
    
    taxa_in <- "rows"
    
  } else if (isFALSE(md)) {
    
    # Just a simple numeric matrix. rownames = samples, colnames = taxa
    attr(mtx, 'taxa_in') <- taxa_in
    return (mtx)
    
  } else {
    
    # Don't pivot, but convert matrix to data.frame with 'Sample' column
    df <- data.frame(
      check.names      = FALSE, 
      stringsAsFactors = FALSE, 
      '.sample'        = rownames(mtx), 
      mtx
    )
  }
  
  
  #________________________________________________________
  # Add Metadata
  #________________________________________________________
  if (!isFALSE(md)) {
    
    md <- metadata(biom)[,unique(md),drop=F]
    
    if (isTRUE(long)) {
      df <- merge(df, md, by.x = ".sample", by.y = "row.names")
    } else {
      df <- merge(md, df, by.x = "row.names", by.y = ".sample")
      names(df)[1] <- ".sample"
      rownames(df) <- df[['.sample']]
    }
    taxa_in <- "rows"
  }
  
  if (isTRUE(long) && length(unique(df[[".taxa"]])) > 1)
    attr(df, 'facet') <- ".taxa"
  
  
  attr(df, 'response') <- ".abundance"
  attr(df, 'taxa_in')  <- taxa_in
  
  set_cache_value(cache_file, df)
  return (df)
}








#' Generate a data.frame of taxa abundance for each sample at the specified taxonomic rank.
#' 
#' @name taxa_table
#' 
#' @inherit taxa_rollup params
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
    biom, rank = 'OTU', taxa = NULL, map = NULL, 
    lineage = FALSE, md = FALSE, unc = "singly", 
    other = FALSE ) {
  
  df <- taxa_rollup(
    biom    = biom, 
    rank    = rank, 
    taxa    = taxa, 
    map     = map, 
    lineage = lineage, 
    md      = md, 
    unc     = unc, 
    other   = other, 
    long    = TRUE )
  
  return (df)
}



#' Generate a matrix of samples by taxa, at the specified taxonomic rank.
#' 
#' @name taxa_matrix
#' 
#' @inherit taxa_rollup params
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
    biom, rank = 'OTU', taxa = NULL, map = NULL, 
    lineage = FALSE, sparse = FALSE, unc = "singly", 
    other = FALSE ) {
  
  mtx <- taxa_rollup(
    biom    = biom, 
    rank    = rank, 
    taxa    = taxa, 
    map     = map, 
    lineage = lineage, 
    sparse  = sparse, 
    unc     = unc, 
    other   = other )
  
  return (mtx)
}
