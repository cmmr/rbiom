#' Get the taxonomy table.
#' 
#' @param biom  A \code{BIOM} object, as returned from \link{read_biom}.
#' 
#' @param ranks  The taxonomic ranks to return in the matrix, or \code{NULL} 
#'        for all of them. Default: \code{NULL}
#'        
#' @param unc  How to handle unclassified, uncultured, and similarly ambiguous
#'        taxa names. Default: \code{"asis"}
#'        \itemize{
#'          \item{\code{"asis"} - }{ Don't check/modify any taxa names. }
#'          \item{\code{"singly"} - }{ Replace with "Unc. <OTU ID>". }
#'          \item{\code{"grouped"} - }{ Replace with "Unc. <Higher Rank>". }
#'          \item{\code{"drop"} - }{ Don't include in the returned matrix. }
#'        }
#' 
#' @return A character matrix with taxa/OTU IDs as row names and taxa ranks as
#'         column names. An 'OTU' column is always added as the last column and
#'         matches the row names.
#'         
#' @family accessors
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     taxonomy(hmp50)[1:4,]
#'     taxonomy(hmp50, c("Family", "Genus"))[1:4,]
#'     head(taxonomy(hmp50, "Genus")[,])
#'     
#'     # Sometimes taxonomic names are incomplete
#'     taxonomy(hmp50)[c(53,107,139), 2:6]
#'     
#'     # rbiom can insert more descriptive placeholders
#'     taxonomy(hmp50, unc = "singly")[c(53,107,139), 3:6]
#'     taxonomy(hmp50, unc = "grouped")[c(53,107,139), 4:6]
#'

taxonomy <- function (biom, ranks = NULL, unc = "asis") {
  
  #________________________________________________________
  # Check that `biom` look right.
  #________________________________________________________
  map <- tryCatch(
    expr = local({
  
      stopifnot(is(biom, "BIOM") || is(biom, "matrix"))
      if (is(biom, 'BIOM')) biom <- biom[['taxonomy']]
      stopifnot(identical(typeof(biom), "character"))
      stopifnot(!is.null(rownames(biom)))
      stopifnot(!is.null(colnames(biom)))
    
      return (biom)
    }), 
    error = function (e) {
      stop(
        "Invalid argument for `biom` in rbiom::taxonomy(). ",
        "Must be a BIOM object or a character matrix with ",
        "OTU ids as row names and taxa ranks as column names.",
        "\n\nError: ", e, "\n\n")
    })
  
  
  #________________________________________________________
  # Make sure the last column is 'OTU'.
  #________________________________________________________
  if ("OTU" %in% colnames(map))
    map <- map[,which(colnames(map) != "OTU"),drop=FALSE]
  map <- cbind(map, 'OTU' = rownames(map))
  
  
  
  #________________________________________________________
  # Sanity check `ranks` and `unc` arguments.
  #________________________________________________________
  if (is_null(ranks))      ranks <- seq_len(ncol(map))
  if (is_character(ranks)) ranks <- pmatch(tolower(ranks), tolower(colnames(map)))
  stopifnot(!any(is.na(ranks)))
  stopifnot(is_integerish(ranks) && all(ranks > 0) && all(ranks <= ncol(map)))
  
  unc <- match.arg(tolower(unc), c("asis", "singly", "grouped", "drop"))
  
  
  
  #________________________________________________________
  # Transform the taxa names.
  #________________________________________________________
  if (unc != "asis")
    map <- tryCatch(
      expr = local({
      
        # Discard technical prefixes/suffixes.
        #________________________________________________________
        map <- sub("^.__", "", map) # Remove leading p__ c__ etc
        map <- sub("^_+",  "", map) # Remove leading underscores
        map <- sub(";$",   "", map) # Remove trailing semicolons
        
        
        # "g" => NA; "Unknown Order" => NA
        #________________________________________________________
        regex <- ".*(unknown|uncultured|unclassified|unidentified|incertae.sedis).*"
        map[which(nchar(map) < 2)] <- NA
        map[grep(regex, map, ignore.case = TRUE)] <- NA
        
        
        # "R_7_group" => "R_7"
        #________________________________________________________
        map[] <- sub("_group$", "", map, ignore.case = TRUE)
        
        
        # "Family XIII" => "Clostridiales XIII"
        #________________________________________________________
        map <- t(apply(map, 1L, function (x) {
          regex <- "^family\\s"
          if (!is.null(prefixed <- grep(regex, x, ignore.case = TRUE)))
            for (i in prefixed)
              if (!is.null(ideal <- grep("^[A-Z][a-z]+$", x[seq_len(i - 1)], value = TRUE)))
                x[i] <- sub(regex, paste0(tail(ideal, 1), " "), x[i], ignore.case = TRUE)
          
          return (x)
        }))
        
        
        # Replace NA with "Unc. <OTU ID>".
        #________________________________________________________
        if (identical(unc, "singly")) {
          x <- which(is.na(map))
          map[x] <- paste("Unc.", rownames(map)[row(map)[x]])
        }
        
        
        # Replace NA with "Unc. <Higher Rank>".
        #________________________________________________________
        if (identical(unc, "grouped"))
          for (i in which(!complete.cases(map)))
            for (j in rev(which(is.na(map[i,]))))
              if (!is.null(x <- na.omit(c("N/A", map[i,seq_len(j - 1)]))))
                map[i,j] <- paste("Unc. ", tail(x, 1))
        
        
        # Drop any row with an NA in one of the ranks column.
        #________________________________________________________
        if (identical(unc, "drop"))
          map <- map[complete.cases(map[,ranks,drop=FALSE]),,drop=FALSE]
        
        
        return (map)
      }), 
      
      error = function (e) {
        stop("Error in renaming taxa with taxonomy() call: ", e)
      })
  
  
  map <- map[,ranks,drop=FALSE]
  
  return (map)
}




#' Set the taxonomy table.
#' 
#' @param x  A \code{BIOM} object, as returned from \link{read_biom}.
#' 
#' @param value  A character matrix with rownames \code{taxa_names(x)}. If
#'        there are more rownames than taxa names, the matrix will be subset.
#'         
#' @family setters
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     taxonomy(hmp50)[1:4,]
#'     
#'     # Sometimes taxonomic names are incomplete
#'     taxonomy(hmp50)[c(107,170,194), 'Genus', drop=FALSE]
#'     
#'     # rbiom can insert more descriptive placeholders
#'     taxonomy(hmp50, unc = "grouped")[c(107,170,194), 'Genus', drop=FALSE]
#'

`taxonomy<-` <- function (x, value) {
  
  stopifnot(is(x, 'BIOM'))
  stopifnot(is.matrix(value))
  stopifnot(typeof(value) == "character")
  stopifnot(all(taxa_names(x) %in% rownames(value)))
  
  x[['taxonomy']] <- value[taxa_names(x),]
  
  return (x) 
}


