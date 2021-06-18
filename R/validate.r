
#===============================================
# Return a list of recognized options
#===============================================
ord_metrics   <- function (biom=NULL) c("PCoA", "tSNE", "NMDS")
adiv_metrics  <- function (biom=NULL) c("Depth", "OTUs", "Shannon", "Chao1", "Simpson", "InvSimpson")
bdiv_metrics  <- function (biom=NULL) c("Manhattan", "Euclidean", "Bray-Curtis", "Jaccard", "UniFrac")
rank_metrics  <- function (biom=NULL) c(taxa.ranks(biom), 'OTU')
taxon_metrics <- function (biom=NULL) unique(c(as.character(taxonomy(biom)), taxa.names(biom)))
meta_metrics  <- function (biom=NULL) colnames(metadata(biom))
other_metrics <- function (biom=NULL) c("Rarefied", "Reads", "Samples", "PctReads", ".") %>% structure(., mode=.)
all_metrics   <- function (biom=NULL) {
  v <- unlist(sapply(
    USE.NAMES = FALSE, 
    X         = c("ord", "adiv", "bdiv", "rank", "taxon", "meta", "other"), 
    FUN       = function (i) {
      k <- do.call(paste0(i, "_metrics"), list(biom=biom))
      n <- if (is.null(attr(k, 'mode'))) rep_len(i, length(k)) else attr(k, 'mode')
      setNames(n, as.vector(k))
    }))
  structure(names(v), mode=unname(v))
}


#===============================================
# Cleanup dynamic options that we can recognize
#===============================================
validate_metrics <- function (biom, metrics, mode="all", multi=FALSE, mixed=FALSE) {
  
  opts <- do.call(paste0(mode, "_metrics"), list(biom=biom))
  okay <- pmatch(tolower(sub("^[!=]=", "", metrics)), tolower(opts))
  vals <- opts[okay]
  
  missing <- which(is.na(okay))
  if (length(missing) > 0)
    stop("Invalid or ambiguous metric(s): ", paste(collapse = ", ", metrics[missing]))
  
  if (is.null(attr(opts, 'mode'))) { attr(vals, 'mode') <- rep_len(mode, length(vals))
  } else                           { attr(vals, 'mode') <- attr(opts, 'mode')[okay] }
  
  
  #-----------------------------------------------
  # Sanity checks
  #-----------------------------------------------
  if (!multi && length(vals) > 1) 
    stop("Only a single metric is allowed. Found: ", paste0(collaspe=", ", vals))
  
  uModes <- unique(attr(vals, 'mode'))
  if (!mixed && length(uModes) > 1)
    stop("All metrics must be the same type. Found: ", paste0(collaspe=", ", uModes))
  
  
  #-----------------------------------------------
  # Solo metadata columns get further attributes
  #-----------------------------------------------
  if (identical(attr(vals, 'mode'), "meta")) {
    
    #-----------------------------------------------
    # Look for '==' or '!=' prefixes
    #-----------------------------------------------
    attr(vals, 'op') <- attr(metrics, 'op')
    if (substr(metrics, 1, 2) %in% c("==", "!="))
      attr(vals, 'op') <- substr(metrics, 1, 2)
    
    #-----------------------------------------------
    # Further classify as 'factor' or 'numeric'
    #-----------------------------------------------
    cl <- class(metadata(biom, vals))
    if (any(cl %in% c('factor', 'character', 'logical'))) {
      attr(vals, 'mode') <- "factor"
      
    } else if (any(cl %in% c('numeric', 'date', 'integer', 'complex', 'Date'))) {
      attr(vals, 'mode') <- "numeric"
      
    } else {
      attr(vals, 'mode') <- head(cl, 1)
    }
  }
  
  return (vals)
}

