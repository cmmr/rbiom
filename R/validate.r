
#===============================================
# Return a list of recognized options
#===============================================
ord_metrics   <- function (biom=NULL) c("PCoA", "tSNE", "NMDS")
adiv_metrics  <- function (biom=NULL) c("Depth", "OTUs", "Shannon", "Chao1", "Simpson", "InvSimpson")
bdiv_metrics  <- function (biom=NULL) c("Manhattan", "Euclidean", "Bray-Curtis", "Jaccard", "UniFrac")
dist_metrics  <- function (biom=NULL) c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
rank_metrics  <- function (biom=NULL) c(taxa.ranks(biom), 'OTU')
taxon_metrics <- function (biom=NULL) unique(c(as.character(taxonomy(biom)), taxa.names(biom)))
meta_metrics  <- function (biom=NULL) colnames(metadata(biom))
clust_metrics <- function (biom=NULL) c("heatmap", "average", "UPGMA", "ward", "mcquitty", "WPGMA", "single", "median", "WPGMC", "complete", "centroid", "UPGMC")
other_metrics <- function (biom=NULL) c("Rarefied", "Reads", "Samples", ".", "stacked") %>% structure(., mode=.)
all_metrics   <- function (biom=NULL) {
  v <- unlist(sapply(
    USE.NAMES = FALSE, 
    X         = c("ord", "adiv", "bdiv", "rank", "taxon", "meta", "clust", "other"), 
    FUN       = function (i) {
      k <- do.call(paste0(i, "_metrics"), list(biom=biom))
      n <- attr(k, 'mode', exact = TRUE)
      if (is.null(n)) n <- rep_len(i, length(k))
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
  
  if (is.null(attr(opts, 'mode', exact = TRUE))) {
    attr(vals, 'mode') <- rep_len(mode, length(vals))
  } else {
    attr(vals, 'mode') <- attr(opts, 'mode', exact = TRUE)[okay]
  }
  
  
  #-----------------------------------------------
  # Sanity checks
  #-----------------------------------------------
  if (!multi && length(vals) > 1) 
    stop("Only a single metric is allowed. Found: ", paste0(collaspe=", ", vals))
  
  uModes <- unique(attr(vals, 'mode', exact = TRUE))
  if (!mixed && length(uModes) > 1)
    stop("All metrics must be the same type. Found: ", paste0(collaspe=", ", uModes))
  
  
  #-----------------------------------------------
  # Solo metadata columns get further attributes
  #-----------------------------------------------
  if (identical(attr(vals, 'mode', exact = TRUE), "meta")) {
    
    #-----------------------------------------------
    # Look for '==' or '!=' prefixes
    #-----------------------------------------------
    attr(vals, 'op') <- attr(metrics, 'op', exact = TRUE)
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

