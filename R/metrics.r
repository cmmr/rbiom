#' List all the options for each type of metric.
#' 
#' @name metrics
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#' 
#' @param mode   One of the following options:
#' \itemize{
#'   \item{\bold{ord} - }{ Ordination }
#'   \item{\bold{adiv} - }{ Alpha Diversity }
#'   \item{\bold{dist} - }{ Distance (Beta Diversity) }
#'   \item{\bold{clust} - }{ Clustering }
#'   \item{\bold{rank} - }{ Taxonomic Rank }
#'   \item{\bold{taxon} - }{ Taxa Names }
#'   \item{\bold{meta} - }{ Metadata Fields }
#'   \item{\bold{all} - }{ All of the Above }
#' }
#'        
#' @return A character vector of supported values. 
#'         For \code{mode = "all"}, a named \code{list()} of character vectors.
#' @export
#' @seealso \code{\link{adiv_table}} \code{\link{bdiv_table}} \code{\link{ordinate}}
#' @examples
#'     library(rbiom)
#'     
#'     metrics(hmp50, 'adiv')
#'     metrics(hmp50, 'dist')
#'
metrics <- function (biom, mode = "all", tree=NULL) {
  
  mode  <- tolower(mode)[[1]]
  
  if (!is(biom, 'BIOM') && mode %in% c('rank', 'taxon', 'meta', 'all'))
    stop("Please provide a BIOM object when using mode='", mode, "'")
  
  modes <- c('ord', 'adiv', 'dist', 'clust', 'rank', 'taxon', 'meta', 'bdiv')
  if (mode == 'all')
    return (sapply(X = modes, FUN = rbiom::metrics, biom=biom, tree=tree))
  
  
  
  if        (mode == 'ord')   { c("PCoA", "tSNE", "NMDS", "UMAP") 
  } else if (mode == 'adiv')  { c("OTUs", "Shannon", "Chao1", "Simpson", "InvSimpson") 
  } else if (mode == 'dist')  { c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski") 
  } else if (mode == 'clust') { c("average", "ward", "mcquitty", "single", "median", "complete", "centroid") 
  } else if (mode == 'rank')  { unique(c(taxa_ranks(biom), 'OTU'))
  } else if (mode == 'taxon') { unique(c(as.character(taxonomy(biom)), taxa_names(biom)))
  } else if (mode == 'meta')  { colnames(metadata(biom))
  } else if (mode == 'bdiv')  {
    hasTree <- ifelse(is(biom, 'BIOM'), has_phylogeny(biom), FALSE) || is(tree, 'phylo')
    if (hasTree) { c("Manhattan", "Euclidean", "Bray-Curtis", "Jaccard", "UniFrac")
    } else       { c("Manhattan", "Euclidean", "Bray-Curtis", "Jaccard") }
  } else { NULL }
}


#===============================================
# Cleanup dynamic options that we can recognize
#===============================================
validate_metrics <- function (biom, metrics, mode=NULL, multi=FALSE, mixed=FALSE, ...) {
  
  #-----------------------------------------------
  # Return a list of recognized options
  #-----------------------------------------------
  ord_metrics   <- function (biom, ...) metrics(biom, 'ord',   ...)
  adiv_metrics  <- function (biom, ...) metrics(biom, 'adiv',  ...) %>% c("Depth")
  bdiv_metrics  <- function (biom, ...) metrics(biom, 'bdiv',  ...)
  dist_metrics  <- function (biom, ...) metrics(biom, 'dist',  ...)
  rank_metrics  <- function (biom, ...) metrics(biom, 'rank',  ...) %>% c("Rank")
  taxon_metrics <- function (biom, ...) metrics(biom, 'taxon', ...)
  meta_metrics  <- function (biom, ...) metrics(biom, 'meta',  ...)
  clust_metrics <- function (biom, ...) metrics(biom, 'clust', ...) %>% c("heatmap", "UPGMA", "WPGMA", "WPGMC", "UPGMC")
  other_metrics <- function (biom, ...) c("Rarefied", "Reads", "Samples", ".", "stacked") %>% structure(., mode=.)
  all_metrics   <- function (biom, ...) {
    v <- unlist(sapply(
      USE.NAMES = FALSE, 
      X         = c("ord", "bdiv", "rank", "adiv", "taxon", "meta", "clust", "other"), 
      FUN       = function (i) {
        k <- do.call(paste0(i, "_metrics"), list(biom=biom))
        n <- attr(k, 'mode', exact = TRUE)
        if (is.null(n)) n <- rep_len(i, length(k))
        setNames(n, as.vector(k))
      }))
    structure(names(v), mode=unname(v))
  }
  
  # Have we already validated this value?
  if (length(unique(attr(metrics, 'mode', exact = TRUE))) == 1)
    mode %||=% attr(metrics, 'mode', exact = TRUE)[[1]]
  mode %||=% "all"
  
  
  opts <- do.call(paste0(mode, "_metrics"), c(list(biom=biom), list(...)))
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


