#' Display a dendrogram of the phylogenetic tree.
#' 
#' @name tree_plot
#' 
#' @family plotting
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#' 
#' @param layout   A layout option supported by \link{ggtree::ggtree}:
#'        \code{"rectangular"}, \code{"dendrogram"}, \code{"slanted"}, 
#'        \code{"ellipse"}, \code{"roundrect"}, \code{"fan"}, 
#'        \code{"circular"}, \code{"inward_circular"}, \code{"radial"}, 
#'        \code{"equal_angle"}, \code{"daylight"} or \code{"ape"}.
#'        Default: \code{"rectangular"}.
#' 
#' @param tiplab   Label the tree leafs with taxa names. Default: \code{TRUE}.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     tree_plot(hmp50)
#'     
#'
tree_plot <- function (biom, layout = "rectangular", tiplab = FALSE, rank = NULL) {
  
  # Sanity checks
  #________________________________________________________
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  if (!is.null(rank) && !all(rank %in% metrics(biom, 'rank')))
    stop(
      "`rank` argument to tree_plot must be one of: ", 
      paste(collapse = ", ", metrics(biom, 'rank')) )
  
  
  # Package and print this call
  #________________________________________________________
  params  <- as.list(environment())
  history <- sprintf("tree_plot(%s)", as.args(params, fun = tree_plot))
  
  
  
  # Use ggtree to render a simple tree
  #________________________________________________________
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  
  tree <- as.cmd(phylogeny(biom))
  
  gglayers <- list()
  gglayers %<>% ggpush(ggtree(tr = tree, layout = layout))
  
  if (tiplab) {
    gglayers %<>% ggpush(geom_tiplab())
    gglayers %<>% ggpush(vexpand(.02, 1))
    gglayers %<>% ggpush(vexpand(.02, -1))
  }
  
  
  # Label clades at the specified rank.
  # Annotates the largest clade with 100% same taxon.
  #________________________________________________________
  if (!is.null(rank)) {
    
    map    <- taxonomy(biom, rank)
    result <- sapply(unique(map), function (i) {
      list(node = NA, size = 0) }, simplify = FALSE)
    
    n <- length(tree$tip.label)
    
    clade_leafs <- function (i) {
      
      if (i <= n) {
        otu   <- tree$tip.label[i]
        taxon <- unname(map[otu])
        if (is.na(result[[taxon]][['node']])) {
          result[[taxon]][['node']] <<- i
          result[[taxon]][['size']] <<- 1
        }
        return (otu)
      }
      
      children <- tree$edge[which(tree$edge[,1] == i),2]
      otus     <- unlist(lapply(children, clade_leafs))
      
      if (any(is.na(otus)))               return (NA)
      if (length(unique(map[otus])) != 1) return (NA)
      
      taxon <- unname(map[otus[[1]]])
      if (result[[taxon]][['size']] < length(otus)) {
        result[[taxon]][['node']] <<- i
        result[[taxon]][['size']] <<- length(otus)
      }
      return (otus)
      
    }
    invisible(clade_leafs(n + 1))
    
    for (taxon in names(result)) {
      node <- result[[taxon]][['node']]
      if (is.na(node)) next
      gglayers %<>% ggpush(geom_cladelab(node = node, label = taxon))
    }
    
    gglayers %<>% ggpush(hexpand(.2, 1))
  }
  
  
  p <- ggbuild(gglayers)
  p$plot_env <- emptyenv()
  
  
  # Attach history of biom modifications and this call
  #________________________________________________________
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  
  return (p)
}
