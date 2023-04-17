#' Display a dendrogram of the phylogenetic tree.
#' 
#' @name tree_plot
#' 
#' @family plotting
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#' 
#' @param layout   Any layout option supported by \link{ggtree::ggtree}:
#'        \code{"rectangular"}, \code{"dendrogram"}, \code{"slanted"}, 
#'        \code{"ellipse"}, \code{"roundrect"}, \code{"fan"}, 
#'        \code{"circular"}, \code{"inward_circular"}, \code{"radial"}, 
#'        \code{"equal_angle"}, \code{"daylight"} or \code{"ape"}.
#'        Default: \code{"rectangular"}.
#' 
#' @param tiplab   Label the tree leafs with taxa names. Options are
#'        \code{NULL} (no tip labels) or a taxonomic rank (for example
#'        \code{tiplab = "Genus"}). Default: \code{NULL}.
#' 
#' @param color.by   How to color the tree. Currently only supports the options
#'        \code{NULL} (no coloring) or \code{'.reads'} (color by number of taxa 
#'        observations). Default: \code{NULL}.
#' 
#' @param label,cladelab   Label monophyletic clades. You can specify the
#'        same or different taxonomic ranks for internal (label) and external 
#'        (cladelab) annotations. Default: \code{NULL}.
#' 
#' @param top,right,bottom,left   Add additional space around the tree.
#'        Sometimes necessary for wide text annotations. Set as fraction of
#'        tree's width/height. For instance, \code{right = 1} reserves the 
#'        right half of the plotting area for non-tree elements.
#'        Default: \code{NULL}.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     tree_plot(hmp50)
#'     
#'     hmp1 <- select(hmp50, samples = "HMP10")
#'     tree_plot(hmp1, cladelab = "Phylum", layout = "roundrect")
#'     tree_plot(hmp1, tiplab = "Genus", layout = "fan", color.by = ".reads")
#'
tree_plot <- function (
    biom, layout = "rectangular", 
    tiplab = NULL, color.by = NULL, label = NULL, cladelab = NULL, 
    top = NULL, right = NULL, bottom = NULL, left = NULL, ...) {
  
  
  # Sanity checks
  #________________________________________________________
  local({
    if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
    
    if (!has_phylogeny(biom)) stop("No phylogenetic data present.")
    
    if(!is_null(c(tiplab, label, cladelab))) {
      ranks <- metrics(biom, 'rank')
      for (i in c('tiplab', 'label', 'cladelab')) {
        rank <- get(i)
        if (is_null(rank)) next
        if (!is.character(rank)) stop(i, " must be NULL or character")
        if (length(rank) > 1)    stop(i, " must be NULL or length 1")
        if (!rank %in% ranks)
          stop(i, " must be one of NULL, '", paste(collapse = "', '", ranks), "'")
      }
    }
  })
  
  
  # Package and print this call
  #________________________________________________________
  params          <- c(as.list(environment()), list(...))
  params[['...']] <- NULL
  history <- sprintf("tree_plot(%s)", as.args(params, fun = tree_plot))
  
  
  
  # Assemble the plot layers with initLayer / setLayer.
  #________________________________________________________
  layers <- list()
  attr(layers, 'params')   <- params
  attr(layers, 'function') <- tree_plot
  
  
  
  # Use ggtree to render a phylogenetic tree
  #________________________________________________________
  tr <- tree_data(biom)
  attr(tr, 'display') <- "tree_data(biom)"
  setLayer("ggtree", layout = layout, tr = tr)
  
  
  # Coloring by reads requires an aes mapping.
  #________________________________________________________
  if (identical(color.by, '.reads')) {
    setLayer("ggtree", mapping = as.cmd(aes(color=log(reads))))
    setLayer("scale_color_continuous", low='darkgreen', high='red')
  }
  
  
  # Add OTU names at the tree tips.
  #________________________________________________________
  if (!is_null(tiplab)) {
    initLayer("tiplab")
    
    if (!identical(tiplab, 'OTU')) {
      
      df  <- attr(layers[['ggtree']][['tr']], 'data',  exact = TRUE)
      phy <- attr(layers[['ggtree']][['tr']], 'phylo', exact = TRUE)
      df[['tiplab']] <- taxonomy(biom, tiplab)[phy$tip.label,1][1:nrow(df)]
      attr(layers[['ggtree']][['tr']], 'data') <- df
      remove("df", "phy")
      
      setLayer(
        layer   = "tiplab", 
        data    = as.cmd(function (x) { subset(x, !is.na(tiplab)) }), 
        mapping = aes(label=tiplab) )
    }
  }
  
  
  # Label clades internally (inside) of the tree.
  #________________________________________________________
  if (!is_null(label)) {
    
    mapping <- aes(x=branch, label=!!as.name(label))
    if (identical(color.by, '.reads'))
      mapping <- aes(x=branch, label=!!as.name(label), color=log(reads))
    
    setLayer(
      layer   = "label",
      data    = as.cmd(function (x) { subset(x, !is.na(y)) }, list(y = as.name(label))), 
      mapping = mapping, 
      fill    = 'white' )
  }
  
  
  # Label clades externally (at right) of the tree.
  #________________________________________________________
  if (!is_null(cladelab)) {
    
    mapping <- aes(node=node, label=!!as.name(cladelab))
    if (identical(color.by, '.reads'))
      mapping <- aes(node=node, label=!!as.name(cladelab), color=log(reads))
    
    setLayer(
      layer       = "cladelab",
      data        = as.cmd(function (x) { subset(x, !is.na(y)) }, list(y = as.name(cladelab))), 
      mapping     = mapping,
      offset      = 0.02,
      offset.text = 0.02 )
  }
  
  
  # Add horizontal/vertical padding for text to bleed into.
  #________________________________________________________
  if (!is_null(c(label, cladelab, tiplab))) {
    top    %<>% if.null(0.02)
    bottom %<>% if.null(0.02)
  }
  if (!is_null(cladelab))
    right %<>% if.null(0.2)
  
  if (!is_null(top))    setLayer('top',    .fn = 'vexpand', direction =  1, ratio = top)
  if (!is_null(right))  setLayer('right',  .fn = 'hexpand', direction =  1, ratio = right)
  if (!is_null(bottom)) setLayer('bottom', .fn = 'vexpand', direction = -1, ratio = bottom)
  if (!is_null(left))   setLayer('left',   .fn = 'hexpand', direction = -1, ratio = left)
  
  
  p <- suppressMessages(ggbuild(layers))
  
  
  # Attach history of biom modifications and this call
  #________________________________________________________
  attr(p, 'history') <- c(attr(biom, 'history'), history)
  
  
  
  return (p)
}


#' Provides a 'treedata' S4 object for use in ggtree functions.
#' 
#' @name tree_data
#' 
#' @param biom   A BIOM object, as returned from \link{read_biom}.
#' 
#' @param reads   Include a 'reads' column indicating the sum of taxa 
#'        observations belonging to each node/leaf. Default: \code{TRUE}.
#' 
#' @param clades   Notate top-most monophyletic nodes. The default,
#'        \code{TRUE}, adds a column for every rank in the biom object. A
#'        character vector of ranks can also be passed in. If the vector is
#'        named, then those names are used for naming the columns in the
#'        returned treedata object. Set to \code{NULL} to not return any 
#'        clade notations.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     tree_data(hmp50)
#'     
#'
tree_data <- function (biom, reads = TRUE, clades = TRUE) {
  
  
  # Sanity checks
  #________________________________________________________
  if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
  stopifnot(is.logical(reads) && length(reads) == 1)
  
  tree   <- phylogeny(biom)
  ranks  <- taxa_ranks(biom)
  labels <- tree$tip.label
  nTips  <- length(labels)
  df     <- data.frame(node = seq_len(max(tree$edge)))
  
  
  # Convert `clades` param to named `ranks`.
  #________________________________________________________
  ranks <- local({
    
    stopifnot(is.logical(clades) || is.character(clades))
    
    if (isFALSE(clades)) return (NULL)
    if (isTRUE(clades))  return (taxa_ranks(biom) %>% setNames(., .))
    
    if (length(bad <- setdiff(clades, taxa_ranks(biom))) > 0)
      stop("Invalid taxonomic rank: ", paste(collapse = ", ", bad))
      
    if (is_null(names(clades))) return (setNames(clades, clades))
    
    if (length(i <- which(names(clades) == "")) > 0)
      names(clades)[i] <- as.vector(clades)[i]
    
    return (clades)
  })
  
  
  # We can pre-populate all the data for leaf nodes.
  #________________________________________________________
  if (isTRUE(reads))
    df[1:nTips, 'reads'] <- as.vector(taxa_sums(biom)[labels])
  
  for (i in seq_along(ranks))
    df[1:nTips, names(ranks)[[i]]] <- as.vector(taxonomy(biom, ranks[[i]])[labels,1])
  
  df[1:nTips, 'OTU'] <- labels
  
  
  # Recurse through the tree.
  #________________________________________________________
  traverse <- function (i) {
    
    if (i <= nTips) return (i)
    
    children    <- tree$edge[which(tree$edge[,1] == i),2]
    descendents <- unlist(lapply(children, traverse))
    
    if (isTRUE(reads))
      df[i, 'reads'] <<- sum(df[children, 'reads'])
    
    for (rank in names(ranks)) {
      if (length(unique(df[children, rank])) == 1) {
        df[i, rank]        <<- df[children[[1]], rank]
        df[children, rank] <<- NA
      }
    }
    
    return (c(i, descendents))
  }
  invisible(traverse(nTips + 1))
  
  
  # Convert phylo object to treedata object.
  #________________________________________________________
  tree <- treeio::full_join(tree, df, by = "node")
  
  return (tree)
}
