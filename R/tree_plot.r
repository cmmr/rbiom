#' Display a dendrogram of the phylogenetic tree.
#' 
#' @noRd
#' @keywords internal
#' 
#' @inherit documentation_default
#' 
#' @family phylogeny
#' @family visualization
#' 
#' @param layout   Any layout option supported by [ggtree::ggtree()]:
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
#' @examples
#'     library(rbiom)
#'     
#'     hmp1 <- sample_select(hmp50, samples = "HMP10")
#'     
#'     # Needs the optional ggtree package
#'     if (nzchar(system.file(package = "ggtree")))
#'       tree_plot(hmp50)
#'       
#'     if (nzchar(system.file(package = "ggtree")))
#'       tree_plot(hmp1, cladelab = "Phylum", layout = "roundrect")
#'     
#'     if (nzchar(system.file(package = "ggtree")))
#'       tree_plot(hmp1, tiplab = "Genus", layout = "fan", color.by = ".reads")

tree_plot <- function (
    biom, layout = "rectangular", 
    tiplab = NULL, color.by = NULL, label = NULL, cladelab = NULL, 
    top = NULL, right = NULL, bottom = NULL, left = NULL, ...) {
  
  if (nchar(system.file(package = "ggtree")) == 0)
    stop("Bioconductor R package 'ggtree' must be installed to use tree_plot().")
  
  validate_biom(clone = FALSE)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- eval_envir(environment(), ...)
  history    <- append_history('fig ', params)
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  stopifnot(has_tree(biom))
  validate_rank("tiplab",   null_ok = TRUE)
  validate_rank("label",    null_ok = TRUE)
  validate_rank("cladelab", null_ok = TRUE)
  
  
  
  #________________________________________________________
  # Assemble the plot layers with init_layers / set_layer.
  #________________________________________________________
  layers <- list()
  attr(layers, 'params')   <- params
  attr(layers, 'function') <- tree_plot
  
  
  
  #________________________________________________________
  # Use ggtree to render a phylogenetic tree
  #________________________________________________________
  tr <- tree_data(biom)
  attr(tr, 'display') <- "tree_data(biom)"
  set_layer(params, 'ggtree', layout = layout, tr = tr)
  
  
  #________________________________________________________
  # Coloring by reads requires an aes mapping.
  #________________________________________________________
  if (eq(color.by, '.reads')) {
    set_layer(params, 'ggtree', mapping = as.cmd(aes(color=log(reads))))
    set_layer(params, 'scale_color_continuous', low='darkgreen', high='red')
  }
  
  
  #________________________________________________________
  # Add OTU names at the tree tips.
  #________________________________________________________
  if (!is_null(tiplab)) {
    
    add_layer(params, 'tiplab')
    
    if (!tiplab == '.otu') {
      
      df  <- attr(layers[['ggtree']][['tr']], 'data',  exact = TRUE)
      phy <- attr(layers[['ggtree']][['tr']], 'phylo', exact = TRUE)
      
      if (nrow(df) == 0)
        df <- attr(tr, 'extraInfo', exact = TRUE)[,'node']
      
      df[['tiplab']] <- otu_taxonomy(biom, tiplab)[phy$tip.label][df[['node']]]
      
      attr(layers[['ggtree']][['tr']], 'data') <- df
      remove("df", "phy")
      
      set_layer(
        params = params, 
        layer  = 'tiplab', 
        'data'    = as.cmd(function (x) { subset(x, !is.na(tiplab)) }), 
        'mapping' = aes(label=tiplab) )
    }
  }
  
  
  #________________________________________________________
  # Label clades internally (inside) of the tree.
  #________________________________________________________
  if (!is_null(label)) {
    
    mapping <- aes(x=branch, label=!!as.name(label))
    if (eq(color.by, '.reads'))
      mapping <- aes(x=branch, label=!!as.name(label), color=log(reads))
    
    set_layer(
      params = params, 
      layer  = 'label',
      'data'    = as.cmd(function (x) { subset(x, !is.na(y)) }, list(y = as.name(label))), 
      'mapping' = mapping, 
      'fill'    = 'white' )
  }
  
  
  #________________________________________________________
  # Label clades externally (at right) of the tree.
  #________________________________________________________
  if (!is_null(cladelab)) {
    
    mapping <- aes(node=node, label=!!as.name(cladelab))
    if (eq(color.by, '.reads'))
      mapping <- aes(node=node, label=!!as.name(cladelab), color=log(reads))
    
    set_layer(
      params = params, 
      layer  = 'cladelab',
      'data'        = as.cmd(function (x) { subset(x, !is.na(y)) }, list(y = as.name(cladelab))), 
      'mapping'     = mapping,
      'offset'      = 0.02,
      'offset.text' = 0.02 )
  }
  
  
  #________________________________________________________
  # Add horizontal/vertical padding for text to bleed into.
  #________________________________________________________
  if (!is_null(c(label, cladelab, tiplab))) {
    top    %<>% if.null(0.02)
    bottom %<>% if.null(0.02)
  }
  if (!is_null(cladelab))
    right %<>% if.null(0.2)
  
  if (!is_null(top))    set_layer(params, 'top',    .fn = 'vexpand', direction =  1, ratio = top)
  if (!is_null(right))  set_layer(params, 'right',  .fn = 'hexpand', direction =  1, ratio = right)
  if (!is_null(bottom)) set_layer(params, 'bottom', .fn = 'vexpand', direction = -1, ratio = bottom)
  if (!is_null(left))   set_layer(params, 'left',   .fn = 'hexpand', direction = -1, ratio = left)
  
  
  p <- suppressMessages(ggbuild(layers))
  
  
  attr(p, 'history') <- history
  set_cache_value(cache_file, p)
  
  return (p)
}



#' Provides a 'treedata' S4 object for use in ggtree functions.
#' 
#' @family phylogeny
#' 
#' @param biom   An rbiom object, as returned from [read_biom()].
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
#' @examples
#'     library(rbiom)
#'     
#'     # Needs the optional ggtree package
#'     if (nzchar(system.file(package = "ggtree"))) {
#'     
#'       tree_data(hmp50)
#'       
#'     }

tree_data <- function (biom, reads = TRUE, clades = TRUE) {
  
  if (nchar(system.file(package = "ggtree")) == 0)
    stop("Bioconductor R package 'ggtree' must be installed to use tree_data().")
  
  validate_biom(clone = FALSE)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- eval_envir(environment())
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  # Sanity checks
  #________________________________________________________
  stopifnot(is_scalar_logical(reads) && !is_na(reads))
  
  tree   <- otu_tree(biom)
  ranks  <- taxa_ranks(biom)
  labels <- tree$tip.label
  nTips  <- length(labels)
  df     <- data.frame(node = seq_len(max(tree$edge)))
  
  
  # Convert `clades` param to named `ranks`.
  #________________________________________________________
  ranks <- local({
    
    stopifnot(is.logical(clades) || is.character(clades))
    
    if (isFALSE(clades)) return (NULL)
    if (isTRUE(clades))  return (ranks %>% setNames(., .))
    
    if (length(bad <- setdiff(clades, ranks)) > 0)
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
    df[1:nTips, names(ranks)[[i]]] <- as.vector(otu_taxonomy(biom, ranks[[i]])[labels])
  
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
  tree <- tidytree::full_join(tree, df, by = "node")
  
  set_cache_value(cache_file, tree)
  return (tree)
}


#' Names of a phylogenetic tree's tips/leafs.
#' 
#' @family phylogeny
#' 
#' @param tree  A phylo object, as returned from [read_tree()].
#' @return A character vector with the leaf names.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "newick.tre", package = "rbiom")
#'     tree  <- read_tree(infile)
#'     leafs <- tree_tips(tree)
#'     head(leafs)
#'
tree_tips <- function (tree) {
  stopifnot(is(tree, 'phylo'))
  tree$tip.label
}



#' Create a subtree by specifying tips to keep.
#' 
#' @name tree_subset
#' 
#' @family phylogeny
#' 
#' @param tree  A phylo object, as returned from [read_tree()].
#' @param tips  A character, numeric, or logical vector of tips to keep.
#' @return A \code{phylo} object for the subtree.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "newick.tre", package = "rbiom")
#'     tree <- read_tree(infile)
#'     tree
#'     
#'     subtree <- tree_subset(tree, tips = head(tree_tips(tree)))
#'     subtree
#'
tree_subset <- function (tree, tips) {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- eval_envir(environment())
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Did we get a tree?
  #________________________________________________________
  
  if (!is(tree, "phylo"))
    stop("Must provide a phylo object as the tree.")
  
  nTips <- length(setdiff(tree$edge[,2], tree$edge[,1]))
  
  
  #________________________________________________________
  # Tips specified as strings c("SteSpe52", "HrbSpe65", ...)
  #________________________________________________________
  
  if (is.character(tips)) {
    
    if (!all(tips %in% tree$tip.label)) {
      missing <- tips[!tips %in% tree$tip.label]
      if (length(missing) > 5) missing <- c(missing[1:5], "...")
      missing <- paste(collapse=", ", missing)
      stop(sprintf("Tips missing from tree: %s", missing))
    }
    
    tips <- which(tree$tip.label %in% tips)
  }
  
  
  #________________________________________________________
  # Tips specified as logicals c(TRUE, TRUE, FALSE, TRUE, ...)
  #________________________________________________________
  
  if (is.logical(tips)) {
    
    if (length(tips) != nTips)
      stop("Logical vector must be same length as tips.")
    
    tips <- which(tips)
  }
  
  
  #________________________________________________________
  # Tips specified as numeric c(5, 1, 2, 31, 3, 12, ...)
  #________________________________________________________
  
  if (is.numeric(tips)) {
    if (max(tips) > nTips)  stop("There aren't that many tips in the tree.")
    if (min(tips) < 1)      stop("Tips with index < 1 are not allowed.")
    if (any(tips %% 1 > 0)) stop("Tip indices must be whole numbers.")
    
  } else {
    stop("Tips must be given as character, logical, or numeric.")
  }
  
  
  #________________________________________________________
  # Keeping zero or all tips?
  #________________________________________________________
  
  tips <- unique(tips)
  
  if (length(tips) < 2)
    stop("At least two tips must be specified.")
  
  if (length(tips) == nTips) {
    set_cache_value(cache_file, tree)
    return (tree)
  }
  
  
  #________________________________________________________
  # Prune the tree
  #________________________________________________________
  
  root      <- setdiff(tree$edge[,1], tree$edge[,2])
  childAt   <- order(c(tree$edge[,2], root))
  nChildren <- c(rep(0, nTips), unname(table(tree$edge[,1])))
  
  eLength   <- tree$edge.length
  if (is_null(eLength))
    eLength <- rep(0, nrow(tree$edge))
  
  for (tip in setdiff(1:nTips, tips)) repeat {
    
    tipIdx <- childAt[tip]
    parent <- tree$edge[tipIdx,1]
    
    nChildren[parent] <- nChildren[parent] - 1
    eLength[tipIdx] <- NA
    
    if (nChildren[parent] > 0) break
    tip <- parent
  }
  
  
  #________________________________________________________
  # Merge branches
  #________________________________________________________
  
  for (tip in tips) repeat {
    
    if (tip == root) break
    
    tipIdx <- childAt[tip]
    parent <- tree$edge[tipIdx,1]
    
    if (nChildren[parent] == 1 && parent != root) {
      
      parentIdx                   <- childAt[parent]
      tree$edge[parentIdx, 2]     <- tree$edge[tipIdx, 2]
      eLength[parentIdx] <- eLength[parentIdx] + eLength[tipIdx]
      eLength[tipIdx]    <- NA
      
      nChildren[parent] <- 0
    }
    
    tip <- parent
  }
  
  
  #________________________________________________________
  # Move the root
  #________________________________________________________
  
  if (nChildren[root] == 1)
    eLength[which(tree$edge[,1] == root)] <- NA
  
  
  #________________________________________________________
  # Discard nodes flagged with NA; make the new subtree
  #________________________________________________________
  
  tree$edge  <- tree$edge[!is.na(eLength),]
  tree$Nnode <- length(unique(tree$edge[,1]))
  
  if (!is_null(tree$edge.length)) tree$edge.length <- eLength[!is.na(eLength)]
  if (!is_null(tree$tip.label))   tree$tip.label   <- tree$tip.label[tips]
  if (!is_null(tree$node.label))  tree$node.label  <- tree$tip.label[unique(sort(tree$edge[,1]))]
  
  tree$edge <- matrix(as.numeric(factor(as.vector(tree$edge))), ncol=2)
  
  
  set_cache_value(cache_file, tree)
  return (tree)
}




