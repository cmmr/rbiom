# Needs to be re-written to not use tidytree package, which causes this mess of 
# warnings quite often:
# #> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
# #> Also defined by 'tidytree'


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
#' @param n   Number of (most abundant) OTUs to include. Set to `Inf` to see 
#'        all OTUs. Default: `50`
#' 
#' @param layout   Any layout option supported by [ggtree::ggtree()]:
#'        \code{"rectangular"}, \code{"dendrogram"}, \code{"slanted"}, 
#'        \code{"ellipse"}, \code{"roundrect"}, \code{"fan"}, 
#'        \code{"circular"}, \code{"inward_circular"}, \code{"radial"}, 
#'        \code{"equal_angle"}, \code{"daylight"} or \code{"ape"}.
#'        Default: `"rectangular"`
#' 
#' @param color.by   How to color the tree. Currently only supports the options 
#'        `NULL` (no coloring) or `".mean"` (color by number of OTU 
#'        observations). Default: `".mean"`
#' 
#' @param tip   Label the tree leafs with taxa names. Options are `NULL` (no 
#'        tip labels) or a taxonomic rank (see `biom$ranks` for options). 
#'        Default: `".otu"`
#' 
#' @param align   Left-align tip labels. Default: `FALSE`
#' 
#' @param clade,branch   Label monophyletic clades externally (clade) or 
#'        internally (branch). Options are a rank listed in `biom$ranks` or an 
#'        index thereof. Default: `clade = 2, branch = NULL`
#' 
#' @param ab_color,ab_size   Adjust colors and/or sizes of tree components to 
#'        reflect OTU abundances. Default: `ab_color = TRUE, ab_size = FALSE`
#' 
#' @param margin_top,margin_right,margin_bottom,margin_left   Add additional 
#'        space around the tree. Sometimes necessary for wide text annotations. 
#'        Set as fraction of tree's width/height. For instance, 
#'        \code{right = 1} reserves the right half of the plotting area for 
#'        non-tree elements. Default: `NULL`
#' 
#' @examples
#'     library(rbiom)
#'     
#'     # Needs the optional ggtree package
#'     if (nzchar(system.file(package = "ggtree")))
#'       tree_plot(hmp50)
#'       
#'     if (nzchar(system.file(package = "ggtree")))
#'       tree_plot(hmp50, layout = "roundrect")
#'     
#'     if (nzchar(system.file(package = "ggtree")))
#'       tree_plot(hmp50, tip = "Genus", layout = "fan")

tree_plot <- function (
    biom, n = 25, layout = 'rectangular', 
    tip = ".otu", clade = 2, branch = NULL, align = FALSE, 
    ab_color = TRUE, ab_size = FALSE, limit.by = NULL, 
    unc = 'singly', caption = TRUE, 
    margin_top = NULL, margin_right = NULL, 
    margin_bottom = NULL, margin_left = NULL, ...) {
  
  if (nchar(system.file(package = "ggtree")) == 0)
    stop("Bioconductor R package 'ggtree' must be installed to use tree_plot().")
  
  biom <- as_rbiom(biom)
  if (is.null(biom$tree))
    cli_abort("biom must include a phylogentic tree")
  
  
  params <- eval_envir(environment(), ...)
  cmd    <- sprintf("tree_plot(%s)", as.args(params, 2, tree_plot))
  remove(list = intersect(env_names(params), ls()))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('tree_plot', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  
  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  with(params, {
    
    validate_bool('caption')
    validate_bool('align')
    validate_bool('ab_color')
    validate_bool('ab_size')
    
    validate_var_range('n', range = c(1, Inf), n = 1)
    
    validate_rank('tip',    null_ok = TRUE)
    validate_rank('clade',  null_ok = TRUE)
    validate_rank('branch', null_ok = TRUE)
    
    if (eq(clade,  ".otu")) { tip <- ".otu"; clade  <- NULL }
    if (eq(branch, ".otu")) { tip <- ".otu"; branch <- NULL }
    
    validate_meta_aes('limit.by', null_ok = TRUE, max = Inf)
    
    sync_metadata() # Also clones biom
  })
  
  
  
  #________________________________________________________
  # Drop low abundance OTUs.
  #________________________________________________________
  with(params, {
    
    biom$counts %<>% rescale_cols()
    
    if (n >= biom$n_otus) {
      
      if (isTRUE(caption))
        caption <- glue("Displaying all {biom$n_otus} OTUs.")
      
    } else {
      
      if (isTRUE(caption))
        caption <- glue("Displaying top {n} of {biom$n_otus} OTUs.")
      
      biom$counts <- biom$counts[1:n,] # Internally sorted by row sums
    }
  })
  
  
  
  #________________________________________________________
  # Generate the treedatato plot.
  #________________________________________________________
  with(params, {
    .tr <- tree_data(
      biom        = biom, 
      unc         = unc, 
      data_ranks  = setdiff(tip, ".otu"), 
      extra_ranks = unique(c(branch, clade)) )
  })
  
  
  #________________________________________________________
  # Initial call to ggtree() with tree data.
  #________________________________________________________
  p <- with(params, {
    ggtree(
      .indent = 4, 
      layout  = layout,
      tr      = .tr,
      mapping = local({
        if (ab_color && ab_size) return (aes(color = .mean, size = .mean))
        if (ab_color)            return (aes(color = .mean))
        if (ab_size)             return (aes(size = .mean))
        return (NULL)
      }) )
  })
  
  cmds <- attr(p, 'display', exact = TRUE)
  
  
  #________________________________________________________
  # Add OTU or other rank names at the tree tips.
  #________________________________________________________
  if (length(params$tip) > 0) {
    
    if (eq(params$tip, '.otu')) {
      args <- list()
      
    } else {
      args <- list(
        .indent = 4,
        data    = as.cmd(function (x) { subset(x, !is.na(y)) }, list(y = as.name(params$tip))), 
        mapping = aes(label = !!as.name(params$tip)) )
    }
    
    args[['align']] <- if (isTRUE(params$align)) TRUE else NULL
    
    x    <- do.call(geom_tiplab, args)
    p    <- p + x
    cmds <- c(cmds, attr(x, 'display', exact = TRUE))
  }
  
  
  
  #________________________________________________________
  # Label clades internally (inside) of the tree.
  #________________________________________________________
  if (length(params$branch) > 0) {
    
    if (isTRUE(params$ab_color)) {
      mapping <- aes(
        x     = branch, 
        label = !!as.name(params$branch), 
        color = .mean )
      
    } else {
      mapping <- aes(
        x     = branch, 
        label = !!as.name(params$branch) )
    }
    
    x <- geom_label(
      .indent = 4, 
      data    = as.cmd(function (x) { subset(x, !is.na(y)) }, list(y = as.name(params$branch))), 
      mapping = mapping, 
      fill    = "white" )
    
    p    <- p + x
    cmds <- c(cmds, attr(x, 'display', exact = TRUE))
  }
  
  
  #________________________________________________________
  # Label clades externally (at right) of the tree.
  #________________________________________________________
  if (length(params$clade) > 0) {
    
    args <- list(
      .indent     = 4, 
      data        = as.cmd(function (x) { subset(x, !is.na(y)) }, list(y = as.name(params$clade))), 
      offset      = ifelse(length(params$tip) == 0, 0.02, 0.2),
      offset.text = 0.02 )
    
    if (isTRUE(params$ab_color)) {
      args[['mapping']] <- aes(
        node  = node, 
        label = !!as.name(params$clade), 
        color = .mean )
      
    } else {
      args[['mapping']] <- aes(
        node  = node, 
        label = !!as.name(params$clade) )
    }
    
    args[['align']] <- if (isTRUE(params$align)) TRUE else NULL
    
    x    <- do.call(geom_cladelab, args)
    p    <- p + x
    cmds <- c(cmds, attr(x, 'display', exact = TRUE))
  }
  
  
  
  #________________________________________________________
  # Use a yellow-orange-black gradient abundance scale.
  #________________________________________________________
  if (isTRUE(params$ab_color)) {
    x    <- scale_color_gradientn(
      .indent = 4,
      name    = "Mean\nAbundance", 
      colors  = get_palette("lajolla")[c(3,6,7,8)],
      labels  = as.cmd(scales::percent) )
    
    p    <- p + x
    cmds <- c(cmds, attr(x, 'display', exact = TRUE))
  }
  
  
  
  #________________________________________________________
  # Add caption to the plot.
  #________________________________________________________
  if (!isFALSE(params$caption)) {
    x    <- labs(caption = params$caption)
    p    <- p + x
    cmds <- c(cmds, attr(x, 'display', exact = TRUE))
  }
  
  
  
  #________________________________________________________
  # Faceting. Doesn't work for geom_cladelab.
  # See: https://github.com/YuLab-SMU/ggtree/issues/603
  #________________________________________________________
  # if (!is.null(params$facet.by)) {
  #   x    <- facet_wrap( ~.id, scale="free")
  #   p    <- p + x
  #   cmds <- c(cmds, attr(x, 'display', exact = TRUE))
  # }
  
  
  #________________________________________________________
  # Add horizontal/vertical padding for text to bleed into.
  #________________________________________________________
  expand <- with(params, {
    
    if (length(clade) > 0) margin_right %<>% if.null(0.2)
    if (length(tip)   > 0) margin_right %<>% if.null(0.1)
    
    if (length(c(branch, clade, tip)) > 0) {
      margin_top    %<>% if.null(0.02)
      margin_bottom %<>% if.null(0.02)
    }
    
    list(
      if (is_null(margin_top))    NULL else vexpand(ratio = margin_top,    direction =  1),  
      if (is_null(margin_right))  NULL else hexpand(ratio = margin_right,  direction =  1),  
      if (is_null(margin_bottom)) NULL else vexpand(ratio = margin_bottom, direction = -1),  
      if (is_null(margin_left))   NULL else hexpand(ratio = margin_left,   direction = -1) )
  })
  
  for (x in expand[!sapply(expand, is.null)]) {
    p    <- p + x
    cmds <- c(cmds, attr(x, 'display', exact = TRUE))
  }
  
  
  
  #______________________________________________________________
  # Finalize custom plot attributes
  #______________________________________________________________
  attr(p, 'code') <- paste0(
      "library(ggtree)\n\n", 
      paste(cmds, collapse=" +\n  ") ) %>% 
    add_class('rbiom_code')
  
  p %<>% add_class('rbiom_plot')
  p$plot_env <- emptyenv()
  
  
  
  set_cache_value(cache_file, p)
  
  return (p)
}



#' Provides a 'treedata' S4 object for use in ggtree functions.
#' 
#' @noRd
#' @keywords internal
#' 
#' @param biom   An rbiom object, as returned from [read_biom()].
#' 
#' @param data_ranks,extra_ranks   Map interior and leaf nodes to these ranks, 
#'        populating @data or @extraInfo in the returned object. Options are
#'        `biom$ranks`, or `".all"` for all ranks. Default: `NULL`
#' 
#' @examples
#'     library(rbiom)
#'     
#'     # Needs the optional tidytree package
#'     if (nzchar(system.file(package = "tidytree"))) {
#'     
#'       tree_data(hmp50)
#'       
#'     }

tree_data <- function (biom, unc = 'singly', data_ranks = NULL, extra_ranks = NULL) {
  
  if (nchar(system.file(package = "ggtree")) == 0)
    stop("Bioconductor R package 'ggtree' must be installed to use tree_data().")
  
  biom <- as_rbiom(biom)
  if (is.null(biom$tree))
    cli_abort("biom must include a phylogentic tree")
  
  params <- eval_envir(environment())
  cmd    <- sprintf("rbiom:::tree_data(%s)", as.args(params, 0, tree_data))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('tree_data', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  validate_rank('data_ranks',  null_ok = TRUE, max = Inf, all_option = ".all")
  validate_rank('extra_ranks', null_ok = TRUE, max = Inf, all_option = ".all")
  
  
  tree   <- biom$tree
  labels <- tree$tip.label
  nTips  <- length(labels)
  df     <- data.frame(node = seq_len(max(tree$edge)))
  
  
  
  # We can pre-populate all the data for leaf nodes.
  #________________________________________________________
  df[1:nTips, '.mean'] <- as.vector(taxa_means(biom)[labels])
  
  for (rank in extra_ranks)
    df[1:nTips, rank] <- as.vector(taxa_map(biom, rank, unc)[labels])
  
  
  
  # Recurse through the tree.
  #________________________________________________________
  traverse <- function (i) {
    
    if (i <= nTips) return (i)
    
    children    <- tree$edge[which(tree$edge[,1] == i),2]
    descendents <- unlist(lapply(children, traverse))
    
    df[i, '.mean'] <<- mean(df[children, '.mean'])
    
    for (rank in extra_ranks) {
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
  
  
  if (length(data_ranks) > 0) {
    
    data  <- attr(tree, 'data',  exact = TRUE)
    phylo <- attr(tree, 'phylo', exact = TRUE)
    
    if (nrow(data) == 0)
      data <- attr(tree, 'extraInfo', exact = TRUE)[,'node']
    
    for (rank in data_ranks)
      data[[rank]] <- as.vector(taxa_map(biom, rank, unc)[phylo$tip.label][data[['node']]])
    
    attr(tree, 'data') <- data
  }
  
  
  attr(tree, 'cmd')     <- cmd
  attr(tree, 'display') <- cmd
  
  set_cache_value(cache_file, tree)
  return (tree)
}




