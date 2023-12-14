# Needs to be re-written to not use tidytree package, which causes this mess of 
# warnings quite often:
# #> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
# #> Also defined by 'tidytree'


# #' Display a dendrogram of the phylogenetic tree.
# #' 
# #' @noRd
# #' @keywords internal
# #' 
# #' @inherit documentation_default
# #' 
# #' @family phylogeny
# #' @family visualization
# #' 
# #' @param layout   Any layout option supported by [ggtree::ggtree()]:
# #'        \code{"rectangular"}, \code{"dendrogram"}, \code{"slanted"}, 
# #'        \code{"ellipse"}, \code{"roundrect"}, \code{"fan"}, 
# #'        \code{"circular"}, \code{"inward_circular"}, \code{"radial"}, 
# #'        \code{"equal_angle"}, \code{"daylight"} or \code{"ape"}.
# #'        Default: \code{"rectangular"}.
# #' 
# #' @param tiplab   Label the tree leafs with taxa names. Options are
# #'        `NULL` (no tip labels) or a taxonomic rank (for example
# #'        \code{tiplab = "Genus"}). Default: `NULL`.
# #' 
# #' @param color.by   How to color the tree. Currently only supports the options
# #'        `NULL` (no coloring) or \code{'.reads'} (color by number of taxa 
# #'        observations). Default: `NULL`.
# #' 
# #' @param label,cladelab   Label monophyletic clades. You can specify the
# #'        same or different taxonomic ranks for internal (label) and external 
# #'        (cladelab) annotations. Default: `NULL`.
# #' 
# #' @param top,right,bottom,left   Add additional space around the tree.
# #'        Sometimes necessary for wide text annotations. Set as fraction of
# #'        tree's width/height. For instance, \code{right = 1} reserves the 
# #'        right half of the plotting area for non-tree elements.
# #'        Default: `NULL`.
# #' 
# #' @examples
# #'     library(rbiom)
# #'     
# #'     hmp1        <- hmp50$clone()
# #'     hmp1$counts <- hmp1$counts[,"HMP1"]
# #'     
# #'     # Needs the optional ggtree package
# #'     if (nzchar(system.file(package = "ggtree")))
# #'       tree_plot(hmp50)
# #'       
# #'     if (nzchar(system.file(package = "ggtree")))
# #'       tree_plot(hmp1, cladelab = "Phylum", layout = "roundrect")
# #'     
# #'     if (nzchar(system.file(package = "ggtree")))
# #'       tree_plot(hmp1, tiplab = "Genus", layout = "fan", color.by = ".reads")
# 
# tree_plot <- function (
#     biom, layout = "rectangular", 
#     tiplab = NULL, color.by = NULL, label = NULL, cladelab = NULL, 
#     top = NULL, right = NULL, bottom = NULL, left = NULL, ...) {
#   
#   if (nchar(system.file(package = "ggtree")) == 0)
#     stop("Bioconductor R package 'ggtree' must be installed to use tree_plot().")
#   
#   biom <- as_rbiom(biom)
#   
#   
#   #________________________________________________________
#   # See if this result is already in the cache.
#   #________________________________________________________
#   params     <- eval_envir(environment(), ...)
#   cache_file <- get_cache_file()
#   if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
#     return (readRDS(cache_file))
#   
#   
#   #________________________________________________________
#   # Sanity checks
#   #________________________________________________________
#   stopifnot(!is.null(biom$tree))
#   validate_rank("tiplab",   null_ok = TRUE)
#   validate_rank("label",    null_ok = TRUE)
#   validate_rank("cladelab", null_ok = TRUE)
#   
#   
#   
#   #________________________________________________________
#   # Assemble the plot layers with init_layers / set_layer.
#   #________________________________________________________
#   layers <- list()
#   attr(layers, 'params')   <- params
#   attr(layers, 'function') <- tree_plot
#   
#   
#   
#   #________________________________________________________
#   # Use ggtree to render a phylogenetic tree
#   #________________________________________________________
#   tr <- tree_data(biom)
#   attr(tr, 'display') <- "tree_data(biom)"
#   set_layer(params, 'ggtree', layout = layout, tr = tr)
#   
#   
#   #________________________________________________________
#   # Coloring by reads requires an aes mapping.
#   #________________________________________________________
#   if (eq(color.by, '.reads')) {
#     set_layer(params, 'ggtree', mapping = as.cmd(aes(color=log(reads))))
#     set_layer(params, 'scale_color_continuous', low='darkgreen', high='red')
#   }
#   
#   
#   #________________________________________________________
#   # Add OTU names at the tree tips.
#   #________________________________________________________
#   if (!is_null(tiplab)) {
#     
#     add_layer(params, 'tiplab')
#     
#     if (!tiplab == '.otu') {
#       
#       df  <- attr(layers[['ggtree']][['tr']], 'data',  exact = TRUE)
#       phy <- attr(layers[['ggtree']][['tr']], 'phylo', exact = TRUE)
#       
#       if (nrow(df) == 0)
#         df <- attr(tr, 'extraInfo', exact = TRUE)[,'node']
#       
#       df[['tiplab']] <- taxa_map(biom, tiplab)[phy$tip.label][df[['node']]]
#       
#       attr(layers[['ggtree']][['tr']], 'data') <- df
#       remove("df", "phy")
#       
#       set_layer(
#         params = params, 
#         layer  = 'tiplab', 
#         'data'    = as.cmd(function (x) { subset(x, !is.na(tiplab)) }), 
#         'mapping' = aes(label=tiplab) )
#     }
#   }
#   
#   
#   #________________________________________________________
#   # Label clades internally (inside) of the tree.
#   #________________________________________________________
#   if (!is_null(label)) {
#     
#     mapping <- aes(x=branch, label=!!as.name(label))
#     if (eq(color.by, '.reads'))
#       mapping <- aes(x=branch, label=!!as.name(label), color=log(reads))
#     
#     set_layer(
#       params = params, 
#       layer  = 'label',
#       'data'    = as.cmd(function (x) { subset(x, !is.na(y)) }, list(y = as.name(label))), 
#       'mapping' = mapping, 
#       'fill'    = 'white' )
#   }
#   
#   
#   #________________________________________________________
#   # Label clades externally (at right) of the tree.
#   #________________________________________________________
#   if (!is_null(cladelab)) {
#     
#     mapping <- aes(node=node, label=!!as.name(cladelab))
#     if (eq(color.by, '.reads'))
#       mapping <- aes(node=node, label=!!as.name(cladelab), color=log(reads))
#     
#     set_layer(
#       params = params, 
#       layer  = 'cladelab',
#       'data'        = as.cmd(function (x) { subset(x, !is.na(y)) }, list(y = as.name(cladelab))), 
#       'mapping'     = mapping,
#       'offset'      = 0.02,
#       'offset.text' = 0.02 )
#   }
#   
#   
#   #________________________________________________________
#   # Add horizontal/vertical padding for text to bleed into.
#   #________________________________________________________
#   if (!is_null(c(label, cladelab, tiplab))) {
#     top    %<>% if.null(0.02)
#     bottom %<>% if.null(0.02)
#   }
#   if (!is_null(cladelab))
#     right %<>% if.null(0.2)
#   
#   if (!is_null(top))    set_layer(params, 'top',    .fn = 'vexpand', direction =  1, ratio = top)
#   if (!is_null(right))  set_layer(params, 'right',  .fn = 'hexpand', direction =  1, ratio = right)
#   if (!is_null(bottom)) set_layer(params, 'bottom', .fn = 'vexpand', direction = -1, ratio = bottom)
#   if (!is_null(left))   set_layer(params, 'left',   .fn = 'hexpand', direction = -1, ratio = left)
#   
#   
#   p <- suppressMessages(ggbuild(layers))
#   
#   
#   set_cache_value(cache_file, p)
#   
#   return (p)
# }
# 
# 
# 
# #' Provides a 'treedata' S4 object for use in ggtree functions.
# #' 
# #' @noRd
# #' @keywords internal
# #' 
# #' @param biom   An rbiom object, as returned from [read_biom()].
# #' 
# #' @param reads   Include a 'reads' column indicating the sum of taxa 
# #'        observations belonging to each node/leaf. Default: `TRUE`.
# #' 
# #' @param clades   Notate top-most monophyletic nodes. The default,
# #'        `TRUE`, adds a column for every rank in the biom object. A
# #'        character vector of ranks can also be passed in. If the vector is
# #'        named, then those names are used for naming the columns in the
# #'        returned treedata object. Set to `NULL` to not return any 
# #'        clade notations.
# #' 
# #' @examples
# #'     library(rbiom)
# #'     
# #'     # Needs the optional ggtree package
# #'     if (nzchar(system.file(package = "ggtree"))) {
# #'     
# #'       tree_data(hmp50)
# #'       
# #'     }
# 
# tree_data <- function (biom, reads = TRUE, clades = TRUE) {
#   
#   if (nchar(system.file(package = "ggtree")) == 0)
#     stop("Bioconductor R package 'ggtree' must be installed to use tree_data().")
#   
#   biom <- as_rbiom(biom)
#   
#   
#   #________________________________________________________
#   # See if this result is already in the cache.
#   #________________________________________________________
#   params     <- eval_envir(environment())
#   cache_file <- get_cache_file()
#   if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
#     return (readRDS(cache_file))
#   
#   
#   # Sanity checks
#   #________________________________________________________
#   stopifnot(is_scalar_logical(reads) && !is_na(reads))
#   
#   tree   <- biom$tree
#   ranks  <- biom$ranks
#   labels <- tree$tip.label
#   nTips  <- length(labels)
#   df     <- data.frame(node = seq_len(max(tree$edge)))
#   
#   
#   # Convert `clades` param to named `ranks`.
#   #________________________________________________________
#   ranks <- local({
#     
#     stopifnot(is.logical(clades) || is.character(clades))
#     
#     if (isFALSE(clades)) return (NULL)
#     if (isTRUE(clades))  return (ranks %>% setNames(., .))
#     
#     if (length(bad <- setdiff(clades, ranks)) > 0)
#       stop("Invalid taxonomic rank: ", paste(collapse = ", ", bad))
#       
#     if (is_null(names(clades))) return (setNames(clades, clades))
#     
#     if (length(i <- which(names(clades) == "")) > 0)
#       names(clades)[i] <- as.vector(clades)[i]
#     
#     return (clades)
#   })
#   
#   
#   # We can pre-populate all the data for leaf nodes.
#   #________________________________________________________
#   if (isTRUE(reads))
#     df[1:nTips, 'reads'] <- as.vector(taxa_sums(biom)[labels])
#   
#   for (i in seq_along(ranks))
#     df[1:nTips, names(ranks)[[i]]] <- as.vector(taxa_map(biom, ranks[[i]])[labels])
#   
#   df[1:nTips, 'OTU'] <- labels
#   
#   
#   # Recurse through the tree.
#   #________________________________________________________
#   traverse <- function (i) {
#     
#     if (i <= nTips) return (i)
#     
#     children    <- tree$edge[which(tree$edge[,1] == i),2]
#     descendents <- unlist(lapply(children, traverse))
#     
#     if (isTRUE(reads))
#       df[i, 'reads'] <<- sum(df[children, 'reads'])
#     
#     for (rank in names(ranks)) {
#       if (length(unique(df[children, rank])) == 1) {
#         df[i, rank]        <<- df[children[[1]], rank]
#         df[children, rank] <<- NA
#       }
#     }
#     
#     return (c(i, descendents))
#   }
#   invisible(traverse(nTips + 1))
#   
#   
#   # Convert phylo object to treedata object.
#   #________________________________________________________
#   tree <- tidytree::full_join(tree, df, by = "node")
#   
#   set_cache_value(cache_file, tree)
#   return (tree)
# }




