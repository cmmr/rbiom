# #' Plot the phylogenetic reference tree.
# #' 
# #' @noRd
# #' @keywords internal
# #' 
# #' @inherit documentation_default
# #' @inherit documentation_rank.2
# #' @inherit documentation_plot_return
# #' 
# #' @family taxonomy
# #' @family phylogeny
# #' @family visualization
# #' 
# #' @param layout  Tree style. Options are: `'unrooted'`, `'dendrogram'`, or 
# #'        `'circular'` (or abbreviations thereof). Default: `'unrooted'`
# #' 
# #' @param fixed  Set all branch lengths to 1. If `FALSE`, uses branch lengths 
# #'        from the reference tree, if present. Default: `TRUE`
# #' 
# #' @param points  Represent OTU abundances as circles. Default: `TRUE`
# #' 
# #' @param drop  Drop any taxa not matched by `taxa`. Default: `FALSE`
# #' 
# #' @examples
# #'     library(rbiom)
# #'     
# #'     tree_plot(hmp50)
# #'     
# #'     tree_plot(hmp50, layout = 'circular', fixed = FALSE)
# #'
# 
# tree_plot <- function (
#     biom, rank = 2, taxa = 6,
#     layout = 'unrooted', fixed = TRUE, points = TRUE, colors = 'okabe',
#     unc = "singly", lineage = FALSE, drop = FALSE, tree = NULL ) {
# 
# 
#   #________________________________________________________
#   # Take care not to cache filepath to tree.
#   #________________________________________________________
#   biom <- as_rbiom(biom)
#   validate_tree(null_ok = TRUE)
# 
#   params <- eval_envir(environment())
#   cmd    <- sprintf("tree_plot(%s)", as.args(params, fun = tree_plot))
# 
# 
#   #________________________________________________________
#   # See if this result is already in the cache.
#   #________________________________________________________
#   cache_file <- get_cache_file('tree_plot', params)
#   # if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
#   #   return (readRDS(cache_file))
#   remove("params")
# 
# 
#   #________________________________________________________
#   # Sanity checks
#   #________________________________________________________
#   if (!is.null(tree))     biom$tree <- tree
#   if (is.null(biom$tree)) cli_abort("biom must include a phylogentic tree")
# 
#   validate_rank()
#   validate_taxa(null_ok = TRUE)
#   validate_var_choices('layout', c('unrooted', 'dendrogram', 'circular'))
#   validate_bool('fixed')
#   validate_bool('points')
#   validate_bool('drop')
# 
# 
#   #________________________________________________________
#   # Draw a hierarchical tree using ggraph.
#   #________________________________________________________
# 
# 
#   # OTU colors (subset OTUs)
#   if (is.null(taxa)) {
# 
#     data <- tidygraph::as_tbl_graph(biom$tree)
# 
#   } else {
# 
#     otu_colors <- taxa_map(biom, rank, taxa, unc, lineage, other = FALSE)
# 
#     if (isTRUE(drop))
#       biom$counts <- biom$counts[names(otu_colors),]
# 
#     data <- tidygraph::as_tbl_graph(biom$tree)
# 
#     data <- local({
# 
#       node_name <- pull(tidygraph::activate(data, nodes), name)
#       edge_from <- pull(tidygraph::activate(data, edges), from)
#       edge_to   <- pull(tidygraph::activate(data, edges), to)
#       node_taxa <- unname(otu_colors[node_name])
#       edge_taxa <- factor(rep(NA, length(edge_to)), levels = levels(node_taxa))
# 
#       f <- function (node) {
#         taxa <- sapply(edge_to[edge_from == node], f)
#         if (node == 1) return (invisible(NULL))
#         if (length(unique(taxa)) == 0) taxa <- node_taxa[node]
#         if (length(unique(taxa)) >= 2) taxa <- NA
#         edge_taxa[edge_to == node] <<- taxa[[1]]
#         return (taxa)
#       }
#       f(1)
# 
#       data %>%
#         tidygraph::activate(nodes) %>% mutate(taxa = node_taxa) #%>%
#         #tidygraph::activate(edges) %>% mutate(taxa = edge_taxa)
# 
#     })
# 
#   }
# 
# 
#   # OTU sizes
#   if (isTRUE(points)) {
# 
#     otu_sizes <- taxa_sums(biom, rank = 0, unc = 'asis')
# 
#     data <- data %>%
#       tidygraph::activate(nodes) %>%
#       mutate(size = unname(otu_sizes[name]))
#   }
# 
# 
#   # Layout
#   if (layout == 'circular') {
# 
#     if (fixed) {
#       p    <- ggraph(data, 'dendrogram', circular = TRUE) + ggplot2::coord_fixed()
#       code <- c("ggraph::ggraph(data, 'dendrogram', circular = TRUE)", "ggplot2::coord_fixed()")
#     } else {
#       p    <- ggraph(data, 'dendrogram', length = length, circular = TRUE) + ggplot2::coord_fixed()
#       code <- c("ggraph::ggraph(data, 'dendrogram', length = length, circular = TRUE)", "ggplot2::coord_fixed()")
#     }
# 
#   } else {
# 
#     if (fixed) {
#       p    <- ggraph(data, layout)
#       code <- glue("ggraph::ggraph(data, '{layout}')")
#     } else {
#       p    <- ggraph(data, layout, length = length)
#       code <- glue("ggraph::ggraph(data, '{layout}', length = length)")
#     }
# 
#   }
# 
# 
#   # Edges
#   if (layout == 'unrooted') {
# 
#     if (is.null(taxa)) {
#       p <- p + geom_edge_link()
#       code %<>% c("ggraph::geom_edge_link()")
#     } else {
#       p <- p + geom_edge_link(ggplot2::aes(color = taxa))
#       code %<>% c("ggraph::geom_edge_link(ggplot2::aes(color = taxa))")
#     }
# 
#   } else {
# 
#     if (is.null(taxa)) {
#       p <- p + geom_edge_elbow()
#       code %<>% c("ggraph::geom_edge_elbow()")
#     } else {
#       p <- p + geom_edge_elbow2(ggplot2::aes(color = node.taxa))
#       code %<>% c("ggraph::geom_edge_elbow2(ggplot2::aes(color = node.taxa))")
#     }
# 
#   }
# 
# 
#   # Tips
#   if (points) {
# 
#     if (is.null(taxa)) {
# 
#       p <- p + geom_node_point(
#         mapping     = ggplot2::aes(size = size),
#         data        = ~ filter(.x, leaf),
#         alpha       = 0.4,
#         show.legend = FALSE )
# 
#       code %<>% c(paste(collapse = "\n    ", c(
#         "ggraph::geom_node_point(",
#         "mapping     = ggplot2::aes(size = size),",
#         "data        = ~ filter(.x, leaf),",
#         "alpha       = 0.4,",
#         "show.legend = FALSE )" )))
# 
#     } else {
# 
#       p <- p + geom_node_point(
#         mapping     = ggplot2::aes(size = size, color = taxa),
#         data        = ~ filter(.x, leaf),
#         alpha       = 0.4,
#         show.legend = FALSE )
# 
#       code %<>% c(paste(collapse = "\n    ", c(
#         "ggraph::geom_node_point(",
#         "mapping     = ggplot2::aes(size = size, color = taxa),",
#         "data        = ~ filter(.x, leaf),",
#         "alpha       = 0.4,",
#         "show.legend = FALSE )" )))
#     }
# 
#   }
# 
# 
#   if (!is.null(taxa)) {
# 
#     palette <- color_palette(pal = colors, keys = pull(data, taxa))
#     p <- p + ggplot2::scale_colour_manual(values = palette, na.value = '#000000')
# 
#     code %<>% c(paste(collapse = "\n", c(
#       "ggplot2::scale_colour_manual(",
#       "    values   = c(%s)," %>% sprintf(as.args(as.list(palette), indent = 6L)),
#       "    na.value = '#000000' )" )))
#   }
# 
# 
#   code <- paste(collapse = "", c(
#     "library(ggraph)\n",
#     # "biom <- ", biom$src, "\n",
#     "plot <- ", cmd, "\n",
#     "data <- plot$data\n",
#     paste(code, collapse = " +\n  ")
#   )) %>% add_class('rbiom_code')
# 
#   attr(p, 'cmd')  <- cmd
#   attr(p, 'code') <- code
#   p <- add_class(p, 'rbiom_plot')
# 
#   # set_cache_value(cache_file, p)
#   return (p)
# }
# 
# 
