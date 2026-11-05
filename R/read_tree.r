#' Read a newick formatted phylogenetic tree.
#'
#' A phylogenetic tree is required for computing UniFrac distance matrices.
#' You can load a tree from a file or by providing the tree string directly.
#' This tree must be in Newick format, also known as parenthetic format and
#' New Hampshire format.
#'
#' @family phylogeny
#'
#' @inherit documentation_default
#'
#' @param src   Input data as either a file path, URL, or Newick string.
#'        Compressed (gzip or bzip2) files are also supported.
#'
#' @details
#' Labels may be enclosed in single quotes, following the Newick standard:
#' inside quotes every character is literal (including `( ) [ ] : ; ,` and
#' whitespace), and `''` stands for one apostrophe. Quoted labels are never
#' altered by `underscores`. Outside quotes, `[...]` is a comment and is
#' ignored. [write_tree()] quotes labels the same way, so trees round-trip
#' exactly.
#'
#' @return A \code{phylo} class object representing the tree.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     infile <- system.file("extdata", "newick.tre", package = "rbiom")
#'     tree <- read_tree(infile)
#'     print(tree)
#'
#'     tree <- read_tree("
#'         (A:0.99,((B:0.87,C:0.89):0.51,(((D:0.16,(E:0.83,F:0.96)
#'         :0.94):0.69,(G:0.92,(H:0.62,I:0.85):0.54):0.23):0.74,J:0.1
#'         2):0.43):0.67);")
#'
#'     # Plotting a tree requires the 'ape' package.
#'     if (requireNamespace("ape", quietly = TRUE))
#'       ape::plot.phylo(tree)
#'
read_tree <- function (src, underscores = FALSE) {

  if (!is_scalar_character(src) || is_na(src))
    cli_abort("`src` must be a file path, URL, or Newick string, not {.type {src}}.")

  # URLs and (possibly compressed) files are read here; literal
  # newick strings are passed straight through.
  is_newick <- grepl("[(;]", src)

  if (is_url(src) || (!is_newick && file.exists(src))) {
    fpath <- as_filepath(src)
    on.exit(fpath$cleanup(), add = TRUE)
    src <- paste(strip_bom(readLines(fpath$path, warn = FALSE, encoding = "UTF-8")), collapse = "")

  } else if (!is_newick) {
    cli_abort("`src` is neither an existing file nor a newick string: {.val {src}}")
  }

  ecodive::read_tree(newick = src, underscores = underscores)
}





#' Create a subtree by specifying tips to keep.
#'
#' Internal nodes left with a single descendant are collapsed (their branch
#' lengths are summed), node labels are kept, and the result is rooted at
#' the most recent common ancestor of the retained tips. Tips are kept in
#' their original order regardless of the order of `tips`.
#'
#' @family phylogeny
#'
#' @inherit documentation_default
#'
#' @param tree  A phylo object, as returned from [read_tree()].
#'
#' @param tips  A character, numeric, or logical vector of tips to keep.
#'
#' @return A \code{phylo} object for the subtree.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     infile <- system.file("extdata", "newick.tre", package = "rbiom")
#'     tree <- read_tree(infile)
#'     tree
#'
#'     subtree <- tree_subset(tree, tips = head(tree$tip.label))
#'     subtree
#'
tree_subset <- function (tree, tips, underscores = FALSE) {

  validate_tree(underscores = underscores)
  params <- eval_envir(environment())


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('tree_subset', params)
  if (!is.null(res <- get_cache_value(cache_file)))
    return (res)


  #________________________________________________________
  # Did we get a tree?
  #________________________________________________________

  if (!inherits(tree, "phylo"))
    cli_abort("Can't convert {.type {tree}} to a `phylo` object.")

  nTips <- length(tree$tip.label)


  #________________________________________________________
  # Tips specified as strings c("SteSpe52", "HrbSpe65", ...)
  #________________________________________________________

  if (is.character(tips)) {

    if (length(x <- setdiff(tips, tree$tip.label)))
      cli_abort("Tips missing from tree: {x}.")

    tips <- which(tree$tip.label %in% tips)
  }


  #________________________________________________________
  # Tips specified as logicals c(TRUE, TRUE, FALSE, TRUE, ...)
  #________________________________________________________

  if (is.logical(tips)) {

    if (length(tips) != nTips)
      cli_abort("Logical vector `tips` must be same length as `tree$tip.label`.")

    tips <- which(tips)
  }


  #________________________________________________________
  # Tips specified as numeric c(5, 1, 2, 31, 3, 12, ...)
  #________________________________________________________

  if (is.numeric(tips)) {
    if (length(tips) > 0) {
      if (anyNA(tips))        cli_abort("`tips` can't contain NA.")
      if (max(tips) > nTips)  cli_abort("`tips` index {.val {max(tips)}} doesn't exist.")
      if (min(tips) < 1)      cli_abort("`tips` index {.val {min(tips)}} doesn't exist.")
      if (any(tips %% 1 > 0)) cli_abort("`tips` indices must be whole numbers.")
    }

  } else {
    cli_abort("`tips` must character, logical, or numeric, not {.type {tips}}.")
  }


  #________________________________________________________
  # Keeping zero or all tips?
  #________________________________________________________

  tips <- sort(unique(as.integer(tips)))

  if (length(tips) < 2)
    cli_abort("At least two tips must be specified.")

  if (length(tips) == nTips)
    return (set_cache_value(cache_file, tree))


  #________________________________________________________
  # Prune the tree
  #________________________________________________________

  res <- prune_tree(tree, tips)

  set_cache_value(cache_file, res)
}



#________________________________________________________
# For each node, the indices of the edges leading to its children
# (an empty integer vector for tips).
#________________________________________________________
tree_children <- function (parent, child, n_node = max(parent, child)) {
  kids <- split(seq_along(parent), factor(parent, levels = seq_len(n_node)))
  unname(kids)
}


#________________________________________________________
# Node ids in preorder (parents before children, siblings in edge
# order), starting from `root`. Iterative; no recursion.
#________________________________________________________
tree_preorder <- function (root, kids, child) {

  n     <- length(kids)
  out   <- integer(n)
  stack <- integer(n)
  top   <- 1L
  k     <- 0L
  stack[1L] <- root

  while (top > 0L) {
    v   <- stack[top]
    top <- top - 1L
    k   <- k + 1L
    out[k] <- v
    cc  <- child[kids[[v]]]
    nc  <- length(cc)
    if (nc > 0L) {
      stack[top + seq_len(nc)] <- rev(cc)
      top <- top + nc
    }
  }

  out[seq_len(k)]
}


#________________________________________________________
# Keep only the tips indexed by `keep` (sorted, unique, >= 2 of them).
#________________________________________________________
prune_tree <- function (tree, keep) {

  edge    <- tree$edge
  parent  <- edge[,1]
  child   <- edge[,2]
  n_tip   <- length(tree$tip.label)
  n_node  <- max(edge)
  el      <- tree$edge.length
  has_len <- !is.null(el)

  root <- setdiff(parent, child)
  if (length(root) != 1) cli_abort("Tree must have exactly one root node.")

  kids    <- tree_children(parent, child, n_node)
  edge_of <- integer(n_node)
  edge_of[child] <- seq_along(child)


  #________________________________________________________
  # Which nodes have a kept descendant, and how many kept children.
  #________________________________________________________
  keep_node <- logical(n_node)
  keep_node[keep] <- TRUE
  n_keep <- integer(n_node)

  for (v in rev(tree_preorder(root, kids, child))) {
    ke <- kids[[v]]
    if (length(ke) > 0L) {
      k <- keep_node[child[ke]]
      n_keep[v]    <- sum(k)
      keep_node[v] <- n_keep[v] > 0L
    }
  }


  #________________________________________________________
  # The new root is the MRCA of the kept tips.
  #________________________________________________________
  new_root <- root
  while (n_keep[new_root] == 1L) {
    cc <- child[kids[[new_root]]]
    new_root <- cc[keep_node[cc]]
  }


  #________________________________________________________
  # Nodes that survive: kept tips and internal nodes with >= 2 kept
  # children. Single-child internal nodes are collapsed into the edge.
  #________________________________________________________
  is_final <- keep_node & (seq_len(n_node) <= n_tip | n_keep >= 2L)

  order      <- tree_preorder(new_root, kids, child)
  new_parent <- integer(n_node)
  acc        <- numeric(n_node)

  for (v in order) {

    if (!keep_node[v]) next
    ke <- kids[[v]]
    if (length(ke) == 0L) next

    cc <- child[ke]
    k  <- keep_node[cc]
    ke <- ke[k]
    cc <- cc[k]

    add <- if (has_len) el[ke] else 0

    if (is_final[v]) {
      new_parent[cc] <- v
      acc[cc]        <- add
    } else {
      new_parent[cc] <- new_parent[v]
      acc[cc]        <- acc[v] + add
    }
  }


  #________________________________________________________
  # Renumber: tips 1..k (original order), then internal nodes in preorder.
  #________________________________________________________
  final_pre <- order[is_final[order]]
  int_final <- final_pre[final_pre > n_tip]

  new_id <- integer(n_node)
  new_id[keep]      <- seq_along(keep)
  new_id[int_final] <- length(keep) + seq_along(int_final)

  cs <- final_pre[final_pre != new_root]

  res <- list(
    edge      = cbind(new_id[new_parent[cs]], new_id[cs]),
    Nnode     = length(int_final),
    tip.label = tree$tip.label[keep] )

  storage.mode(res$edge) <- "integer"

  if (has_len)                     res$edge.length <- unname(acc[cs])
  if (!is.null(tree$node.label))   res$node.label  <- tree$node.label[int_final - n_tip]

  structure(res, class = "phylo", order = "cladewise")
}
