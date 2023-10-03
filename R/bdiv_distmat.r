#' Make a distance matrix of samples vs samples.
#' 
#' @name bdiv_distmat
#' 
#' @family beta_diversity
#' 
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'        object (e.g. as returned from [read_biom()]). For matrices, the 
#'        rows and columns are assumed to be the taxa and samples, respectively.
#'        
#' @param bdiv  The beta diversity distance algorithm to use. Options are:
#'        \code{"Bray-Curtis"}, \code{"Manhattan"}, \code{"Euclidean"}, 
#'        \code{"Jaccard"}, and \code{"UniFrac"}. Non-ambiguous abbreviations 
#'        of the algorithm names are also accepted. For \code{"UniFrac"}, a 
#'        phylogenetic tree must be present in \code{biom} or explicitly 
#'        provided via \code{tree=}. Default: \code{"Bray-Curtis"}.
#'     
#' @param weighted  Take relative abundances into account. When 
#'        \code{weighted=FALSE}, only presence/absence is considered.
#'        Default: \code{TRUE}
#'        
#' @param tree  A \code{phylo} object representing the phylogenetic
#'        relationships of the taxa in \code{biom}. Will be taken from the tree
#'        defined in \code{biom} if not explicitly specified. Only required 
#'        when computing UniFrac distance matrices. Default: \code{NULL}
#'        
#' @return A distance matrix.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- sample_select(hmp50, 1:10)
#'     dm   <- bdiv_distmat(biom, 'unifrac')
#'     as.matrix(dm)[1:4,1:4]
#'     plot(hclust(dm))
#'

bdiv_distmat <- function (biom, bdiv="Bray-Curtis", weighted=TRUE, tree=NULL) {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("bdiv_distmat", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  stopifnot(is_logical(weighted))
  stopifnot(is_character(bdiv))
  stopifnot(is_null(tree) || is(tree, 'phylo'))
  
  
  #________________________________________________________
  # Get the input into a simple_triplet_matrix
  #________________________________________________________
  
  if (is(biom, "simple_triplet_matrix")) { counts <- biom
  } else if (is(biom, "BIOM"))           { counts <- biom$counts
  } else if (is(biom, "matrix"))         { counts <- slam::as.simple_triplet_matrix(biom)
  } else {
    stop("biom must be a matrix, simple_triplet_matrix, or BIOM object.")
  }
  
  
  #________________________________________________________
  # Sanity check the matrix
  #________________________________________________________
  if (!is.numeric(counts[['v']]))     stop("The abundance matrix must be numeric.")
  if (!all(is.finite(counts[['v']]))) stop("Non-finite values in abundance matrix.")
  
  
  #________________________________________________________
  # Enable abbreviations of metric names.
  #________________________________________________________
  bdiv %<>% validate_arg(biom, 'bdiv', tree=tree, n = 1L)
  
  
  #________________________________________________________
  # Find the UniFrac tree
  #________________________________________________________
  
  if (bdiv == "UniFrac") {
    
    # Find a tree for the UniFrac algorithm
    if (!is(tree, "phylo")) {
      if (is(biom, "BIOM")) {
        if (is(biom$phylogeny, "phylo")) {
          tree <- biom$phylogeny
        }
      }
      if (is(tree, "character")) {
        if (file.exists(tree)) {
          tree <- read_tree(tree)
        }
      }
      if (!is(tree, "phylo")) {
        stop("No tree provided to bdiv_distmat().")
      }
    }
    
    
    # Make sure the matrix has rownames set
    if (is_null(rownames(counts)))
      stop("The abundance matrix does not have rownames set to Taxa IDs.")
    
    
    # Abundance matrix's Taxa IDs aren't all in the tree
    if (length(missing <- setdiff(rownames(counts), tree$tip.label)) > 0) {
      
      # Try swapping spaces/underscores in tip labels
      if (any(grepl("_", missing, fixed = TRUE))) {
        tree$tip.label %<>% gsub(" ", "_", ., fixed = TRUE)
      } else if (any(grepl(" ", missing, fixed = TRUE))) {
        tree$tip.label %<>% gsub("_", " ", ., fixed = TRUE)
      }
      
      if (!all(rownames(counts) %in% tree$tip.label))missing %>%
        glue_collapse(sep = ", ", width = 30, last = " and ") %>%
        paste("OTUs missing from reference tree:", .) %>%
        stop()
    }
    
    # Prune the tree down to only what we need
    if (length(setdiff(tree$tip.label, rownames(counts))) > 0)
      tree <- tree_subset(tree, rownames(counts))
    
    counts <- counts[as.character(tree$tip.label),]
  }
  
  
  #________________________________________________________
  # Order the sparse matrix's values by sample, then by taxa
  #________________________________________________________
  
  ord      <- order(counts$j, counts$i)
  counts$i <- counts$i[ord]
  counts$j <- counts$j[ord]
  counts$v <- counts$v[ord]
  
  
  #________________________________________________________
  # Run C++ implemented dissimilarity algorithms multithreaded
  #________________________________________________________
  
  if (bdiv == "UniFrac") {
    
    dm <- par_unifrac(counts, tree, ifelse(weighted, 1L, 0L))
    
    
  } else {
    
    counts <- t(as.matrix(counts))
    dm <- par_beta_div(counts, tolower(bdiv), ifelse(weighted, 1L, 0L))
    dm <- as.dist(dm)
    attr(dm, 'Labels') <- rownames(counts)
    
  }
  
  
  #________________________________________________________
  # Return the dist object
  #________________________________________________________
  set_cache_value(cache_file, dm)
  return (dm)
}




#' Make a data.frame of distances between samples.
#' 
#' @inherit bdiv_distmat params
#' 
#' @family beta_diversity
#'        
#' @param md  Include metadata in the output data frame? Options are: 
#'        \itemize{
#'          \item{\code{FALSE} - }{ Don't include metadata. }
#'          \item{\code{TRUE} - }{ Include all metadata. }
#'          \item{\emph{character vector} - }{ Include only the specified metadata 
#'          columns. Column names can be prefixed with \bold{==} or \bold{!=} to 
#'          indicate that only within or between groupings, respectively, are to 
#'          be kept. See examples below. }
#'        }
#'        Default: \code{FALSE}
#'        
#' @param within,between   Metadata column name(s) for intra- or inter- sample 
#'        comparisons. Default: \code{within=NULL, between=NULL}
#'        
#' @return A data.frame with first three columns named ".sample1", ".sample2", 
#'         and ".distance".
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     # Return in long format with metadata
#'     biom <- sample_select(hmp50, 18:21)
#'     bdiv_table(biom, 'unifrac', md = c("Body Site", "Sex"))
#'     
#'     # Only look at distances among the stool samples
#'     bdiv_table(biom, 'unifrac', md = c("==Body Site", "Sex"))
#'     
#'     # Or between males and females
#'     bdiv_table(biom, 'unifrac', md = c("Body Site", "!=Sex"))
#'

bdiv_table <- function (
    biom, bdiv="Bray-Curtis", weighted=TRUE, tree=NULL, 
    md=FALSE, within=NULL, between=NULL) {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("bdiv_table", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  history <- paste0(collapse = "\n", c(
    attr(biom, 'history', exact = TRUE),
    sprintf("data <- bdiv_table(%s)", as.args(params, fun = bdiv_table)) ))
  
  
  #________________________________________________________
  # Remove "!=" and "==" prefixes.
  #________________________________________________________
  if (isFALSE(md)) md <- NULL
  if (isTRUE(md))  md <- names(sample_metadata(biom))
  
  within  <- unique(sub("^[!=]=", "", c(within,  md[startsWith(md, "==")])))
  between <- unique(sub("^[!=]=", "", c(between, md[startsWith(md, "!=")])))
  md      <- unique(sub("^[!=]=", "", c(md, within, between)))
  
  if (any(within %in% between))
    stop(
      "Metadata name '", 
      paste0(intersect(within, between), collapse = ", "), 
      "' cannot be set as both a within (==) and between (!=) grouping.")
  
  
  #________________________________________________________
  # Compute the distance matrix
  #________________________________________________________
  dm <- bdiv_distmat(biom = biom, bdiv = bdiv, weighted = weighted, tree = tree)
  
  
  #________________________________________________________
  # Convert to long form
  #________________________________________________________
  df <- as.matrix(dm)
  df <- data.frame(
    stringsAsFactors = FALSE,
    '.sample1'       = rownames(df)[row(df)],
    '.sample2'       = colnames(df)[col(df)],
    '.distance'      = as.numeric(df)
  )
  df <- subset(df, .sample1 < .sample2)
  
  
  #________________________________________________________
  # Add metadata columns
  #________________________________________________________
  for (col in md) {
    
    map <- sample_metadata(biom, col)
    v1  <- as.character(map[df$.sample1])
    v2  <- as.character(map[df$.sample2])
    
    
    # Limit to only within or between comparisons.
    #________________________________________________________
    if (col %in% within)  df <- df[v1 == v2,,drop=FALSE]
    if (col %in% between) df <- df[v1 != v2,,drop=FALSE]
    v1 <- as.character(map[df$.sample1])
    v2 <- as.character(map[df$.sample2])
    
    
    # Change "Male vs Female" to "Female vs Male" (alphabetical).
    #________________________________________________________
    df[[col]] <- paste(
      ifelse(v1 < v2, v1, v2), 
      "vs", 
      ifelse(v1 < v2, v2, v1) )
    
    
    # Change "Male vs Male" to "Male".
    #________________________________________________________
    df[[col]] <- ifelse(v1 == v2, v1, df[[col]])
    
    
    # Keep factors as factors when possible
    #________________________________________________________
    if (is.factor(map)) {
      if (col %in% within) {
        df[[col]] <- factor(df[[col]], levels = levels(map))
      } else {
        df[[col]] <- factor(df[[col]])
      }
    }
    
  }
  
  
  attr(df, 'response') <- ".distance"
  attr(df, 'history')  <- history
  
  
  set_cache_value(cache_file, df)
  return (df)
}





#' Run statistics on a distance matrix vs a categorical or numeric variable.
#' 
#' @name distmat_stats
#' 
#' @family beta_diversity
#' 
#' @param dm  The distance matrix as a \code{dist} object, or a list of
#'        \code{dist} objects.
#'        
#' @param groups  A named vector of grouping values. The names should 
#'        correspond to \code{attr(dm, 'Labels')}.
#'        
#' @param test  The name of the function from the \code{vegan} R package to 
#'        use. Options are: \code{"adonis2"} and \code{"mrpp"}. Non-ambiguous 
#'        abbreviations of the function names are also accepted. 
#'        Default: \code{"adonis2"}
#'        
#' @param permutations  Number of random permutations to use. 
#'        Default: \code{999}
#'        
#' @param seed  Random seed for permutations. Default: \code{0}
#'
#' @param p.adj   Method to use for multiple comparisons adjustment of 
#'        p-values. Run \code{p.adjust.methods} for a list of available 
#'        options. Default: \code{"fdr"}.
#'        
#' @return A data.frame with summary statistics from [vegan::permustats()]. The columns are:
#'        \itemize{
#'          \item{\emph{stat} - }{ The observed statistic. For mrpp, this is the overall weighted mean of group mean distances. }
#'          \item{\emph{z} - }{ The difference of observed statistic and mean of permutations divided by the standard deviation of permutations (also known as z-values). Evaluated from permuted values without observed statistic. }
#'          \item{\emph{p.val} - }{ Include only the specified metadata columns. }
#'          \item{\emph{adj.p} - }{ P-value after correcting for the false discovery rate with the \code{p.adj} method. }
#'          \item{\emph{test} - }{ The \code{test} and \code{p.adj} methods used. }
#'        }
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- sample_select(hmp50, 1:10)
#'     dm   <- bdiv_distmat(biom, 'unifrac')
#'     distmat_stats(dm, groups = sample_metadata(biom, 'Body Site'))
#'     distmat_stats(dm, groups = sample_metadata(biom, 'Age'))
#'

distmat_stats <- function (dm, groups, test = "adonis2", seed = 0, permutations = 999, p.adj = "fdr") {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("distmat_stats", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Sanity checking
  #________________________________________________________
  test  <- match.arg(test, c("adonis2", "mrpp"))
  p.adj <- match.arg(p.adj, p.adjust.methods)
  stopifnot(is_scalar_integerish(seed) && !is.na(seed))
  stopifnot(is_scalar_integerish(permutations) && !is.na(permutations))
  stopifnot(!is.null(names(groups)))
  stopifnot(is(dm, 'dist') || is.list(dm))
  
  
  if (is(dm, 'dist')) {
    
    
    #________________________________________________________
    # Just a single distance matrix.
    #________________________________________________________
    
    grouping <- groups[attr(dm, 'Labels')]
    set.seed(seed)
    
    pstats <- tryCatch(
      expr  = {
        ptest <- switch(
          EXPR = test,
          adonis2 = vegan::adonis2(formula = dm ~ grouping, permutations = permutations),
          mrpp    = vegan::mrpp(dat = dm, grouping = grouping, permutations = permutations) )
        summary(vegan::permustats(ptest)) },
      error = function (e) list(statistic=NA, z=NA, p=NA) )
    
    results <- with(pstats, signif(data.frame(.stat=statistic, .z=z, .p.val=p), 3))
    results[['.adj.p']] <- p.adjust(results[['.p.val']], method = p.adj)
    results[['.test']]  <- paste(sep = "; ", test, p.adj)
    
    rownames(results) <- NULL
    
    
    attr(results, 'cmd') <- c(
      sprintf("set.seed(%i)", seed),
      "grouping <- groups[attr(dm, 'Labels')]",
      sprintf(permutations, fmt = switch(
        EXPR = test,
        adonis2 = "ptest    <- vegan::adonis2(formula = dm ~ grouping, permutations = %i)",
        mrpp    = "ptest    <- vegan::mrpp(dat = dm, grouping = grouping, permutations = %i)" )),
      "pstats   <- tryCatch(",
      "  expr  = summary(vegan::permustats(ptest)),",
      "  error = function (e) list(statistic=NA, z=NA, p=NA) )",
      "with(pstats, data.frame(.stat=statistic, .z=z, .p.val=p))",
      sprintf("results[['.adj.p']] <- p.adjust(results[['.p.val']], method = '%s')", p.adj),
      sprintf("results[['.test']]  <- %s)", double_quote(paste(sep = "; ", test, p.adj))) )
    
    
    
  } else {
    
    #________________________________________________________
    # Multiple distance matrices in a list().
    #________________________________________________________
    
    
    results <- plyr::ldply(dm, function (dm2) {
      
      grouping <- groups[attr(dm2, 'Labels')]
      set.seed(seed)
      
      pstats <- tryCatch(
        expr  = {
          ptest <- switch(
            EXPR = test,
            adonis2 = vegan::adonis2(formula = dm2 ~ grouping, permutations = permutations),
            mrpp    = vegan::mrpp(dat = dm2, grouping = grouping, permutations = permutations) )
          summary(vegan::permustats(ptest)) },
        error = function (e) list(statistic=NA, z=NA, p=NA) )
      
      with(pstats, signif(data.frame(.stat=statistic, .z=z, .p.val=p), 3))
    })
    results[['.adj.p']] <- p.adjust(results[['.p.val']], method = p.adj)
    results[['.test']]  <- paste(sep = "; ", test, p.adj)
    
    
    attr(results, 'cmd') <- c(
      "results <- plyr::ldply(dm, function (dm2) {",
      sprintf("  set.seed(%i)", seed),
      "  grouping <- groups[attr(dm2, 'Labels')]",
      sprintf(permutations, fmt = switch(
        EXPR = test,
        adonis2 = "  ptest    <- vegan::adonis2(formula = dm2 ~ grouping, permutations = %i)",
        mrpp    = "  ptest    <- vegan::mrpp(dat = dm2, grouping = grouping, permutations = %i)" )),
      "  pstats   <- tryCatch(",
      "    expr  = summary(vegan::permustats(ptest)),",
      "    error = function (e) list(statistic=NA, z=NA, p=NA) )",
      "  with(pstats, data.frame(.stat=statistic, .z=z, .p.val=p))",
      "}",
      sprintf("results[['.adj.p']] <- p.adjust(results[['.p.val']], method = '%s')", p.adj),
      sprintf("results[['.test']]  <- %s)", double_quote(paste(sep = "; ", test, p.adj))) )
  }
  
  
  #________________________________________________________
  # Return the statistics table.
  #________________________________________________________
  set_cache_value(cache_file, results)
  return (results)
  
}







#' INTERNAL simple case ordination handler. Use bdiv_ord_table() instead.
#' 
#' @name ordinate
#' @noRd
#' 
#' @param dm   A \code{dist}-class distance matrix, as returned from 
#'        [bdiv_distmat()] or [stats::dist()]. Required.
#'        
#' @param ord    Method for reducing dimensionality. Options are:
#'        \itemize{
#'            \item{\code{"UMAP"} - }{ Uniform manifold approximation and projection; [uwot::umap()]. }
#'            \item{\code{"PCoA"} - }{ Principal coordinate analysis; [ape::pcoa()]. }
#'            \item{\code{"NMDS"} - }{ Nonmetric multidimensional scaling; [vegan::metaMDS()]. }
#'            \item{\code{"tSNE"} - }{ t-distributed stochastic neighbor embedding; [tsne::tsne()]. }
#'        }
#'        Default: \code{"UMAP"} \cr\cr
#'        Non-ambiguous abbreviations are also accepted.
#'        
#' @param k   Number of ordination dimensions to return. Either \code{2L} or 
#'        \code{3L}. Default: \code{2L}
#'        
#' @param ...  Additional arguments to pass on to [uwot::umap()], 
#'        [ape::pcoa()], [vegan::metaMDS()], or [tsne::tsne()].
#'        
#' @return A data.frame with columns \code{.sample}, \code{.x}, \code{.y}, 
#'         and (optionally) \code{.z}.
#' 
#' @examples
#'     library(rbiom)
#'     
#'     dm  <- bdiv_distmat(hmp50, "bray")
#'     ord <- rbiom:::distmat_ordinate(dm, "PCoA")
#'     head(ord)
#'     

distmat_ordinate <- function (dm, ord = "PCoA", k = 2L, ...) {
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  dots <- list(...)
  ord  <- match.arg(ord, c("PCoA", "tSNE", "NMDS", "UMAP"))
  stopifnot(is_integerish(k) && k >= 2 && k <= 3)
  stopifnot(is(dm, 'dist'))
  
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("ordinate", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  
  #________________________________________________________
  # Call 3rd party function with proper arguments.
  #________________________________________________________
  if (ord == "PCoA") {
    args <- c(dots, list(D = dm))
    res  <- do.call(ape::pcoa, args)[['vectors']][,1:k]
    
  } else if (ord == "tSNE") {
    args <- c(dots, list(X = dm, k = k))
    res  <- suppressMessages(do.call(tsne::tsne, args))
    
  } else if (ord == "NMDS") {
    args <- c(dots, list(comm = dm, k = k))
    args[['trace']] %<>% if.null(0)
    res  <- do.call(vegan::metaMDS, args)[['points']]
    rownames(res) <- attr(dm, "Labels", exact = TRUE)
    
  } else if (ord == "UMAP") {
    args <- c(dots, list(X = dm))
    args[['n_components']] %<>% if.null(k)
    args[['n_neighbors']]  %<>% if.null(max(2, min(100, as.integer(attr(dm, 'Size') / 3))))
    res  <- do.call(uwot::umap, args)
  }
  
  res <- signif(as.data.frame(res), 3)
  colnames(res) <- head(c(".x", ".y", ".z"), ncol(res))
  res[[".sample"]] <- rownames(res)
  rownames(res) <- NULL
  res %<>% keep_cols(".sample", ".x", ".y", ".z")
  
  
  set_cache_value(cache_file, res)
  return (res)
}




