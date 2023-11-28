#' Calculate PCoA and other ordinations, including taxa biplots and statistics.
#' 
#' The biplot parameters (\code{taxa}, \code{unc}, \code{p.top}, and
#' \code{p.adj}) only only have an effect when \code{rank} is not \code{NULL}.
#' 
#' @inherit documentation_cmp
#' @inherit documentation_dist_test
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' @family ordination
#' 
#'        
#' @param ...  Additional arguments to pass on to [uwot::umap()], 
#'        [ape::pcoa()], [vegan::metaMDS()], or [tsne::tsne()].
#'             
#' @return A data.frame with columns \code{.sample}, \code{.weighted}, 
#'         \code{.bdiv}, \code{.ord}, \code{.x}, \code{.y}, and (optionally) 
#'         \code{.z}. Any columns given by \code{md}, \code{split.by}, and 
#'         \code{stat.by} are included as well.\cr
#'         If \code{stat.by} is given, then \code{$stats} and 
#'         \code{$stats$code)} are set.\cr
#'         If \code{rank} is given, then \code{$taxa_coords}, 
#'         \code{$taxa_stats}, and \code{$taxa_stats$code} are set.
#'         
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     ord <- bdiv_ord_table(hmp50, "bray", "pcoa", stat.by="Body Site", rank="g")
#'     head(ord)
#'     
#'     ord$stats
#'     
#'     ord$taxa_stats
#'     
#'     

bdiv_ord_table <- function (
    biom, bdiv = "Bray-Curtis", ord = "UMAP", weighted = TRUE, md = NULL, k = 2, 
    split.by = NULL, stat.by = NULL, tree = NULL, within = NULL, between = NULL,
    test = "adonis2", seed = 0, permutations = 999, rank = -1, taxa = 5, 
    p.top = Inf, p.adj = 'fdr', unc = "singly", ...) {
  
  validate_biom(clone = FALSE)
  validate_tree(null_ok = TRUE)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- eval_envir(environment(), ...)
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  
  
  #________________________________________________________
  # Validate arguments only relevant for rbiom objects.
  #________________________________________________________
  biom %<>% as_percent()
  tree %<>% aa(display = "tree")
  validate_ord(max = Inf)
  validate_bdiv(max = Inf)
  validate_bool("weighted", max = Inf)
  
  validate_meta('md',       null_ok = TRUE, cmp = TRUE, max = Inf)
  validate_meta('stat.by',  null_ok = TRUE, cmp = TRUE)
  validate_meta('split.by', null_ok = TRUE, cmp = TRUE, col_type = "cat")
  
  # Validates and appends to `within` and `between`.
  validate_meta_cmp(c('md', 'stat.by', 'split.by'))
  
  
  
  #________________________________________________________
  # Biplot - determine taxonomic rank
  #________________________________________________________
  if (missing(rank) && is.character(taxa)) {
    rank <- names(which.max(lapply(otu_taxonomy(biom), function (x) sum(x %in% taxa))))
    
  } else {
    validate_rank(max = Inf)
  }
  
  
  
  results <- blply(
    biom   = biom, 
    vars   = split.by, 
    iters  = list(weighted = weighted, bdiv = bdiv),
    prefix = TRUE,
    FUN    = function (b, weighted, bdiv) {
      
      
      dm <- bdiv_distmat(biom = b, bdiv = bdiv, weighted = weighted, tree = tree)
      
      
      #________________________________________________________
      # Sample Ordination(s). x/y/z coordinates.
      #________________________________________________________
      args          <- c(list(dm = dm, ord = ord, k = k), params$.dots)
      sample_coords <- do.call(distmat_ord_table, args)
      attr(dm, 'sample_coords') <- sample_coords
      
      
      #________________________________________________________
      # Statistics: Sample Groups vs Distance Matrix
      #________________________________________________________
      if (!is.null(stat.by)) {
        
        sample_stats <- distmat_stats(
          dm           = dm, 
          groups       = sample_metadata(b, stat.by), 
          test         = test, 
          seed         = seed,
          permutations = permutations )
        
        attr(dm, 'sample_stats') <- sample_stats
      }
      
      
      
      #________________________________________________________
      # Biplot statistics and x/y/z coordinates.
      #________________________________________________________
      if (!is.null(rank)) {
        
        rank_results <- plyr::llply(setNames(rank, rank), function (r) {
          
          taxa_mtx    <- taxa_matrix(biom = b, rank = r, taxa = taxa, unc = unc)
          taxa_stats  <- biplot_taxa_stats(dm, taxa_mtx, test, seed, permutations, p.adj, p.top)
          taxa_coords <- plyr::ddply(sample_coords, ".ord", function (df) {
            biplot_taxa_coords(sample_coords = df, taxa_mtx)
          })
          
          return (list(taxa_stats, taxa_coords))
        })
        
        attr(dm, 'taxa_stats')  <- plyr::ldply(rank_results, .id = ".rank", `[[`, 1)
        attr(dm, 'taxa_coords') <- plyr::ldply(rank_results, .id = ".rank", `[[`, 2)
      }
      
      
      return (dm)
  })
  
  sample_coords <- plyr::ldply(results, attr, 'sample_coords') %>% as_rbiom_tbl()
  sample_stats  <- plyr::ldply(results, attr, 'sample_stats')  %>% as_rbiom_tbl()
  taxa_coords   <- plyr::ldply(results, attr, 'taxa_coords')   %>% as_rbiom_tbl()
  taxa_stats    <- plyr::ldply(results, attr, 'taxa_stats')    %>% as_rbiom_tbl()
  
  if (plyr::empty(sample_coords)) sample_coords <- NULL
  if (plyr::empty(sample_stats))  sample_stats  <- NULL
  if (plyr::empty(taxa_coords))   taxa_coords   <- NULL
  if (plyr::empty(taxa_stats))    taxa_stats    <- NULL
  
  
  if (!is.null(sample_stats)) {
    sample_stats[['.adj.p']] <- signif(p.adjust(sample_stats[['.p.val']], method = p.adj), 3)
    attr(sample_stats, 'tbl_sum') <- c(
      'Test' = paste0(test, " ~ ", coan(stat.by), ". ", permutations, " permutations."))
  }
  
  if (!is.null(taxa_stats)) {
    taxa_stats[['.adj.p']] <- signif(p.adjust(taxa_stats[['.p.val']], method = p.adj), 3)
    if (p.top < 1)
      taxa_stats <- taxa_stats[taxa_stats[['.adj.p']] <= p.top,,drop=FALSE]
    attr(taxa_stats, 'tbl_sum') <- c(
      'Test' = paste0(test, " ~ taxa. ", permutations, " permutations."))
  }
  
  #________________________________________________________
  # Rearrange columns.
  #________________________________________________________
  # sample_coords[['.weighted']] %<>% ifelse("Weighted", "Unweighted")
  # sample_coords %<>% within(.distance <- paste(.weighted, .method))
  # sample_coords %<>% keep_cols(c(".sample", ".distance", ".ord", ".x", ".y", ".z"))
  
  
  for (i in unique(c(split.by, stat.by, md)))
    if (!hasName(sample_coords, i))
      sample_coords[[i]] <- sample_metadata(biom, i)[as.character(sample_coords[['.sample']])]
  
  
  
  
  #________________________________________________________
  # R code for reproducing statistics.
  #________________________________________________________
  
  if (!is.null(sample_stats) || !is.null(taxa_stats)) {
    
    biom <- structure(NA, 'display' = "biom")
    
    
    
    #________________________________________________________
    # Sample metadata groups
    #________________________________________________________
    
    if (!is.null(sample_stats)) {
      
      attr(sample_stats, 'code') <- local({
        
        
        # Simple case: no iterations over the stats section(s).
        #________________________________________________________
        
        if (length(results) == 1)
          return (paste(
            sep = "\n",
            "dm     <- %s" %>% fmt_cmd(bdiv_distmat, biom, bdiv, weighted, tree),
            "groups <- sample_metadata(biom, %s)[attr(dm, 'Labels')]" %>% sprintf(double_quote(stat.by)),
            "set.seed(%i)" %>% sprintf(seed),
            sprintf(permutations, fmt = switch(
              EXPR = test,
              adonis2 = "ptest  <- vegan::adonis2(formula = dm ~ groups, permutations = %i)",
              mrpp    = "ptest  <- vegan::mrpp(dat = dm, grouping = groups, permutations = %i)" )),
            "pstats <- summary(vegan::permustats(ptest))",
            "stats <- with(pstats, signif(data.frame(.stat=statistic, .z=z, .p.val=p, .adj.p=p), 3))",
            "print(stats)" ))
        
        
        # Complex case: stats assembled from multiple iterations.
        #________________________________________________________
        
        iters   <- list(weighted = weighted, bdiv = bdiv)
        dm_args <- list(
          biom   = biom, 
          vars   = split.by, 
          iters  = aa(iters, display = "iters"), 
          prefix = TRUE, 
          FUN    = aa(bdiv_distmat, display = "bdiv_distmat") )
        dm_args[['tree']] <- tree
        
        paste(collapse = "\n", c(
          "iters   <- list(%s)"  %>% sprintf(as.args(iters)),
          "dm_list <- blply(%s)" %>% sprintf(as.args(dm_args, fun = blply)),
          "stats   <- plyr::ldply(dm_list, function (dm) {",
          "  groups <- sample_metadata(biom, %s)[attr(dm, 'Labels')]" %>% sprintf(double_quote(stat.by)),
          "  set.seed(%i)" %>% sprintf(seed),
          sprintf(permutations, fmt = switch(
            EXPR = test,
            adonis2 = "  ptest  <- vegan::adonis2(formula = dm ~ groups, permutations = %i)",
            mrpp    = "  ptest  <- vegan::mrpp(dat = dm, grouping = groups, permutations = %i)" )),
          "  pstats <- summary(vegan::permustats(ptest))",
          "  with(pstats, signif(data.frame(.stat=statistic, .z=z, .p.val=p), 3))",
          "})",
          sprintf("stats[['.adj.p']] <- signif(p.adjust(stats[['.p.val']], %s), 3))", double_quote(p.adj)),
          "print(stats)" ))
        
      }) %>%
        add_class('rbiom_code')
      
    }
    
    
    
    #________________________________________________________
    # Taxa as groups
    #________________________________________________________
    
    if (!is.null(taxa_stats)) {
      
      attr(taxa_stats, 'code') <- local({
        
        
        # Simple case: no iterations over the stats section(s).
        #________________________________________________________
        
        if (length(results) == 1 && length(rank) < 2)
          return (paste(
            sep = "\n",
            "dm         <- %s" %>% fmt_cmd(bdiv_distmat, biom, bdiv, weighted, tree),
            "taxa_mtx   <- t(%s)" %>% fmt_cmd(taxa_matrix, biom, rank, taxa, unc),
            "taxa_stats <- plyr::adply(taxa_mtx, 2L, .id = '.taxa', function (abundances) {",
            "  abundances <- abundances[attr(dm, 'Labels')]",
            "  set.seed(%i)" %>% sprintf(seed),
            sprintf(permutations, fmt = switch(
              EXPR = test,
              adonis2 = "  ptest  <- vegan::adonis2(formula = dm ~ abundances, permutations = %i)",
              mrpp    = "  ptest  <- vegan::mrpp(dat = dm, grouping = abundances, permutations = %i)" )),
            "  pstats <- summary(vegan::permustats(ptest))",
            "  with(pstats, signif(data.frame(.stat=statistic, .z=z, .p.val=p), 3))",
            "})",
            sprintf("taxa_stats[['.adj.p']] <- signif(p.adjust(taxa_stats[['.p.val']], %s), 3)", double_quote(p.adj)),
            "print(taxa_stats)" ))
        
        
        # Complex case: stats assembled from multiple iterations.
        #________________________________________________________
        
        iters   <- list(weighted = weighted, bdiv = bdiv)
        dm_args <- list(
          biom   = biom, 
          vars   = split.by, 
          iters  = aa(iters, display = "iters"), 
          prefix = TRUE, 
          FUN    = aa(bdiv_distmat, display = "bdiv_distmat") )
        dm_args[['tree']] <- tree
        
        paste(collapse = "\n", c(
          "iters      <- list(%s)"  %>% sprintf(as.args(iters)),
          "dm_list    <- blply(%s)" %>% sprintf(as.args(dm_args, fun = blply)),
          "ranks      <- %s"        %>% sprintf(as.args(list(rank))),
          "taxa_stats <- plyr::ldply(dm_list, function (dm) {",
          "  plyr::ldply(setNames(ranks, ranks), .id = '.rank', function (rank) {",
          "    taxa_mtx <- t(%s)" %>% fmt_cmd(taxa_matrix, biom_, rank_, taxa, unc),
          "    plyr::adply(taxa_mtx, 2L, .id = '.taxa', function (abundances) {",
          "      abundances <- abundances[attr(dm, 'Labels')]",
          "      set.seed(%i)" %>% sprintf(seed),
          sprintf(permutations, fmt = switch(
            EXPR = test,
            adonis2 = "      ptest  <- vegan::adonis2(formula = dm ~ abundances, permutations = %i)",
            mrpp    = "      ptest  <- vegan::mrpp(dat = dm, grouping = abundances, permutations = %i)" )),
          "      pstats <- summary(vegan::permustats(ptest))",
          "      with(pstats, signif(data.frame(.stat=statistic, .z=z, .p.val=p), 3))",
          "    })",
          "  })",
          "})",
          sprintf("taxa_stats[['.adj.p']] <- signif(p.adjust(taxa_stats[['.p.val']], %s), 3)", double_quote(p.adj)),
          if (p.top < 1)          { sprintf("taxa_stats <- subset(taxa_stats, .adj.p <= %f)", p.top)
          } else if (p.top < Inf) { sprintf("taxa_stats <- subset(taxa_stats, rank(.p.val) <= %i)", p.top)
          } else                  { NULL },
          "print(taxa_stats)" ))
        
      }) %>%
        add_class('rbiom_code')
    }
    
  }
  
  
  
  attr(sample_coords, 'stats')       <- sample_stats
  attr(sample_coords, 'taxa_coords') <- taxa_coords
  attr(sample_coords, 'taxa_stats')  <- taxa_stats
  
  
  set_cache_value(cache_file, sample_coords)
  return (sample_coords)
}




biplot_taxa_stats <- function (dm, taxa_mtx, test, seed, permutations, p.adj, p.top) {
  
  taxa_stats <- plyr::adply(
    .data        = t(taxa_mtx)[attr(dm, 'Labels'),], 
    .margins     = 2L, 
    .id          = ".taxa", 
    .fun         = distmat_stats,
    dm           = dm, 
    test         = test, 
    seed         = seed,
    permutations = permutations )
  
  if (!is.data.frame(taxa_stats) || plyr::empty(taxa_stats)) return (NULL)
  
  taxa_stats[['.adj.p']] <- p.adjust(taxa_stats[['.p.val']], p.adj)
  
  if (p.top >= 1 && p.top < Inf)
    taxa_stats <- taxa_stats[rank(taxa_stats[['.p.val']]) <= p.top,,drop=FALSE]
  
  if (p.top < 1)
    taxa_stats <- taxa_stats[taxa_stats[['.adj.p']] <= p.top,,drop=FALSE]
  
  
  return (taxa_stats)
}




biplot_taxa_coords <- function (sample_coords, taxa_mtx) {
  
  xyz     <- intersect(names(sample_coords), c('.x', '.y', '.z')) %>% setNames(.,.)
  origins <- lapply(xyz, function (i) { round(mean(sample_coords[[i]]), 10) })
  names(origins) %<>% paste0("0")
  taxa_mtx <- t(taxa_mtx)[as.character(sample_coords[['.sample']]),]
  
  #________________________________________________________
  # Find weighted average x, y, z, and size for each taxon.
  #________________________________________________________
  plyr::adply(taxa_mtx, 2L, .id = ".taxa", function (vals) {
    
    taxa_size <- mean(vals)
    if (taxa_size == 0) return (NULL)
    
    taxa_xyz <- lapply(xyz, function (i) { mean(sample_coords[[i]] * vals) / taxa_size })
    
    data.frame(check.names = FALSE, '.size' = taxa_size, taxa_xyz, origins) %>%
      signif(3)
  })
  
}


