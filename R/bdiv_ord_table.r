#' Calculate PCoA and other ordinations, including taxa biplots and statistics.
#' 
#' The biplot parameters (\code{taxa}, \code{unc}, \code{p.top}, and
#' \code{p.adj}) only only have an effect when \code{rank} is not `NULL`.
#' 
#' @inherit documentation_dist_test
#' @inherit documentation_rank.NULL
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
#'         \code{stat.by} are included as well.
#'         If \code{stat.by} is given, then \code{$stats} and 
#'         \code{$stats$code)} are set.
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
    biom, bdiv = "Bray-Curtis", ord = "PCoA", weighted = TRUE, md = NULL, k = 2, 
    stat.by = NULL, split.by = NULL, tree = NULL, 
    test = "adonis2", seed = 0, permutations = 999, rank = NULL, taxa = 6, 
    p.top = Inf, p.adj = 'fdr', unc = "singly", underscores = FALSE, ...) {
  
  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE, underscores = underscores)
  

  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  
  biom <- biom$clone()
  biom$counts %<>% rescale_cols()
  tree %<>% aa(display = "tree")
  
  validate_ord(max = Inf)
  validate_bdiv(max = Inf)
  validate_bool("weighted", max = Inf)
  
  validate_biom_field('md',       null_ok = TRUE, max = Inf)
  validate_biom_field('stat.by',  null_ok = TRUE)
  validate_biom_field('split.by', null_ok = TRUE, col_type = "cat")
  
  na.omit(biom, clone = FALSE, fields = c(stat.by, split.by))
  
  
  #________________________________________________________
  # Biplot - determine taxonomic rank
  #________________________________________________________
  if (is.character(taxa)) {
    rank <- names(which.max(lapply(biom$taxonomy, function (x) sum(x %in% taxa))))
    
  } else {
    validate_rank(max = Inf, null_ok = TRUE)
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
      sample_coords <- distmat_ord_table(dm = dm, ord = ord, k = k, seed = seed, ...)
      attr(dm, 'sample_coords') <- sample_coords
      
      
      #________________________________________________________
      # Statistics: Sample Groups vs Distance Matrix
      #________________________________________________________
      if (!is.null(stat.by)) {
        
        sample_stats <- distmat_stats(
          dm           = dm, 
          groups       = pull(b, stat.by), 
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
    sample_stats[['.adj.p']] <- p.adjust(sample_stats[['.p.val']], method = p.adj)
    attr(sample_stats, 'tbl_sum') <- c(
      'Test' = paste0(test, " ~ ", coan(stat.by), ". ", permutations, " permutations."))
  }
  
  if (!is.null(taxa_stats)) {
    taxa_stats[['.adj.p']] <- p.adjust(taxa_stats[['.p.val']], method = p.adj)
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
      sample_coords[[i]] <- pull(biom, i)[as.character(sample_coords[['.sample']])]
  
  
  
  
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
            "groups <- pull(biom, %s)[attr(dm, 'Labels')]" %>% sprintf(double_quote(stat.by)),
            "set.seed(%i)" %>% sprintf(seed),
            sprintf(permutations, fmt = switch(
              EXPR = test,
              adonis2 = "ptest  <- vegan::adonis2(formula = dm ~ groups, permutations = %i)",
              mrpp    = "ptest  <- vegan::mrpp(dat = dm, grouping = groups, permutations = %i)" )),
            "pstats <- summary(vegan::permustats(ptest))",
            "stats  <- with(pstats, data.frame(statistic, z, p))" ))
        
        
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
          "  groups <- pull(biom, %s)[attr(dm, 'Labels')]" %>% sprintf(double_quote(stat.by)),
          "  set.seed(%i)" %>% sprintf(seed),
          sprintf(permutations, fmt = switch(
            EXPR = test,
            adonis2 = "  ptest  <- vegan::adonis2(formula = dm ~ groups, permutations = %i)",
            mrpp    = "  ptest  <- vegan::mrpp(dat = dm, grouping = groups, permutations = %i)" )),
          "  pstats <- summary(vegan::permustats(ptest))",
          "  with(pstats, data.frame(statistic, z, p))",
          "})" ))
        
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
            "  with(pstats, data.frame(statistic, z, p))",
            "})" ))
        
        
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
        
        tm_args <- list(
          biom = aa(NA, display = 'biom'), 
          rank = aa(NA, display = 'rank'), 
          taxa = taxa, 
          unc  = unc )
        
        paste(collapse = "\n", c(
          "iters      <- list(%s)"  %>% sprintf(as.args(iters)),
          "dm_list    <- blply(%s)" %>% sprintf(as.args(dm_args, fun = blply)),
          "ranks      <- %s"        %>% sprintf(as.args(list(rank))),
          "taxa_stats <- plyr::ldply(dm_list, function (dm) {",
          "  plyr::ldply(setNames(ranks, ranks), .id = '.rank', function (rank) {",
          "    taxa_mtx <- t(taxa_matrix(%s))" %>% sprintf(as.args(tm_args, fun = taxa_matrix)),
          "    plyr::adply(taxa_mtx, 2L, .id = '.taxa', function (abundances) {",
          "      abundances <- abundances[attr(dm, 'Labels')]",
          "      set.seed(%i)" %>% sprintf(seed),
          sprintf(permutations, fmt = switch(
            EXPR = test,
            adonis2 = "      ptest  <- vegan::adonis2(formula = dm ~ abundances, permutations = %i)",
            mrpp    = "      ptest  <- vegan::mrpp(dat = dm, grouping = abundances, permutations = %i)" )),
          "      pstats <- summary(vegan::permustats(ptest))",
          "      with(pstats, data.frame(statistic, z, p))",
          "    })",
          "  })",
          "})",
          sprintf("taxa_stats[['adj.p']] <- p.adjust(taxa_stats[['p']], '%s')", p.adj),
          if (p.top < 1)          { sprintf("taxa_stats <- subset(taxa_stats, adj.p <= %f)", p.top)
          } else if (p.top < Inf) { sprintf("taxa_stats <- subset(taxa_stats, rank(p) <= %i)", p.top)
          } else                  { NULL } ))
        
      }) %>%
        add_class('rbiom_code')
    }
    
  }
  
  
  
  attr(sample_coords, 'stats')       <- sample_stats
  attr(sample_coords, 'taxa_coords') <- taxa_coords
  attr(sample_coords, 'taxa_stats')  <- taxa_stats
  
  
  return (sample_coords)
}




biplot_taxa_stats <- function (dm, taxa_mtx, test, seed, permutations, p.adj, p.top) {
  
  if (test == "none") return (NULL)
  
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


