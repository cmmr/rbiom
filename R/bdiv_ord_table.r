#' Calculate PCoA and other ordinations, including taxa biplots and statistics.
#' 
#' @name bdiv_ord_table
#' 
#' @family ordination
#' @family beta_diversity
#' 
#' @inherit bdiv_distmat params
#' 
#' @param biom   A \code{BIOM} object, as returned from [read_biom()].
#'        Alternatively, a {dist} class distance matrix can be given, in which
#'        case only the parameters \code{ord}, \code{k}, and \code{...} are
#'        allowed.
#'        
#' @param bdiv   The beta diversity distance algorithm to use. Options are:
#'        \code{"Bray-Curtis"}, \code{"Manhattan"}, \code{"Euclidean"}, 
#'        \code{"Jaccard"}, and \code{"UniFrac"}. A phylogenetic tree must be
#'        present in \code{biom} or explicitly provided via \code{tree=} to use
#'        the UniFrac methods.
#'        Multiple values allowed. Default: \code{"Bray-Curtis"}.
#'        
#' @param ord    Method for reducing dimensionality. Options are:
#'        \itemize{
#'            \item{\code{"UMAP"} - }{ Uniform manifold approximation and projection; [uwot::umap()]. }
#'            \item{\code{"PCoA"} - }{ Principal coordinate analysis; [ape::pcoa()]. }
#'            \item{\code{"NMDS"} - }{ Nonmetric multidimensional scaling; [vegan::metaMDS()]. }
#'            \item{\code{"tSNE"} - }{ t-distributed stochastic neighbor embedding; [tsne::tsne()]. }
#'        }
#'        Default: \code{"UMAP"} \cr\cr
#'        Multiple values allowed. Non-ambiguous abbreviations are allowed.
#'        
#' @param weighted   Take relative abundances into account. When 
#'        \code{weighted=FALSE}, only presence/absence is considered.
#'         Multiple values allowed. Default: \code{TRUE}
#'     
#' @param md   Include metadata in the output data frame? Options are: 
#'        \itemize{
#'          \item{\code{NULL} - }{ Don't include metadata. (Default) }
#'          \item{\code{TRUE} - }{ Include all metadata. }
#'          \item{\emph{character vector} - }{ Include only the specified metadata columns. }
#'        }
#'        Default: \code{NULL}
#'        
#' @param k   Number of ordination dimensions to return. Either \code{2L} or 
#'        \code{3L}. Default: \code{2L}
#'        
#' @param split.by   Name(s) of metadata columns that the data should be split
#'        by prior to calculating distance matrices, ordinations, and statistics.
#'        Equivalent to the \code{facet.by} parameter in [bdiv_ord_plot()].
#'        Default: \code{NULL}
#'        
#' @param stat.by   Name of the metadata column for statistical groups. 
#'        Equivalent to the \code{color.by} parameter in [bdiv_ord_plot()].
#'        Default: \code{NULL}
#'        
#' @param test,seed,permutations   Passthrough parameters to [distmat_stats()] 
#'        for assessing significance of \code{stat.by} and \code{rank} groups. 
#'        Default: \code{test="adonis2", seed=0, permutations=999}
#'        
#' @param rank   [Biplot] What rank of taxa to use for the biplot (e.g. 
#'        "Phylum"), or \code{NULL} for no biplot. Run \code{taxa_ranks()} to 
#'        see all options for a given BIOM object. Default: \code{NULL}
#'        
#' @param taxa   [Biplot] Which taxa to include. An integer value will return 
#'        the top n most abundant taxa. A value 0 <= n < 1 will return any taxa 
#'        with that mean abundance or greater (e.g. 0.1). A character vector of
#'        taxon names will show only those taxa. Default: \code{5}
#'        
#' @param unc   [Biplot] How to handle unclassified, uncultured, and similarly 
#'        ambiguous taxa names. Options are: 
#'        \itemize{
#'          \item{\code{"singly"} - }{ Replaces them with the OTU name. }
#'          \item{\code{"grouped"} - }{ Replaces them with a higher rank's name. }
#'          \item{\code{"drop"} - }{ Excludes them from the result. }
#'          \item{\code{"asis"} - }{ To not check/modify any taxa names. }
#'        }
#'        Default: \code{"singly"} \cr\cr
#'        Non-ambiguous abbreviations are allowed.
#'
#' @param p.top   [Biplot] Only return taxa with the most significant 
#'        differences in abundance. If \code{p.top} is >= 1, then the 
#'        \code{p.top} most significant taxa are displayed. If \code{p.top} 
#'        is less than one, all taxa with an adjusted p-value <= \code{p.top} 
#'        are displayed. Recommended to be used in combination with the 
#'        \code{taxa} parameter to set a lower bound on the mean abundance of 
#'        considered taxa. Default: \code{Inf}
#'
#' @param p.adj   [Biplot] Method to use for multiple comparisons adjustment of 
#'        p-values. Run \code{p.adjust.methods} for a list of available options.
#'        Default: \code{fdr}
#'        
#' @param ...  Additional arguments to pass on to [uwot::umap()], 
#'        [ape::pcoa()], [vegan::metaMDS()], or [tsne::tsne()].
#'             
#' @return A data.frame with columns \code{.sample}, \code{.ord}, \code{.x}, 
#'         \code{.y}, and (optionally) \code{.z}.
#'         
#'         If \code{biom} is a \code{BIOM} object, then \code{.weighted}, 
#'         \code{.bdiv}, and any columns given by \code{md}, \code{split.by}, 
#'         and \code{stat.by} are included as well.
#'         
#'         If \code{stat.by} is given, then \code{attr(, 'sample_stats')}
#'         and \code{attr(, 'sample_stats_cmds')} are set.
#'         
#'         If \code{rank} is given, then \code{attr(, 'taxa_coords')},
#'         \code{attr(, 'taxa_stats')}, and \code{attr(, 'taxa_stats_cmds')} 
#'         are set.
#'         
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     ord <- bdiv_ord_table(hmp50, "bray", "pcoa", stat.by="Body Site", rank="g")
#'     head(ord)
#'     attr(ord, 'sample_stats')
#'     attr(ord, 'taxa_stats')
#'     
#'     

bdiv_ord_table <- function (
    biom, bdiv="Bray-Curtis", ord="UMAP", weighted=TRUE, md=NULL, k=2, 
    split.by=NULL, stat.by=NULL, tree=NULL, 
    test="adonis2", seed=0, permutations=999, rank=NULL, taxa=5, 
    p.top=Inf, p.adj='fdr', unc="singly", ...) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(c(as.list(environment()), list(...)), eval)
  cache_file <- get_cache_file("bdiv_ord_table", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  stopifnot(is(biom, 'dist') || is(biom, 'BIOM'))
  dots <- list(...)
  ord %<>% validate_arg(NULL, 'ord', n = c(1,Inf))
  
  
  
  #________________________________________________________
  # Can't do as much when biom is a dist object
  #________________________________________________________
  if (is(biom, 'dist')) {
    
    extraneous_args <- setdiff(c('ord', 'k', '...'), formalArgs(f))
    invalid_args    <- intersect(extraneous_args, names(match.call()))
    if (length(invalid_args) > 0) {
      invalid_args %<>% paste(collapse = ", ")
      stop("In bdiv_ord_table() `biom` must be a BIOM object to use: ", invalid_args)
    }
    
    results <- plyr::ldply(setNames(ord, ord), .id = ".ord", function (o) {
      do.call(distmat_ordinate, c(list(dm = biom, ord = o, k = k), dots))
    })
    
    set_cache_value(cache_file, results)
    return (results)
  }
  
  
  
  #________________________________________________________
  # Validate arguments only relevant for BIOM objects
  #________________________________________________________
  biom     %<>% as_percent()
  tree     %<>% aa(display = deparse(substitute(tree)))
  bdiv     %<>% validate_arg(biom, 'bdiv',     'bdiv', n = c(1,Inf))
  weighted %<>% validate_arg(NULL, 'weighted',         n = c(1,Inf))
  
  if (isTRUE(md)) md <- colnames(sample_metadata(biom))
  md       %<>% validate_arg(biom, 'md',       'meta')
  split.by %<>% validate_arg(biom, 'split.by', 'meta')
  stat.by  %<>% validate_arg(biom, 'stat.by',  'meta')
  
  
  
  #________________________________________________________
  # Biplot - determine taxonomic rank
  #________________________________________________________
  if (!is_null(rank)) {
    rank %<>% validate_arg(biom, 'rank', n = c(1,10))
  } else if (is.character(taxa)) {
    rank <- names(which.max(apply(otu_taxonomy(biom), 2L, function (x) sum(x %in% taxa))))
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
      sample_coords <- plyr::ldply(setNames(ord, ord), .id = ".ord", function (o) {
        do.call(distmat_ordinate, c(list(dm = dm, ord = o, k = k), dots))
      })
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
          permutations = permutations,
          p.adj        = p.adj )
        
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
  
  sample_coords <- plyr::ldply(results, attr, 'sample_coords')
  sample_stats  <- plyr::ldply(results, attr, 'sample_stats')
  taxa_coords   <- plyr::ldply(results, attr, 'taxa_coords')
  taxa_stats    <- plyr::ldply(results, attr, 'taxa_stats')
  
  if (nrow(sample_coords) == 0) sample_coords <- NULL
  if (nrow(sample_stats)  == 0) sample_stats  <- NULL
  if (nrow(taxa_coords)   == 0) taxa_coords   <- NULL
  if (nrow(taxa_stats)    == 0) taxa_stats    <- NULL
  
  
  if (!is.null(sample_stats)) {
    sample_stats[['.adj.p']] <- signif(p.adjust(sample_stats[['.p.val']], method = p.adj), 3)
  }
  
  if (!is.null(taxa_stats)) {
    taxa_stats[['.adj.p']] <- signif(p.adjust(taxa_stats[['.p.val']], method = p.adj), 3)
    if (p.top < 1)
      taxa_stats <- taxa_stats[taxa_stats[['.adj.p']] <= p.top,,drop=FALSE]
  }
  
  #________________________________________________________
  # Rearrange columns.
  #________________________________________________________
  # sample_coords[['.weighted']] %<>% ifelse("Weighted", "Unweighted")
  # sample_coords %<>% within(.distance <- paste(.weighted, .method))
  # sample_coords %<>% keep_cols(c(".sample", ".distance", ".ord", ".x", ".y", ".z"))
  
  for (i in unique(c(split.by, stat.by, md)))
    sample_coords[[i]] <- sample_metadata(biom, i)[sample_coords[['.sample']]]
  
  
  attr(sample_coords, 'sample_stats') <- sample_stats
  attr(sample_coords, 'taxa_coords')  <- taxa_coords
  attr(sample_coords, 'taxa_stats')   <- taxa_stats
  
  
  
  #________________________________________________________
  # R code for reproducing statistics.
  #________________________________________________________
  
  if (!is.null(sample_stats) || !is.null(taxa_stats)) {
    
    biom <- structure(NA, 'display' = "biom")
    
    
    
    #________________________________________________________
    # Sample metadata groups
    #________________________________________________________
    
    if (!is.null(sample_stats)) {
      
      attr(sample_coords, 'sample_stats_cmds') <- local({
        
        
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
            "sample_stats <- with(pstats, signif(data.frame(.stat=statistic, .z=z, .p.val=p, .adj.p=p), 3))",
            "print(sample_stats)" ))
        
        
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
          "iters        <- list(%s)"  %>% sprintf(as.args(iters)),
          "dm_list      <- blply(%s)" %>% sprintf(as.args(dm_args, fun = blply)),
          "sample_stats <- plyr::ldply(dm_list, function (dm) {",
          "  groups <- sample_metadata(biom, %s)[attr(dm, 'Labels')]" %>% sprintf(double_quote(stat.by)),
          "  set.seed(%i)" %>% sprintf(seed),
          sprintf(permutations, fmt = switch(
            EXPR = test,
            adonis2 = "  ptest  <- vegan::adonis2(formula = dm ~ groups, permutations = %i)",
            mrpp    = "  ptest  <- vegan::mrpp(dat = dm, grouping = groups, permutations = %i)" )),
          "  pstats <- summary(vegan::permustats(ptest))",
          "  with(pstats, signif(data.frame(.stat=statistic, .z=z, .p.val=p), 3))",
          "})",
          sprintf("sample_stats[['.adj.p']] <- signif(p.adjust(sample_stats[['.p.val']], %s), 3))", double_quote(p.adj)),
          "print(sample_stats)" ))
        
      })
      
    }
    
    
    
    #________________________________________________________
    # Taxa as groups
    #________________________________________________________
    
    if (!is.null(taxa_stats)) {
      
      attr(sample_coords, 'taxa_stats_cmds') <- local({
        
        
        # Simple case: no iterations over the stats section(s).
        #________________________________________________________
        
        if (length(results) == 1 && length(rank) < 2)
          return (paste(
            sep = "\n",
            "dm         <- %s" %>% fmt_cmd(bdiv_distmat, biom, bdiv, weighted, tree),
            "taxa_mtx   <- %s" %>% fmt_cmd(taxa_matrix, biom, rank, taxa, unc),
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
            sprintf("taxa_stats[['.adj.p']] <- signif(p.adjust(taxa_stats[['.p.val']], %s), 3))", double_quote(p.adj)),
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
          "    taxa_mtx <- %s" %>% fmt_cmd(taxa_matrix, biom_, rank_, taxa, unc),
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
          sprintf("taxa_stats[['.adj.p']] <- signif(p.adjust(taxa_stats[['.p.val']], %s), 3))", double_quote(p.adj)),
          if (p.top < 1)          { sprintf("taxa_stats <- subset(taxa_stats, .adj.p <= %f)", p.top)
          } else if (p.top < Inf) { sprintf("taxa_stats <- subset(taxa_stats, rank(.p.val) <= %i)", p.top)
          } else                  { NULL },
          "print(taxa_stats)" ))
        
      })
    }
    
  }
  
  
  
  set_cache_value(cache_file, sample_coords)
  return (sample_coords)
}




biplot_taxa_stats <- function (dm, taxa_mtx, test, seed, permutations, p.adj, p.top) {
  
  taxa_stats <- plyr::adply(
    .data        = taxa_mtx[attr(dm, 'Labels'),], 
    .margins     = 2L, 
    .id          = ".taxa", 
    .fun         = distmat_stats,
    dm           = dm, 
    test         = test, 
    seed         = seed,
    permutations = permutations,
    p.adj        = p.adj )
  
  if (!is.data.frame(taxa_stats)) return (NULL)
  
  if (p.top >= 1)
    taxa_stats <- taxa_stats[rank(taxa_stats[['.p.val']]) <= p.top,,drop=FALSE]
  
  if (p.top < 1)
    taxa_stats <- taxa_stats[taxa_stats[['.adj.p']] <= p.top,,drop=FALSE]
  
  
  return (taxa_stats)
}




biplot_taxa_coords <- function (sample_coords, taxa_mtx) {
  
  xyz     <- intersect(names(sample_coords), c('.x', '.y', '.z')) %>% setNames(.,.)
  origins <- lapply(xyz, function (i) { round(mean(sample_coords[[i]]), 10) })
  names(origins) %<>% paste0("0")
  taxa_mtx <- taxa_mtx[sample_coords[['.sample']],]
  
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


