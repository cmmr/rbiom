#' Reduce a distance matrix to two or three dimensions
#' 
#' @name bdiv_ord_table
#' 
#' @param biom   A \code{BIOM} object, as returned from \link{read_biom}.
#'        Alternatively, a distance matrix, such as from \link{bdiv_distmat}.
#'        
#' @param dist   The distance algorithm to use. Options are:
#'        \code{"Bray-Curtis"}, \code{"Manhattan"}, \code{"Euclidean"}, 
#'        \code{"Jaccard"}, and \code{"UniFrac"}. A phylogenetic tree must be
#'        present in \code{biom} or explicitly provided via \code{tree=} to use
#'        the UniFrac methods. Ignored if \code{biom} is a distance matrix.
#'        Multiple values allowed. Default: \code{"Bray-Curtis"}.
#'                
#' @param ord    Method for reducing the dimensionality of a distance matrix.
#'        Multiple values allowed. Options are:
#'        \describe{
#'            \item{\code{"UMAP"} (default)}{
#'              Uniform manifold approximation and projection,
#'              via \code{\link[uwot]{umap}}. }
#'            \item{\code{"PCoA"}}{
#'              Principal coordinate analysis,
#'              via \code{\link[ape]{pcoa}}. }
#'            \item{\code{"NMDS"}}{
#'              Nonmetric multidimensional scaling,
#'              via \code{\link[vegan]{metaMDS}}. }
#'            \item{\code{"tSNE"}}{
#'              t-distributed stochastic neighbor embedding,
#'              via \code{\link[tsne]{tsne}}. }
#'           }
#'        
#' @param weighted   Take relative abundances into account. When 
#'        \code{weighted=FALSE}, only presence/absence is considered.
#'         Ignored if \code{biom} is a distance matrix. Multiple values 
#'         allowed. (Default: \code{TRUE})
#'     
#' @param tree   A \code{phylo} object representing the phylogenetic
#'        relationships of the taxa in \code{biom}. Will be taken from the tree
#'        embedded in the \code{biom} object if not explicitly specified. Only
#'        required for computing UniFrac distance matrices.  Ignored if 
#'        \code{biom} is a distance matrix.
#'     
#' @param md   Include metadata in the output data frame? Ignored if 
#'        \code{biom} is a distance matrix. Options are: 
#'        \describe{
#'          \item{\code{NULL}}{ Don't include metadata. (Default) }
#'          \item{\code{TRUE}}{ Include all metadata. }
#'          \item{\emph{character vector}}{ Include only the specified metadata columns. }
#'        }
#' @param k   Number of dimensions to return. (Default: \code{2})
#'        
#' @param split.by  Name(s) of metadata columns that the data should be split
#'        by prior to calculating distance matrices, ordinations, and statistics.
#'        
#' @param stat.by,seed,perms   Passthrough parameters to \code{bdiv_distmat()} 
#'        for computing adonis statistics. (Default: \code{stat.by=NULL, seed=0, 
#'        perms=1000})
#'        
#' @param rank   [Biplot] What rank of taxa to display (e.g. "Phylum"), or 
#'        \code{NULL} for no biplot. Run \code{taxa_ranks()} to see all options 
#'        for a given BIOM object. (Default: \code{NULL})
#'        
#' @param taxa   [Biplot] Which taxa to display. An integer value will return 
#'        the top n most abundant taxa. A value 0 <= n < 1 will return any taxa 
#'        with that mean abundance or greater (e.g. 0.1). A character vector of
#'        taxon names will show only those taxa. (Default: \code{5})
#'
#' @param p.top   [Biplot] Only returns taxa with the most significant 
#'        differences in abundance. If \code{p.top} is >= 1, then the 
#'        \code{p.top} most significant taxa are displayed. If \code{p.top} 
#'        is less than one, all taxa with an adjusted p-value <= \code{p.top} 
#'        are displayed. Recommended to be used in combination with the 
#'        \code{taxa} parameter to set a lower bound on the mean abundance of 
#'        considered taxa. (Default: \code{Inf})
#'
#' @param p.adj   [Biplot] Method to use for multiple comparisons adjustment of 
#'        p-values. Run \code{p.adjust.methods} for a list of available options.
#'        (Default: \code{fdr})
#'        
#' @param ...  Additional arguments to pass on to \code{\link[ape]{pcoa}}, 
#'             \code{\link[vegan]{metaMDS}}, or \code{\link[tsne]{tsne}}.
#'             
#' @return A data.frame with columns .axis.1, .axis.2, ..., axis.k as well as 
#'         columns given by \code{md}, \code{split.by}, and \code{stat.by}. 
#'         Sample IDs are in column \code{'.sample'} and the metric used in 
#'         column \code{'.metric'}. For pcoa ordinations, 
#'         \code{attr(, 'eig')} will contain the eigenvalues useful for 
#'         construction of ".. % variation explained" labels. 
#'         \code{attr(, 'stats_raw')} and \code{attr(,'stats_tbl')} contain the 
#'         raw and tabular adonis statistics when \code{stat.by} is set. 
#'         Distance matrices are in \code{attr(, 'dm')}. When \code{rank} is
#'         non-NULL, or \code{taxa} is a list of taxon names, a data.frame of
#'         biplot coordinates will be given by \code{attr(,'biplot')}.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     ord <- bdiv_ord_table(hmp50, "bray-curtis", "pcoa")
#'     head(ord)
#'     

bdiv_ord_table <- function (
    biom, dist="Bray-Curtis", ord="UMAP", weighted=TRUE, tree=NULL, 
    md=NULL, k=2, split.by=NULL, stat.by=NULL, seed=0, perms=1000, 
    rank=NULL, taxa=5, p.adj='fdr', p.top=5, ...) {
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(c(as.list(environment()), list(...)), eval)
  cache_file <- get_cache_file("bdiv_ord_table", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  ord_args <- list(...)
  
  ord  %<>% validate_metrics(NULL, ., mode = "ord",  multi=TRUE)
  dist %<>% validate_metrics(biom, ., mode = "bdiv", multi=TRUE)
  
  
  #________________________________________________________
  # Can't do as much when biom is a dist object
  #________________________________________________________
  if (!is(biom, 'BIOM'))
    for (i in c('md', 'split.by', 'stat.by', 'rank'))
      if (!is_null(get(i)))
        stop("BIOM object needed for running bdiv_ord_table() with `", i, "` argument.")
  
  
  #________________________________________________________
  # Verify metadata columns exist
  #________________________________________________________
  if (isTRUE(md))
    md <- colnames(metadata(biom))
  
  if (is.character(md)) {
    missing <- setdiff(c(md, split.by, stat.by), colnames(metadata(biom)))
    if (length(missing) > 0)
      stop("Metadata column(s) ", paste(collapse = ", ", missing), " not found.")
  }
  
  
  #________________________________________________________
  # Subdivide BIOM object according to metadata
  #________________________________________________________
  if (is_null(split.by)) {
    biom_list <- list(biom)
  } else {
    biom_list <- blply(biom, split.by, function (b) b)
  }
  
  
  #________________________________________________________
  # Biplot - determine taxonomic rank
  #________________________________________________________
  if (is_null(rank) && is.character(taxa)) {
    rank <- names(which.max(apply(taxonomy(biom), 2L, function (x) sum(x %in% taxa))))
  } else if (!is_null(rank)) {
    rank <- validate_metrics(biom, rank, mode="rank", multi=TRUE)
  }
  
  
  #________________________________________________________
  # Aggregate iterative output together into final tables
  #________________________________________________________
  add_loop_columns <- function (x, i, w, d, o = NULL) {
    if (is_null(x)) return (NULL)
    x[['.weight']] <- ifelse(w, "Weighted", "Unweighted")
    x[['.dist']]   <- d
    if (!is_null(o)) x[['.ord']]    <- o
    for (j in colnames(attr(biom_list, 'split_labels')))
      x[[j]] <- attr(biom_list, 'split_labels')[i,j]
    return (x)
  }
  
  
  #________________________________________________________
  # Loop over all combinations of variables
  #________________________________________________________
  ord_results       <- NULL
  stats_tbl_results <- NULL
  biplot_results    <- NULL
  
  for (i in seq_along(biom_list))
    for (w in as.logical(weighted))
      for (d in dist) {
        
        #________________________________________________________
        # Generate distance matrix and/or Adonis statistics
        #________________________________________________________
        b <- biom_list[[i]]
        
        if (is(b, 'dist')) {
          dm <- b
          if (!is_null(stat.by)) {
            set.seed(seed)
            stats_raw <- try(vegan::adonis2(dm ~ stat.by, permutations=perms-1), silent = TRUE)
            stats_tbl <- adonis_table(stats_raw)
          } else {
            stats_tbl <- NULL
          }
          
        } else {
          dm <- bdiv_distmat(biom = b, method = d, weighted = w, tree = tree, stat.by=stat.by, seed=seed, perms=perms-1)
          stats_tbl <- attr(dm, 'stats_tbl', exact = TRUE)
          attr(dm, 'stats_raw') <- NULL
          attr(dm, 'stats_tbl') <- NULL
        }
        
        
        for (o in ord) {
          
          
          #________________________________________________________
          # Ordinate samples => ('.axis.1', '.axis.2')
          #________________________________________________________
          if (o == "PCoA") {
            res <- do.call(ape::pcoa, c(list(D = dm), ord_args))
            df  <- res$vectors[,1:k] %>% as.data.frame()
            eig <- res$values$Relative_eig[1:k] %>% as.data.frame()
            
          } else if (o == "tSNE") {
            ord_args[['X']] <-  dm
            ord_args[['k']] <-  k
            res <- suppressMessages(do.call(tsne::tsne, ord_args))
            df  <- res %>% as.data.frame()
            rownames(df) <- attr(dm, "Labels", exact = TRUE)
            eig <- NULL
            
          } else if (o == "NMDS") {
            ord_args[['comm']]   <-  dm
            ord_args[['k']]      <-  k
            ord_args[['trace']] %<>% if.null(0)
            res <- do.call(vegan::metaMDS, ord_args)
            df  <- res$points %>% as.data.frame()
            eig <- NULL
            
          } else if (o == "UMAP") {
            ord_args[['X']]             <-  dm
            ord_args[['n_components']] %<>% if.null(k)
            ord_args[['n_neighbors']]  %<>% if.null(max(2, min(100, as.integer(attr(dm, 'Size') / 3))))
            res <- do.call(uwot::umap, ord_args)
            df  <- res %>% as.data.frame()
            eig <- NULL
            
          } else {
            stop("'", o, "' is not a valid argument for 'ord'.")
          }
          colnames(df)    <- paste0(".axis.", 1:k)
          df[['.sample']] <- rownames(df)
          
          
          #________________________________________________________
          # Calculate 0, 1, or more biplot(s)
          #________________________________________________________
          biplot <- ordinate_biplot(b, df, rank, taxa, p.adj, p.top, perms)
          
          rownames(df) <- NULL
          ord_results    %<>% rbind(add_loop_columns(df,     i, w, d, o))
          biplot_results %<>% rbind(add_loop_columns(biplot, i, w, d, o))
          
        } # end o / ord
        
        stats_tbl_results %<>% rbind(add_loop_columns(stats_tbl, i, w, d))
        
      } # end i,w,d / biom_list, weighted, dist
  
  
  #________________________________________________________
  # Make properly ordered factors
  #________________________________________________________
  strings_to_factors <- function (x) {
    if (is_null(x)) return (NULL)
    for (i in intersect(colnames(x), c('.weight', '.dist', '.ord')))
      x[[i]] %<>% factor(levels = unique(x[[i]]))
    return (x)
  }
  ord_results       %<>% strings_to_factors()
  stats_tbl_results %<>% strings_to_factors()
  biplot_results    %<>% strings_to_factors()
  
  
  #________________________________________________________
  # Add metadata to the ord_results table.
  #________________________________________________________
  sampleIDs <- ord_results[['.sample']]
  for (i in rev(unique(c(md, split.by, stat.by))))
    ord_results[[i]] <- metadata(biom, i)[sampleIDs]
  remove("sampleIDs")
  
  
  attr(ord_results, 'response') <- ".diversity"
  attr(ord_results, 'stats_tbl') <- stats_tbl_results
  attr(ord_results, 'biplot')    <- biplot_results
  
  
  set_cache_value(cache_file, ord_results)
  return (ord_results)
}



ordinate_biplot <- function (biom, coords, ranks, taxa, p.adj, p.top, perms) {
  
  if (is_null(ranks)) return (NULL)
  
  #________________________________________________________
  # Handle multiple ranks.
  #________________________________________________________
  biplots <- plyr::ldply(ranks, function (rank) {
    
    #________________________________________________________
    # Calculate Weighted average for each taxon
    #________________________________________________________
    biplot                 <- distill(biom, rank, md = FALSE)
    biplot                 <- rename_response(biplot, ".abundance")
    biplot[['.rank']]      <- rank
    biplot[['.abundance']] <- biplot[['.abundance']] / sample_sums(biom)[biplot[['.sample']]]
    biplot[['.x']]         <- coords[biplot[['.sample']], '.axis.1']
    biplot[['.y']]         <- coords[biplot[['.sample']], '.axis.2']
    center.x               <- round(mean(biplot[['.x']]), 10)
    center.y               <- round(mean(biplot[['.y']]), 10)
    
    
    #________________________________________________________
    # Subset the taxa. By name, top n, or min abundance.
    #________________________________________________________
    if (!is_null(taxa)) {
      
      if (is.numeric(taxa) && length(taxa) == 1) {
        if (taxa >= 1) {
          taxa <- top_taxa(biom, rank, taxa)
        } else {
          ts   <- as_percent(biom) %>% taxa_means(rank)
          taxa <- which(ts >= taxa) %>% names()
        }
      }
      
      if (is.character(taxa)) {
        biplot <- biplot[biplot[['.taxa']] %in% taxa,,drop=F]
      }
    }
    
    
    #________________________________________________________
    # Get x,y,size and p-value for each taxon
    #________________________________________________________
    biplot <- plyr::ddply(biplot, c('.rank', '.taxa'), function(df) {
      
      # Calculate Weighted average for each taxon
      Abundance <- sum(df[['.abundance']])
      if (Abundance == 0) return (NULL)
      Axis.1 <- sum(df[['.x']] * df[['.abundance']]) / Abundance
      Axis.2 <- sum(df[['.y']] * df[['.abundance']]) / Abundance
      
      data.frame(
        '.axis.1'    = Axis.1, 
        '.axis.2'    = Axis.2,
        '.ori.1'     = center.x,
        '.ori.2'     = center.y,
        '.abundance' = Abundance/nrow(df), 
        '.p.val'     = local({
          
          # Use random reassignment of abundances to calculate 
          # p-value. Known as approximate/random/Monte Carlo
          # permutation test.
          ptest <- boot::boot(
            data      = df, 
            R         = perms, 
            statistic = function(data, ind) {
              r.x <- sum(data[['.x']] * data[ind,'.abundance']) / Abundance
              r.y <- sum(data[['.y']] * data[ind,'.abundance']) / Abundance
              d   <- sqrt((r.x - center.x)^2 + (r.y - center.y)^2)
              return (d)
            })
          d <- sqrt((Axis.1 - center.x)^2 + (Axis.2 - center.y)^2)
          p <- length(which(ptest[['t']] >= d)) / ptest[['R']]
          
          return (p)
          
        }) )
    })
    
    
    #________________________________________________________
    # Limit to taxa with adj. p-value <= p.top or top n signif taxa
    #________________________________________________________
    biplot[['.adj.p']] <- p.adjust(biplot[['.p.val']], method=p.adj)
    if (p.top >= 1) {
      biplot <- biplot[rank(-biplot[['.adj.p']]) <= p.top,,drop=F]
    } else {
      biplot <- biplot[biplot[['.adj.p']] <= p.top,,drop=F]
    }
    
    
    #________________________________________________________
    # Final biplot data frame
    #________________________________________________________
    if (nrow(biplot) == 0)
      return (NULL)
    
    return (biplot)
  })
  
  return (biplots)
}



