
# gradient = Ramped, then set as 'color'.
# color    = Fed directly into pheatmap's 'color' for cell coloring. Overrides 'gradient'.
# color.by = Which metadata columns to add as tracks above the heatmap.
# colors   = How to color the tracks above the heatmap (annotation_colors).

# clustering_method overrides x
# clust = alias for clustering_method
# dist  = alias for clustering_distance_rows and clustering_distance_cols
#   - accepts "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"


# "rank ~ heatmap" ONLY:
#  - normalize.rows = accentuate where each taxa is most abundant
#  - taxa  = number of top taxa (e.g. 9); minimum mean abundance (e.g. 0.02); taxa list
#  - other = rollup remaining taxa into an "Other" taxa and display that too
#  - abbr  = allow abbreviated 'taxa' values, e.g. c('staph', 'lact')


plot_heatmap <- function (
  biom, x, y, 
  color.by = NULL, colors = NULL, gradient = heat.colors(20), 
  normalize.rows = TRUE, taxa = 10, other = FALSE, abbr = TRUE, weighted = TRUE, 
  clust = as.vector(x), clustering_method = clust, 
  dist = "euclidean", clustering_distance_rows = dist, clustering_distance_cols = dist, ...) {
  
  params <- c(as.list(environment()), list(...))
  mode   <- attr(y, 'mode', exact = TRUE)
  
  
  
  #--------------------------------------------------------------
  # Plot's title will include a clean name for clustering method
  #--------------------------------------------------------------
  clustering_method <- validate_metrics(NULL, clustering_method, 'clust')
  
  common2hclust <- c(
    "heatmap" = "complete", "ward"  = "ward.d2",
    "UPGMA"   = "average",  "WPGMA" = "mcquitty",
    "WPGMC"   = "median",   "UPGMC" = "centroid" )
  if (clustering_method %in% names(common2hclust))
    clustering_method <- unname(common2hclust[clustering_method])
  params[['clustering_method']] <- clustering_method
  
  hclust2title <- c(
    "average"  = "average (UPGMA)",  "ward.d2"  = "Ward's", 
    "mcquitty" = "mcquitty (WPGMA)", "single"   = "single linkage", 
    "median"   = "median (WPGMC)",   "complete" = "complete linkage", 
    "centroid" = "centroid (UPGMC)" )
  cl_title <- unname(hclust2title[clustering_method])
  remove("common2hclust", "hclust2title")
  
  
  #--------------------------------------------------------------
  # Set-up specific to 'bdiv' or 'rank' modes
  #--------------------------------------------------------------
  if (mode == "bdiv") {
    dm   <- beta.div(biom = biom, method = y, weighted = weighted)
    mat  <- as.matrix(dm)
    args <- list(
      mat                      = mat,
      clustering_distance_rows = dm, 
      clustering_distance_cols = dm,
      show_rownames            = FALSE,
      show_colnames            = FALSE,
      main                     = paste0(
        ifelse(weighted, "W", "Unw"), "eighted ", y, ", ", 
        cl_title, " clustering" )
    )
    
  } else if (mode == "rank") {
    mat      <- as.matrix(taxa.rollup(biom, rank = y, sparse = T))
    gradient <- rev(gradient)
    
    params[['clustering_distance_rows']] <- validate_metrics(NULL, clustering_distance_rows, 'dist')
    params[['clustering_distance_cols']] <- validate_metrics(NULL, clustering_distance_cols, 'dist')
    
    #--------------------------------------------------------------
    # Limit to top n taxa, n relative abundance, or custom list
    #--------------------------------------------------------------
    if (!is.null(taxa)) {
      
      if (is.numeric(taxa)) {
        rel <- sort(rowMeans(t(t(mat) / colSums(mat))), decreasing = TRUE)
        if (taxa >= 1) { taxa <- head(names(rel), taxa) 
        } else         { taxa <- names(rel)[rel >= taxa] }
        
      } else if (isTRUE(abbr)) {
        
        opts <- tolower(rownames(mat))
        taxa <- tolower(as.character(taxa))
        taxa <- sapply(taxa, startsWith, x = opts)
        taxa <- rownames(mat)[apply(taxa, 1L, any)]
        
      } else {
        taxa <- intersect(as.character(taxa), rownames(mat))
      }
      
      if (length(taxa) == 0)
        stop("No taxa match the criteria: ", capture.output(str(params[['taxa']])))
      
      #--------------------------------------------------------------
      # Sum the remaining taxa into an "Other" category
      #--------------------------------------------------------------
      if (length(taxa) == nrow(mat) || isFALSE(other)) {
        mat <- mat[taxa,,drop=F]
        
      } else {
        other <- matrix(
          data = colSums(mat[!rownames(mat) %in% taxa,,drop=F]), 
          nrow = 1, dimnames = list("Other", colnames(mat)))
        mat <- rbind(mat[taxa,,drop=F], other)
      }
    }
    
    #--------------------------------------------------------------
    # Accentuate in which samples a taxa is most abundant.
    #--------------------------------------------------------------
    if (normalize.rows) {
      mat <- t(apply(mat, 1, function(x)(x-min(x))/(max(x)-min(x))))
      mat[is.nan(mat)] = 0
    }
    
    args <- list(
      mat  = mat,
      main = paste0(
        ifelse(normalize.rows, "Normalized ", "Absolute "),
        y, " Abundances per Sample\n",
        paste(collapse=" x ", unique(params[paste0('clustering_distance_', c('rows', 'cols'))])),
        " distance, ", cl_title, " clustering" ))
    
  }
  
  
  #--------------------------------------------------------------
  # Default arguments for pheatmap common to all modes.
  #--------------------------------------------------------------
  if (!is.null(angle <- params[['xlab.angle']])) {
    opts <- c(0, 45, 90, 270, 315)
    if (!is.numeric(angle)) angle <- 315
    args[['angle_col']] <- opts[which.min(abs(opts - angle))]
  }
  if (!is.null(gradient)) args[['color']]             <- colorRampPalette(gradient)(50)
  if (!is.null(color.by)) args[['annotation_col']]    <- metadata(biom)[,color.by,drop=F]
  if (!is.null(colors))   args[['annotation_colors']] <- colors
  
  
  #--------------------------------------------------------------
  # Pass along pheatmap arguments, e.g. 'show_rownames'=TRUE
  #--------------------------------------------------------------
  for (i in intersect(names(params), formalArgs(pheatmap::pheatmap)))
    args[[i]] <- params[[i]]
  
  
  #--------------------------------------------------------------
  # Prevent "must have n >= 2 objects to cluster" error messages
  #--------------------------------------------------------------
  if (nrow(args[['mat']]) < 3) args[['cluster_rows']] <- FALSE
  if (ncol(args[['mat']]) < 3) args[['cluster_cols']] <- FALSE
  
  
  do.call(pheatmap::pheatmap, args)
}
