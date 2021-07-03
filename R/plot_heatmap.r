
# colors = Ramped, then set as 'color'. Default is blue -> white -> red.
# color  = Fed directly into pheatmap's 'color' for cell coloring. Overrides 'colors'.
# color.by = Which metadata columns to add as tracks above the heatmap.
# annotation_colors = How to color the tracks above the heatmap.

# "rank ~ Samples" ONLY:
#  - normalize.rows = accentuate where each taxa is most abundant
#  - taxa  = number of top taxa (e.g. 9); minimum mean abundance (e.g. 0.02); taxa list
#  - other = rollup remaining taxa into an "Other" taxa and display that too
#  - abbr  = allow abbreviated 'taxa' values, e.g. c('staph', 'lact')


plot_heatmap <- function (biom, x, y, color.by = NULL, colors = c("#2E5473", "#FFFFFF", "#BD3230"), normalize.rows = TRUE, taxa = 10, other = FALSE, abbr = TRUE, ...) {
  
  params <- c(as.list(environment()), list(...))
  mode   <- paste(attr(y, 'mode', exact = TRUE), "~", attr(x, 'mode', exact = TRUE))
  
  
  if (mode == "bdiv ~ Samples") {
    dm   <- beta.div(biom, y)
    mat  <- as.matrix(dm)
    args <- list(
      mat                      = mat,
      clustering_distance_rows = dm, 
      clustering_distance_cols = dm,
      show_rownames            = FALSE,
      show_colnames            = FALSE
    )
    
  } else if (mode == "rank ~ Samples") {
    mat <- as.matrix(taxa.rollup(biom, rank = y, sparse = T))
    
    # Limit to top n taxa, n relative abundance, or specific list of taxa
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
      
      # Sum the remaining taxa into an "Other" category
      if (length(taxa) == nrow(mat) || isFALSE(other)) {
        mat <- mat[taxa,,drop=F]
        
      } else {
        other <- matrix(
          data = colSums(mat[!rownames(mat) %in% taxa,,drop=F]), 
          nrow = 1, dimnames = list("Other", colnames(mat)))
        mat <- rbind(mat[taxa,,drop=F], other)
      }
    }
    
    # Accentuate in which samples a taxa is most abundant.
    if (normalize.rows) {
      mat <- t(apply(mat, 1, function(x)(x-min(x))/(max(x)-min(x))))
      mat[is.nan(mat)] = 0
    }
    
    args <- list(mat = mat)
  }
  
  
  #--------------------------------------------------------------
  # Default arguments for pheatmap common to all modes.
  #--------------------------------------------------------------
  args[['color']] <- ifelse(params[['xlab.angle']] == "auto", 315, params[['xlab.angle']])
  if (!is.null(colors))   args[['color']]          <- colorRampPalette(colors)(50)
  if (!is.null(color.by)) args[['annotation_col']] <- metadata(biom)[,color.by,drop=F]
  
  
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
