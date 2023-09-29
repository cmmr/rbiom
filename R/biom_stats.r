#' Compute p-values for data vs categorical metadata.
#' 
#' @noRd
#' 
#' @param biom   A \code{BIOM} object, as returned from [read_biom()]. 
#'        Alternatively, a data.frame with column names expected by \code{x}, 
#'        \code{y}, and \code{by}. \bold{\emph{Required}}
#'        
#' @param x   A CATEGORICAL metadata column name to use for comparisons.
#'        \bold{\emph{Required}}
#' 
#' @param y   A \emph{character vector} with a value from one of the following 
#'        sets:
#'        \itemize{
#'          \item{\bold{Alpha Diversity Metric - }}{
#'            \code{"OTUs"}, \code{"Shannon"}, \code{"Chao1"}, \code{"Simpson"}, 
#'            and/or \code{"InvSimpson"}.
#'          }
#'          \item{\bold{Beta Diversity Metric - }}{
#'            \code{"Manhattan"}, \code{"Euclidean"}, \code{"Bray-Curtis"}, 
#'            \code{"Jaccard"}, or \code{"UniFrac"}. Use in combination with the
#'            \code{weighted} parameter. Metadata column names can be 
#'            prefixed with \bold{==} or \bold{!=} to limit distance
#'            calculations to \emph{within} or \emph{between}, respectively, 
#'            those categories. See examples below.
#'          }
#'          \item{\bold{Taxonomic Rank - }}{
#'            \code{"Kingdom"}, \code{"Phylum"}, \code{"Class"}, \code{"Order"}, 
#'            \code{"Family"}, \code{"Genus"}, \code{"Species"}, \code{"Strain"}, or 
#'            \code{"OTU"}. Supported ranks will vary by biom. \cr
#'            Run \code{taxa_ranks(biom)} to see the available options.
#'          }
#'        }
#'        \bold{\emph{Required}}
#'        
#' @param by   Additional metadata columns to group by (e.g. facets in a plot).
#'        Default: \code{NULL}
#' 
#' @param pairwise   If \bold{FALSE}, one p-value is generated for each 
#'        \code{by} group, or 1 p-value when \code{by = NULL}. If \bold{TRUE}, 
#'        unique values in \code{x} are compared pairwise against each other. 
#'        Default: \code{FALSE}
#'        
#' @param adj   Adjustment to use for multiple comparisons. See 
#'        \code{stats::p.adjust.methods} for available options. 
#'        Default: \code{"fdr"}
#'        
#' @param weighted   If \bold{TRUE}, run beta diversity metrics in weighted 
#'        mode which take abundances into account. Set to \code{FALSE} to use 
#'        unweighted mode, which only considers presence/absence.
#'        Default: \code{TRUE}
#'        
#' @param digits   Round p-values to this many significant figures. Set to 
#'        \bold{NULL} to disable rounding. Default: \code{3}
#'        
#' @param y.pos   Add an extra column to the results to indicate the upper most 
#'        value of each group. Added for ggpubr::geom_bracket(). Options are 
#'        \bold{NULL} for no y.pos column, \bold{max} for the maximum value, 
#'        \bold{box} for a box plot's the whisker upper bound, or \bold{violin} 
#'        for the highest point of a violin plot. Default: \code{NULL}
#'        
#' @param adj.n   Manually set the number of comparisons to correct for. Only
#'        set this if you know what you are doing! Default: \code{NULL}
#'               
#' @return A data.frame with columns \bold{.p.val} and \bold{.adj.p}, as well 
#'         as columns for tracking \code{x}, \code{y}, and \code{by} values.
#' 
#' @examples
#'     library(rbiom)
#'     
#'     biom <- sample_rarefy(hmp50)
#'     biom <- sample_subset(biom, `Body Site` %in% c('Saliva', 'Stool', 'Buccal mucosa'))
#'       
#'     
#'     # Overall, do the three body sites have different Shannon
#'     #  diversity indices (alpha diversity)? - YES
#'     biom_stats(biom, x = "Body Site", y = "Shannon")
#'     
#'     # Which body sites have the most different Shannon diversity
#'     #  indices? - Buccal mucosa vs. Saliva
#'     biom_stats(biom, x = "Body Site", y = "Shannon", pairwise = TRUE)
#'     
#'     # Do males or females have more alpha diversity variation among 
#'     #  body sites? - Females
#'     biom_stats(biom, x = "Body Site", y = "Shannon", by = "Sex")
#'     
#'     # Which phylum is most different overall between males and females?
#'     #  - Saccharibacteria
#'     biom_stats(biom, x = "Sex", y = "Phylum") %>% head()
#'     
#'     # Which phyla is the most differentially abundant between males and 
#'     #  females on a particular body site? - Tenericutes in Saliva
#'     biom_stats(biom, x = "Sex", y = "Phylum", by = "Body Site") %>% head()
#'     
#'     # Overall, are males and females characterized by their OTU 
#'     #  abundances? - NO
#'     biom_stats(biom, x = "Sex", y = "unifrac")
#'     
#'     # What about on a particular body site? - Yes, Buccal mucosa.
#'     biom_stats(biom, x = "==Sex", y = "unifrac", by = "==Body Site")
#'     
#'     # Recall that for distance metrics such as unifrac, the value of
#'     #  interest (the distance) is calculated from a pair of samples.
#'     #  Sometimes both samples will be from males, sometimes both from
#'     #  females, and sometimes one from each.
#'     # The '==' prefix on '==Sex' and '==Body Site' indicates that only
#'     #  pairs from the same sex and same body site should be considered.
#'     # Using the '!=' prefix does the opposite - considering only sample
#'     #  pairs from different sexes and body sites.
#'     
#'     # Removing the '==' prefix generates the complete set of 
#'     #  comparisons. Note how the test changes from Mann-Whitney to
#'     #  Kruskal-Wallis. This is because instead of a two group comparison
#'     #  ('Male', 'Female') we are now comparing three groups 
#'     #  ('Male', 'Female', 'Female vs Male').
#'     biom_stats(biom, x = "Sex", y = "unifrac", by = "Body Site")
#'
biom_stats <- function (
    biom, x, y, by = NULL, adj = "fdr", 
    pairwise = FALSE, weighted = TRUE, digits = 3, 
    y.pos = NULL, y.pos.facet = "Metric", adj.n = NULL ) {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("biom_stats", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Sanity checks
  #________________________________________________________
  if (length(x) == 0) stop("Please provide a value to 'x'.")
  if (length(y) == 0) stop("Please provide a value to 'y'.")
  
  
  #________________________________________________________
  # Assemble a data frame to run stats on.
  #________________________________________________________
  if (is(biom, 'BIOM')) {
    
    md <- c(x, by)
    md <- md[sub("^[\\!=]=", "", md) %in% colnames(sample_metadata(biom))]
    df <- biom_distill(biom, metric = y, weighted = weighted, long = TRUE, md = md)
    y  <- ".value"
    
    mode <- attr(df, 'mode', exact = TRUE)
    if (mode == 'adiv') {
      by <- c(".adiv", by)
      
    } else if (mode == 'taxa') {
      by <- c(names(df)[[2]], by)
      
    } else if (mode == 'bdiv') {
      by <- sub("^[\\!=]=", "", by)
      x  <- sub("^[\\!=]=", "", x)
    }
    
  } else if (is(biom, 'data.frame')) {
    
    #________________________________________________________
    # The user supplied the data frame to run stats on
    #________________________________________________________
    
    df <- biom
    
    missing <- setdiff(c(x, y, by), names(df))
    if (length(missing) > 0)
      stop(
        "Expected column not present in data frame: ", 
        paste(collapse = ", ", missing) )
    
  } else {
    
    #________________________________________________________
    # The user gave us an unexpected data type to run stats on
    #________________________________________________________
    
    stop(
      "'biom' should be BIOM or data.frame, not ", 
      paste(collapse=", ", class(biom)) )
  }
  
  
  #________________________________________________________
  # Create a `.x` column if there's multiple x values
  #________________________________________________________
  if (length(x) > 1) {
    df[['.x']] <- apply(df[,x,drop=FALSE], 1L, paste, collapse=".")
    x <- ".x"
  }
  
  
  #________________________________________________________
  # No data to run stats on.
  #________________________________________________________
  if (nrow(df) == 0)
    return (NULL)
  
  
  #________________________________________________________
  # Run non-parametric statistics
  #________________________________________________________
  plyby   <- ply_cols(rev(by))
  results <- plyr::ddply(df, plyby, function (z) {
    
    if (isTRUE(pairwise)) {
      
      #________________________________________________________
      # wilcox-style stats
      #________________________________________________________
      groups <- if (is.factor(z[[x]])) levels(z[[x]]) else sort(unique(z[[x]]))
      
      if (length(groups) < 2) {
        
        result <- data.frame(
          check.names = FALSE, 
          n      = nrow(z),
          Test   = NA,
          Group1 = NA,
          Group2 = NA,
          .p.val = 1 )
        
      } else {
        
        pairs <- combn(groups, 2)
        
        result <- data.frame(
          check.names = FALSE,
          n      = apply(pairs, 2L, function (pair) { sum(z[[x]] %in% pair) }),
          Test   = "Mann-Whitney",
          Group1 = pairs[1,],
          Group2 = pairs[2,],
          .p.val = apply(pairs, 2L, function (pair) {
            
            tryCatch({
              wt <- wilcox.test(
                x = z[z[[x]] == pair[[1]], y],
                y = z[z[[x]] == pair[[2]], y] ) %>%
                suppressWarnings()
              wt[['p.value']]
            }, error = function (e) 1)
          }) )
        
      }
      
      
    } else {
      #________________________________________________________
      # kruskal-style stats, possibly interspersed with wilcox
      #________________________________________________________
      
      ngroups <- length(unique(z[[x]]))
      result  <- data.frame(n = nrow(z))
      
      if (ngroups < 2) {
        result[['Test']]   <- NA
        result[['.p.val']] <- 1
        
      } else if (ngroups == 2) {
        result[['Test']]   <- "Mann-Whitney"
        result[['.p.val']] <- tryCatch({
          groups <- as.character(unique(z[[x]]))
          wt <- wilcox.test(
            x = z[z[[x]] == groups[[1]], y],
            y = z[z[[x]] == groups[[2]], y] ) %>%
            suppressWarnings()
          wt[['p.value']]
        }, error = function (e) 1)
        
      } else {
        result[['Test']]   <- "Kruskal-Wallis"
        result[['.p.val']] <- tryCatch({
          kruskal.test(z[[y]], z[[x]])[['p.value']]
        }, error = function (e) 1)
        
      }
    }
    
    
    #________________________________________________________
    # Vertical start position for p-value plot annotations
    #________________________________________________________
    
    if (!is_null(y.pos)) {
      
      if ("vline" %in% names(attributes(biom))) {
        vline <- attr(biom, 'vline', exact = TRUE)
        vline <- setNames(vline[['.ymax']], vline[['.group']])
        z[['.vline']] <- vline[z[['.group']]]
        remove("vline")
      }
      
      result[['y.pos']] <- local({
        
        max(plyr::daply(z, ".group", function (group_df) {
          
          group_df <- group_df[is.finite(group_df[[y]]),,drop=FALSE]
          vals     <- group_df[[y]]
          
          if (length(vals) == 0) return (0)
          
          res <- 0
          for (i in y.pos)
            res <- max(na.rm = TRUE, res, local({
              if (i == "max")           { max(vals)
              } else if (i == "mean")   { mean(vals)
              } else if (i == "box")    { boxplot.stats(vals)$stats[5]
              } else if (i == "violin") { max(density(vals)[['x']])
              } else if (i == "vline")  { max(group_df[['.vline']]) } }))
          
          # Avoid horizontal gridlines.
          res <- max(labeling::extended(0, res, 5, only.loose = TRUE))
          
          return (res)
        }))
      })
    }
    
    return (result)
    
  })
  
  # Remove '.id' column that ddply adds when by=NULL
  if (length(by) == 0)
    results <- results[,-1,drop=F]
  
  # Don't count multiple metrics as multiple comparisons
  if (any(c(".adiv", ".bdiv") %in% by)) {
    metric  <- intersect(c(".adiv", ".bdiv"), by)
    results <- plyr::ddply(results, metric, function (z) {
      z[['.adj.p']] <- p.adjust(
        p      = z[['.p.val']], 
        method = adj   %||% "fdr", 
        n      = adj.n %||% nrow(z) )
      return(z)
    })
  } else {
    results[['.adj.p']] <- p.adjust(
      p      = results[['.p.val']], 
      method = adj   %||% "fdr", 
      n      = adj.n %||% nrow(results) )
  }
  
  # Order by p-value, most to least significant (unless plotting)
  if (is_null(y.pos))
    results <- results[order(results[['.p.val']]),,drop=F]
  
  
  # Round the p.val and adj.p numbers
  if (!is_null(digits)) {
    results[['.p.val']] <- signif(results[['.p.val']], digits = digits)
    results[['.adj.p']] <- signif(results[['.adj.p']], digits = digits)
  }
  
  # Make y.pos the last column
  if (!is_null(y.pos) && "y.pos" %in% names(results)) {
    
    y.pos.facet <- intersect(y.pos.facet, by)
    
    if (is_null(y.pos.facet)) {
      
      results[['y.pos']] <- signif(max(results[['y.pos']]), digits = 2)
      
    } else {
      plyby   <- ply_cols(rev(y.pos.facet))
      results <- plyr::ddply(results, plyby, function (z) {
        z[['y.pos']] <- signif(max(z[['y.pos']]), digits = 2)
        return (z)
      })
    }
    
    results <- results[,order(names(results) == "y.pos"),drop=F]
  }
  
  rownames(results) <- NULL
  
  
  set_cache_value(cache_file, result)
  return (results)
}
