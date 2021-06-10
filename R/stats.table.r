#' Compute p-values for data vs metadata.
#' 
#' @name stats.table
#' @param biom  A \code{BIOM} object, as returned from \link{read.biom}. Alternatively,
#'                a data.frame with column names expected by \code{x}, \code{y}, and \code{by}.
#' @param x     The metadata column name of the category to use for comparisons.
#' @param y     A \emph{character vector} with values from one of the following sets:
#'     \describe{
#'         \item{Alpha Diversity Metrics (one or more)}{
#'           \bold{OTUs}, \bold{Shannon}, \bold{Chao1}, \bold{Simpson}, and/or \bold{InvSimpson}.
#'           You may also set \code{y = "Diversity"} to get all five metrics.
#'         }
#'         \item{Beta Diversity Metrics (one only)}{
#'           \bold{manhattan}, \bold{euclidean}, \bold{bray-curtis}, \bold{jaccard}, 
#'           or \bold{unifrac}. Use in combination with the \code{weighted} parameter.
#'           Metadata column names can be prefixed with \bold{==} or \bold{!=} to limit distance
#'           calculations to \emph{within} or \emph{between}, respectively, those categories. See
#'           examples below. Setting \code{y = "Distance"} will use \bold{unifrac} if a 
#'           phylogenetic tree is present, or \bold{bray-curtis} otherwise.
#'         }
#'         \item{Taxa Abundances (one only)}{
#'           \bold{Kingdom}, \bold{Phylum}, \bold{Class}, \bold{Order}, \bold{Family}, \bold{Genus}, 
#'           \bold{Species}, \bold{Strain}, or \bold{OTU}. Supported ranks will vary by biom. Run
#'           \code{taxa.ranks(biom)} to see the available options. Specifying \code{y = "Taxa"}
#'           will default to the most precise rank possible.
#'         }
#'        }
#' @param by       Additional metadata columns to group by (e.g. facets in a plot).
#' @param pairwise  If \bold{FALSE} (the default), one p-value is generated for each \code{by} 
#'                  group, or 1 p-value when \code{by = NULL}. If \bold{TRUE}, unique values in
#'                  \code{x} are compared pairwise against each other.
#' @param adj      Adjustment to use for multiple comparisons. See \code{stats::p.adjust.methods} 
#'                 for available options. Default: \bold{fdr}.
#' @param weighted   If \bold{TRUE}(the default), run beta diversity metrics in weighted mode. 
#'                   Set to \bold{FALSE} to use unweighted mode.
#' @param digits   Round p-values to this many significant figures (default 3). 
#'                 Set to \bold{NULL} to disable rounding.
#' @param y.pos  Add an extra column to the results to indicate the upper most value of each group.
#'               Added for ggpubr::geom_bracket(). Options are \bold{NULL} for no y.pos column 
#'               (the default), \bold{max} for the maximum value, \bold{box} for a box plot's the 
#'               whisker upper bound, or \bold{violin} for the highest point of a violin plot.
#' @return A data.frame with columns \bold{p.val}, \bold{adj.p}, \bold{adj.sig}, as well as
#'                 columns for tracking for \code{x}, \code{y}, and \code{by} categories.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom   <- rarefy(read.biom(infile))
#'     biom   <- subset(biom, `Body Site` %in% c('Saliva', 'Stool', 'Buccal mucosa'))
#'       
#'     
#'     # Overall, do the three body sites different have different Shannon
#'     #  diversity indices (alpha diversity)? - YES
#'     stats.table(biom, x = "Body Site", y = "Shannon")
#'     
#'     # Which body sites have the most different Shannon diversity
#'     #  indices? - Buccal mucosa vs. Saliva
#'     stats.table(biom, x = "Body Site", y = "Shannon", pairwise = TRUE)
#'     
#'     # Do males or females have more alpha diversity variation among 
#'     #  body sites? - Females
#'     stats.table(biom, x = "Body Site", y = "Shannon", by = "Sex")
#'     
#'     # Which phylum is most different overall between males and females?
#'     #  - Saccharibacteria
#'     stats.table(biom, x = "Sex", y = "Phylum") %>% head()
#'     
#'     # Which phyla is the most differentially abundant between males and 
#'     #  females on a particular body site? - Tenericutes in Saliva
#'     stats.table(biom, x = "Sex", y = "Phylum", by = "Body Site") %>% head()
#'     
#'     # Overall, are males and females characterized by their OTU 
#'     #  abundances? - NO
#'     stats.table(biom, x = "Sex", y = "unifrac")
#'     
#'     # What about on a particular body site? - Yes, Buccal mucosa.
#'     stats.table(biom, x = "==Sex", y = "unifrac", by = "==Body Site")
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
#'     stats.table(biom, x = "Sex", y = "unifrac", by = "Body Site")
#'
stats.table <- function (biom, x, y, by = NULL, adj = "fdr", pairwise = FALSE, weighted = TRUE, digits = 3, y.pos = NULL, y.pos.facet="Metric") {
  
  #--------------------------------------------------------------
  # Sanity checks
  #--------------------------------------------------------------
  if (length(x) == 0) stop("Please provide a value to 'x'.")
  if (length(y) == 0) stop("Please provide a value to 'y'.")
  
  
  #--------------------------------------------------------------
  # Assemble a data frame to run stats on.
  #--------------------------------------------------------------
  if (is(biom, 'BIOM')) {
    
    df <- distill(biom, metric = y, weighted = weighted, long = TRUE, md = c(x, by))
    y  <- attr(df, 'response')
    
    if (y == 'Diversity') {
      by <- c("Metric", by)
      
    } else if (y == 'Abundance') {
      by <- c(names(df)[[2]], by)
      
    } else if (y == 'Distance') {
      by <- sub("^[\\!=]=", "", by)
      x  <- sub("^[\\!=]=", "", x)
      
    }
    
  } else if (is(biom, 'data.frame')) {
    
    #--------------------------------------------------------------
    # The user supplied the data frame to run stats on
    #--------------------------------------------------------------
    
    df <- biom
    
    missing <- setdiff(c(x, y, by), names(df))
    if (length(missing) > 0)
      stop(
        "Expected column not present in data frame: ", 
        paste(collapse = ", ", missing) )
    
  } else {
    
    #--------------------------------------------------------------
    # The user gave us an unexpected data type to run stats on
    #--------------------------------------------------------------
    
    stop(
      "'biom' should be BIOM or data.frame, not ", 
      paste(collapse=", ", class(biom)) )
  }
  
  
  #--------------------------------------------------------------
  # Run non-parametric statistics
  #--------------------------------------------------------------
  plyby <- NULL
  for (i in rev(by)) plyby <- c(plyr::as.quoted(as.name(i)), plyby)
  results <- plyr::ddply(df, plyby, function (z) {
    
    if (isTRUE(pairwise)) {
      
      #--------------------------------------------------------------
      # wilcox-style stats
      #--------------------------------------------------------------
      groups <- if (is.factor(z[[x]])) levels(z[[x]]) else sort(unique(z[[x]]))
      
      if (length(groups) < 2) {
        
        result <- data.frame(
          check.names = FALSE, 
          n      = nrow(z),
          Test   = NA,
          Group1 = NA,
          Group2 = NA,
          p.val  = 1 )
        
      } else {
        
        pairs <- combn(groups, 2)
        
        result <- data.frame(
          check.names = FALSE,
          n      = apply(pairs, 2L, function (pair) {
            sum(z[[x]] %in% pair)
          }),
          Test   = "Mann-Whitney",
          Group1 = pairs[1,],
          Group2 = pairs[2,],
          p.val  = apply(pairs, 2L, function (pair) {
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
      #--------------------------------------------------------------
      # kruskal-style stats, possibly interspersed with wilcox
      #--------------------------------------------------------------
      
      ngroups <- length(unique(z[[x]]))
      result  <- data.frame(n = nrow(z))
      
      if (ngroups < 2) {
        result[['Test']]  <- NA
        result[['p.val']] <- 1
        
      } else if (ngroups == 2) {
        result[['Test']]  <- "Mann-Whitney"
        result[['p.val']] <- tryCatch({
          groups <- as.character(unique(z[[x]]))
          wt <- wilcox.test(
            x = z[z[[x]] == groups[[1]], y],
            y = z[z[[x]] == groups[[2]], y] ) %>%
            suppressWarnings()
          wt[['p.value']]
        }, error = function (e) 1)
        
      } else {
        result[['Test']]  <- "Kruskal-Wallis"
        result[['p.val']] <- tryCatch({
          kruskal.test(z[[y]], z[[x]])[['p.value']]
        }, error = function (e) 1)
        
      }
    }
    
    
    if (isTRUE(any(y.pos %in% c("max", "box", "violin")))) {
      if (y.pos == "max")    result[['y.pos']] <- max(z[[y]][is.finite(z[[y]])])
      if (y.pos == "box")    result[['y.pos']] <- IQR(z[[y]]) * 1.5
      if (y.pos == "violin") result[['y.pos']] <- max(density(z[[y]])[['x']])
    }
    
    return (result)
    
  })
  
  # Remove '.id' column that ddply adds when by=NULL
  if (length(by) == 0)
    results <- results[,-1,drop=F]
  
  results[['adj.p']] <- p.adjust(results[['p.val']], method = adj)
  
  # Order by p-value, most to least significant (unless plotting)
  if (is.null(y.pos))
    results <- results[order(results[['p.val']]),,drop=F]
  
  # Round the p.val and adj.p numbers
  if (!is.null(digits)) {
    results[['p.val']] <- signif(results[['p.val']], digits = digits)
    results[['adj.p']] <- signif(results[['adj.p']], digits = digits)
  }
  
  # Make y.pos the last column
  if (!is.null(y.pos) && "y.pos" %in% names(results)) {
    
    y.pos.facet <- intersect(y.pos.facet, by)
    
    if (is.null(y.pos.facet)) {
      
      results[['y.pos']] <- signif(max(results[['y.pos']]), digits = 2)
      
    } else {
      plyby <- NULL
      for (i in rev(y.pos.facet)) plyby <- c(plyr::as.quoted(as.name(i)), plyby)
      results <- plyr::ddply(results, plyby, function (z) {
        z[['y.pos']] <- signif(max(z[['y.pos']]), digits = 2)
        return (z)
      })
    }
    
    results <- results[,order(names(results) == "y.pos"),drop=F]
  }
  
  rownames(results) <- NULL
  
  return (results)
}
