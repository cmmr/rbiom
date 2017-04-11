#' Make a distance matrix of samples vs samples.
#' 
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'     object, as returned from \link{read.biom}. For matrices, the rows and 
#'     columns are assumed to be the taxa and samples, respectively.
#' @param method  The distance algorithm to use. Options are:
#'     \bold{\dQuote{manhattan}}, \bold{\\dQuote{euclidean}}, 
#'     \bold{\dQuote{bray}}, \bold{\dQuote{jaccard}}, and
#'     \bold{\dQuote{unifrac}}. Non-ambiguous abbrevations of the method 
#'     names are also accepted. A phylogentic tree must be present in 
#'     \code{biom} to use the UniFrac methods.
#' @param weighted  Take relative abundances into account. When 
#'     \code{weighted=FALSE}, only presence/absence is considered.
#' @param tree  A \code{phylo} object representing the phylogenetic
#'     relationships of the taxa in \code{biom}. Will be taken from the tree
#'     embedded in the \code{biom} object if not explicitly specified. Only
#'     required for computing UniFrac distance matrices.
#' @param progressbar  Whether to display a progress bar and status messages
#'     (TRUE/FALSE). Will automatically tie in with \pkg{shiny} if run within a
#'     \pkg{shiny} session.
#' @return A distance matrix.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     biom <- select(biom, 1:10)
#'     
#'     dm <- beta.div(biom, 'unifrac')
#'     
#'     as.matrix(dm)[1:4,1:4]
#'     plot(hclust(dm))
#'

beta.div <- function (biom, method, weighted=TRUE, tree=NULL, progressbar=FALSE) {
  
  
  #--------------------------------------------------------------
  # Enable abbreviations of metric names.
  #--------------------------------------------------------------
  
  methodList <- c("manhattan", "euclidean", "bray", "jaccard", "unifrac")
  method <- methodList[pmatch(method, methodList)]
  
  
  
  #--------------------------------------------------------------
  # Sanity Checks
  #--------------------------------------------------------------
  
  if (length(method) != 1 | identical(method, NA))
    stop(simpleError("Invalid method for beta.div()"))
  
  if (!is(biom, "BIOM"))
    stop(simpleError("beta.div function needs a BIOM object."))
  
  
  
  #--------------------------------------------------------------
  # Use stand-alone function for UniFrac
  #--------------------------------------------------------------
  
  if (method == "unifrac") {
    return (unifrac(biom, weighted=weighted, tree=tree, progressbar=progressbar))
  }
  
  
  
  #--------------------------------------------------------------
  # Jaccard is simply a transformation of bray
  #--------------------------------------------------------------
  
  if (method == "jaccard") {
    dm <- beta.div(biom, "bray", weighted, progressbar)
    dm <- 2 * dm / (1 + dm)
    return (dm)
  }
  
  
  
  #--------------------------------------------------------------
  # Get the input into a simple_triplet_matrix
  #--------------------------------------------------------------
  if (is(biom, "simple_triplet_matrix")) { counts <- biom
  } else if (is(biom, "BIOM"))           { counts <- biom$counts
  } else if (is(biom, "matrix"))         { counts <- slam::as.simple_triplet_matrix(biom)
  } else {
    stop(simpleError("biom must be a matrix, simple_triplet_matrix, or BIOM object."))
  }
  
  
  #--------------------------------------------------------------
  # Lightweight dissimilarity algorithms
  #--------------------------------------------------------------
  
  as.dist(
    slam::crossapply_simple_triplet_matrix(
      x   = counts,
      FUN = switch( as.character(weighted),
                    
        "FALSE" = switch( method,
          "manhattan"   = function(x,y) { sum(x>0)+sum(y>0)-2*sum(x&y) }, 
          "euclidean"   = function(x,y) { sqrt(sum(x>0)+sum(y>0)-2*sum(x&y)) }, 
          "bray"        = function(x,y) { (sum(x>0)+sum(y>0)-2*sum(x&y))/(sum(x>0)+sum(y>0)) }
        ),
        
        "TRUE" = switch( method,
          "manhattan"   = function(x,y) { sum(abs(x-y)) }, 
          "euclidean"   = function(x,y) { sqrt(sum((x-y)^2)) }, 
          "bray"        = function(x,y) { sum(abs(x-y))/sum(x+y) }
        )
        
      )
    )
  )
  
}





