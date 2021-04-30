#' Make a distance matrix of samples vs samples.
#' 
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'     object, as returned from \link{read.biom}. For matrices, the rows and 
#'     columns are assumed to be the taxa and samples, respectively.
#' @param method  The distance algorithm to use. Options are:
#'     \bold{\dQuote{manhattan}}, \bold{\dQuote{euclidean}}, 
#'     \bold{\dQuote{bray-curtis}}, \bold{\dQuote{jaccard}}, and
#'     \bold{\dQuote{unifrac}}. Non-ambiguous abbreviations of the method 
#'     names are also accepted. A phylogenetic tree must be present in 
#'     \code{biom} or explicitly provided via \code{tree=} to use the UniFrac methods.
#' @param weighted  Take relative abundances into account. When 
#'     \code{weighted=FALSE}, only presence/absence is considered.
#' @param tree  A \code{phylo} object representing the phylogenetic
#'     relationships of the taxa in \code{biom}. Will be taken from the tree
#'     embedded in the \code{biom} object if not explicitly specified. Only
#'     required for computing UniFrac distance matrices.
#' @param long  Pivot the returned data to long format?
#'     \describe{
#'       \item{\bold{FALSE}}{ Return a \bold{dist} object. (Default) }
#'       \item{\bold{TRUE}}{  A \bold{data.frame} with columns "Sample1", 
#'                            "Sample2" and "Distance" is returned. }
#'     }
#' @param md  Include metadata in the output data frame? Options are: 
#'     \describe{
#'       \item{\code{FALSE}}{ Don't include metadata. (Default) }
#'       \item{\code{TRUE}}{ Include all metadata. }
#'       \item{\emph{character vector}}{ Include only the specified metadata columns.
#'              Column names can be prefixed with \bold{==} or \bold{!=} to indicate
#'              that only within or between groupings, respectively, are to be kept.
#'              See examples below. }
#'     }
#' @return If both \code{long} and \code{md} are \bold{FALSE}, returns a distance 
#'         matrix, otherwise a data.frame.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     fullbiom <- read.biom(infile)
#'     biom <- select(fullbiom, 1:10)
#'     
#'     dm <- beta.div(biom, 'unifrac')
#'     
#'     as.matrix(dm)[1:4,1:4]
#'     plot(hclust(dm))
#'     
#'     # Return in long format with metadata
#'     biom <- select(fullbiom, 18:21)
#'     beta.div(biom, 'unifrac', md = c("Body Site", "Sex"))
#'     
#'     # Only look at distances among the stool sample
#'     beta.div(biom, 'unifrac', md = c("==Body Site", "Sex"))
#'     
#'     # Or between males and females
#'     beta.div(biom, 'unifrac', md = c("Body Site", "!=Sex"))
#'

beta.div <- function (biom, method, weighted=TRUE, tree=NULL, long=FALSE, md=FALSE) {
  
  #--------------------------------------------------------------
  # Enable abbreviations of metric names.
  #--------------------------------------------------------------
  
  methodList <- c("manhattan", "euclidean", "bray-curtis", "jaccard", "unifrac")
  method <- methodList[pmatch(tolower(method), methodList)]
  
  
  #--------------------------------------------------------------
  # Sanity Checks
  #--------------------------------------------------------------
  
  if (!is.logical(weighted)) stop(simpleError("Weighted must be TRUE/FALSE."))
  if (length(method) != 1)   stop(simpleError("Invalid method for beta.div()"))
  if (is.na(method))         stop(simpleError("Invalid method for beta.div()"))
  
  
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
  # Find the UniFrac tree
  #--------------------------------------------------------------
  
  if (identical(method, "unifrac")) {
    if (!is(tree, "phylo")) {
      if (is(biom, "BIOM")) {
        if (is(biom$phylogeny, "phylo")) {
          tree <- biom$phylogeny
        }
      }
      if (is(tree, "character")) {
        if (file.exists(tree)) {
          tree <- rbiom::read.tree(tree)
        }
      }
      if (!is(tree, "phylo")) {
        stop(simpleError("No tree provided to beta.div()."))
      }
    }
    
    if (length(setdiff(rownames(counts), tree$tip.label)) > 0)
      stop(simpleError("OTUs missing from reference tree."))
    
    if (length(setdiff(tree$tip.label, rownames(counts))) > 0)
      tree <- rbiom::subtree(tree, rownames(counts))
    
    counts <- counts[as.character(tree$tip.label),]
  }
  
  
  #--------------------------------------------------------------
  # Order the sparse matrix's values by sample, then by taxa
  #--------------------------------------------------------------
  
  ord      <- order(counts$j, counts$i)
  counts$i <- counts$i[ord]
  counts$j <- counts$j[ord]
  counts$v <- counts$v[ord]
  
  
  #--------------------------------------------------------------
  # Run C++ implemented dissimilarity algorithms multithreaded
  #--------------------------------------------------------------
  
  if (identical(method, "unifrac")) {
    
    dm <- par_unifrac(counts, tree, ifelse(weighted, 1L, 0L))
    
    
  } else {
    
    counts <- t(as.matrix(counts))
    dm <- par_beta_div(counts, method, ifelse(weighted, 1L, 0L))
    dm <- as.dist(dm)
    attr(dm, 'Labels') <- rownames(counts)
    
  }
  
  if (isFALSE(long) && isFALSE(md))
    return (dm)
  
  
  #--------------------------------------------------------------
  # Convert to long form
  #--------------------------------------------------------------
  df <- as.matrix(dm)
  df <- data.frame(
    stringsAsFactors = FALSE,
    Sample1  = rownames(df)[row(df)],
    Sample2  = colnames(df)[col(df)],
    Distance = as.numeric(df)
  )
  df <- subset(df, Sample1 < Sample2)
  
  
  #--------------------------------------------------------------
  # Add metadata columns
  #--------------------------------------------------------------
  if (!isFALSE(md)) {
    
    if (isTRUE(md))        md <- names(metadata(biom))
    if (!is.character(md)) md <- names(metadata(biom))[md]
    
    for (i in unique(md)) {
      
      op <- substr(i, 1, 2)
      if (op %in% c("==", "!=")) {
        
        # Limit to only within or between comparisons.
        #--------------------------------------------------------------
        i   <- substr(i, 3, nchar(i))
        map <- metadata(biom, i)
        df  <- df[get(op)(map[df$Sample1], map[df$Sample2]),,drop=F]
        
      } else {
        map <- metadata(biom, i)
      }
      
      v1 <- as.character(map[df$Sample1])
      v2 <- as.character(map[df$Sample2])
      
      
      # Change "Male vs Female" to "Female vs Male" (alphabetical).
      #--------------------------------------------------------------
      df[[i]] <- paste(
        ifelse(v1 < v2, v1, v2), 
        "vs", 
        ifelse(v1 < v2, v2, v1) )
      
      
      # Change "Male vs Male" to "Male".
      #--------------------------------------------------------------
      df[[i]] <- ifelse(v1 == v2, v1, df[[i]])
      
    }
  }
  
  return (df)
  
}





