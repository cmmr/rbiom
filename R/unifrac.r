#' Compute Weighted and Unweighted UniFrac distance matrices.
#' 
#' This is the function called internally by \link{beta.div}, but is made
#' visible here so you can use it with matrices and trees without having to
#' first convert them to \code{BIOM} objects.
#'
#' @param biom  A \code{matrix}, \code{simple_triplet_matrix}, or \code{BIOM} 
#'     object, as returned from \link{read.biom}. For matrices, the rows and 
#'     columns are assumed to be the taxa and samples, respectively.
#' @param weighted  Use weighted UniFrac, which takes abundance into account
#'     rather than simply presence/absence.
#' @param tree  A \code{phylo} object providing a phylogenetic tree for the
#'     taxa names in \code{biom}. If \code{tree=NULL}, then the tree will be
#'     loaded from \code{biom}, if encoded there.
#' @param progressbar  Whether to display a progress bar and status messages
#'     (logical). Will automatically tie in with \pkg{shiny} if run within a
#'     \pkg{shiny} session.
#' @return A distance matrix of class \code{dist}.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.biom", package = "rbiom")
#'     biom <- read.biom(infile)
#'     biom <- select(biom, 1:10)
#'     
#'     dm <- unifrac(biom)
#'     plot(hclust(dm), cex=.8)
#'     as.matrix(dm)[1:4,1:4]
#'     
#'     # Using a custom matrix and tree
#'     mtx <- matrix(sample.int(12*20), ncol=20)
#'     dimnames(mtx) <- list(LETTERS[1:12], letters[1:20])
#'     tree <- ape::as.phylo(hclust(dist(mtx)))
#'     
#'     dm <- unifrac(mtx, tree=tree)
#'     as.matrix(dm)[1:4,1:4]
#'     

unifrac <- function (biom, weighted=TRUE, tree=NULL, progressbar=FALSE) {
  
  otus <- if (is(biom, "BIOM")) biom$counts else biom
  if (is(biom, "BIOM") & is.null(tree))  tree <- biom$phylogeny
  if (is(otus, "simple_triplet_matrix")) otus <- as.matrix(otus)
  
  if (!is(otus, "matrix")) stop("No OTU table provided.")
  if (!is(tree, "phylo")) stop("No reference tree provided.")
  
  pb <- progressBar(progressbar)
  
  
  #--------------------------------------------------------------
  # Quick transformations / definitions
  #--------------------------------------------------------------
  
  otuNames    <- intersect(tree$tip.label, rownames(otus))
  otuCount    <- length(otuNames)
  sampleNames <- colnames(otus)
  sampleCount <- length(sampleNames)
  
  if (otuCount < 2 | sampleCount < 2) stop("Invalid OTU table size")
  if (otuCount < nrow(otus)) otus <- otus[rownames(otus) %in% otuNames,]
  if (weighted) otus <- otus / slam::col_sums(otus)
  
  tree <- ape::drop.tip(tree, which(!tree$tip.label %in% otuNames))
  otus <- otus[as.character(tree$tip.label),]
  
  
  
  #--------------------------------------------------------------
  # Map edges/branches to the leaves that are under them.
  #--------------------------------------------------------------
  
  edge2leaves <- list()
  otuDepths   <- setNames(rep(NA, otuCount), otuNames)
  
  edgeTraverse <- function (node, depth=0) {
    
    edgeIndex <- which(tree$edge[,2] == node)
    
    if (length(edgeIndex) == 1)
      depth  <- depth + tree$edge.length[[edgeIndex]]
    
    childNodes <- tree$edge[tree$edge[,1] == node, 2]
    if (length(childNodes) == 0) {
      leaves <- node
      otuDepths[[tree$tip.label[[tree$edge[edgeIndex, 2]]]]] <<- depth
    } else {
      leaves <- unlist(lapply(childNodes, edgeTraverse, depth=depth))
    }
    
    if (length(edgeIndex) == 1)
      edge2leaves[[edgeIndex]] <<- leaves
    
    return (invisible(leaves))
  }
  rootNode <- setdiff(unique(tree$edge[,1]), unique(tree$edge[,2]))
  edgeTraverse(rootNode)
  
  remove("rootNode")
  
  
  
  #--------------------------------------------------------------
  # Map samples to their branch weights
  #--------------------------------------------------------------
  
  branchLengths    <- setNames(tree$edge.length, as.character(1:length(tree$edge.length)))
  sample2branchwts <- sapply(sampleNames, function (sampleName) {

    res <- setNames(unlist(lapply(1:length(edge2leaves), function (edgeIndex) {
      sum(otus[edge2leaves[[edgeIndex]], sampleName])
    })), as.character(1:length(edge2leaves)))

    res <- res[res > 0]
    if (weighted == FALSE) res <- setNames(rep(1, length(res)), names(res))
    res <- res * branchLengths[names(res)]

    return (res)
  }, simplify = FALSE)
  
  
  
  #--------------------------------------------------------------
  # Distrubute workload among cpu cores
  #--------------------------------------------------------------
  
  n  <- sampleCount
  cl <- configCluster(nTasks=(n**2 - n) / 2, pb, "Running UniFrac")
  
  distResults <- {
    
    set <- i <- NULL
    
    foreach (set=cl$sets, .combine='c', .options.snow=cl$opts, .packages='foreach') %dopar% {
      foreach (i=set, .combine='c') %do% {
        
        rowIdx <- ceiling((1/2.) * (- (-8*(i-1) + 4 *n**2 -4*n - 7)**0.5 + 2*n -1) - 1) + 1
        colIdx <- n - (rowIdx * (n - 1 - rowIdx) + (rowIdx*(rowIdx + 1))/2) + i
        
        x  <- sample2branchwts[[sampleNames[[rowIdx]]]]
        y  <- sample2branchwts[[sampleNames[[colIdx]]]]
        
        if (weighted == FALSE) {
          res <- sum(x[!names(x) %in% names(y)], y[!names(y) %in% names(x)])
          res <- res / sum(x, y[!names(y) %in% names(x)])
          
        } else {
          
          z <- intersect(names(x), names(y))
          z <- abs(x[z] - y[z])
          x <- x[!names(x) %in% names(z)]
          y <- y[!names(y) %in% names(z)]
          
          res <- sum(x, y, z)
          res <- res / sum(otuDepths * otus[,rowIdx], otuDepths * otus[,colIdx])
        }
        return (res)
      }
    }
  }
  
  
  
  #--------------------------------------------------------------
  # Convert the distances into a dist object
  #--------------------------------------------------------------
  
  class(distResults)          <- "dist"
  attr(distResults, "Size")   <- n
  attr(distResults, "Upper")  <- FALSE
  attr(distResults, "Diag")   <- FALSE
  attr(distResults, "Labels") <- sampleNames
  
  
  pb$close()
  
  return (distResults)
}

