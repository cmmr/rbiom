#' Reduce a distance matrix to two or three dimensions
#' 
#' @name ordinate
#' @param biom    A \code{BIOM} object, as returned from \link{read.biom}. Alternatively,
#'                a distance matrix, such as from \link{beta.div}.
#' @param ord     Options are:
#'     \describe{
#'         \item{pcoa}{
#'           Principal coordinate analysis, via \code{\link[ape]{pcoa}}.
#'         }
#'         \item{nmds}{
#'           Nonmetric multidimensional scaling, via \code{\link[vegan]{metaMDS}}.
#'         }
#'         \item{tsne}{
#'           t-distributed stochastic neighbor embedding, via \code{\link[tsne]{tsne}}.
#'         }
#'        }
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
#' @param md  Include metadata in the output data frame? Options are: 
#'     \describe{
#'       \item{\code{FALSE}}{ Don't include metadata. (Default) }
#'       \item{\code{TRUE}}{ Include all metadata. }
#'       \item{\emph{character vector}}{ Include only the specified metadata columns. }
#'     }
#' @param k  Number of dimensions to return. (Default: 2)
#' @param ...  Additional arguments to pass on to \code{\link[ape]{pcoa}}, 
#'             \code{\link[vegan]{metaMDS}}, or \code{\link[tsne]{tsne}}.
#' @return A data.frame with columns Axis.1, Axis.2, ..., Axis.k as well as columns
#'         given by \code{md}. Rownames are the sample IDs. For pcoa ordinations, 
#'         \code{attr(, 'eig')} will contain the eigenvalues useful for construction 
#'         of ".. % variation explained" labels.
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom   <- read.biom(infile)
#'     ord    <- ordinate(biom, "pcoa", "bray-curtis")
#'     head(ord)
#'     

ordinate <- function (biom, ord, method, weighted=TRUE, tree=NULL, md=FALSE, k=2, ...) {
  
  dm <- beta.div(biom = biom, method = method, weighted = weighted, tree = tree)
  
  if (ord == "pcoa") {
    res <- ape::pcoa(dm, ...)
    ord <- res$vectors[,1:k]
    attr(ord, 'eig') <- res$values$Relative_eig[1:k]
    
  } else if (ord == "tsne") {
    ord <- suppressMessages(tsne::tsne(dm, k=k, ...))
    rownames(ord) <- attr(dm, "Labels")
    
  } else if (ord == "nmds") {
    res <- vegan::metaMDS(dm, k=k, trace=0, ...)
    ord <- res$points
    
  } else {
    stop("'", ord, "' is not a valid argument for 'ord'.")
  }
  colnames(ord) <- paste0("Axis.", 1:k)
  ord <- as.data.frame(ord)
    
  if (isTRUE(md)) md <- colnames(metadata(biom))
  if (is.character(md) && length(md) > 0)
    for (i in md)
      ord[[i]] <- metadata(biom, i)[rownames(ord)]
  
    
  return (ord)
  
}
