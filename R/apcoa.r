# #' PCoA adjusted for corvariates.
# #' 
# #' Please cite: Shi Y, Zhang L, Do KA, Peterson CB, Jenq RR. aPCoA: covariate
# #' adjusted principal coordinates analysis. Bioinformatics. 2020 Jul 
# #' 1;36(13):4099-4101. doi: 10.1093/bioinformatics/btaa276.
# #' 
# #' @noRd
# #' 
# #' @param distmat  A distance matrix (\code{dist} class object) that you want
# #'      to run the aPCoA on.
# #' 
# #' @param covariates  A data.frame with the confounding covariate(s). The row 
# #'      names of this data frame should match the labels of the distance matrix.
# #' 
# #' @return A numeric matrix with the same row names as \code{covariates}, and
# #'      one column for each of the computed adjusted principal coordinates.
# #' 
# #' @examples
# #'     library(rbiom)
# #'     library(ggplot2)
# #'     
# #'     biom <- rarefy(hmp50)
# #'     
# #'     dm   <- bdiv_distmat(biom, 'unifrac')
# #'     reg_pcoa <- pcoa(dm)[['vectors']]
# #'     adj_pcoa <- apcoa(dm, biom$metadata[,'Sex'])
# #'     
# #'     ids   <- biom$samples
# #'     color <- pull(biom, 'Sex')
# #'     ggplot(mapping=aes(x=reg_pcoa[ids, 1], y=reg_pcoa[ids, 2], color=color)) + geom_point()
# #'     ggplot(mapping=aes(x=adj_pcoa[ids, 1], y=adj_pcoa[ids, 2], color=color)) + geom_point()
# 
# apcoa <- function (distmat, covariates)  {
#   
#   params <- eval_envir(environment())
#   
#   #________________________________________________________
#   # See if this result is already in the cache.
#   #________________________________________________________
#   cache_file <- get_cache_file('apcoa', params)
#   if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
#     return (readRDS(cache_file))
#   remove("params")
#   
#   
#   if (!inherits(distmat,    "dist"))       stop ("distmat must be of class 'dist'.")
#   if (!inherits(covariates, "data.frame")) stop ("covariates must be of class 'data.frame'.")
#   if (ncol(covariates) == 0)         stop ("covariates has no columns'.")
#   
#   # Rename the metadata columns to ensure formula-compatibility
#   names(covariates) <- paste0("CV", seq_len(ncol(covariates)))
#   
#   # Force the order of row names to match
#   rn_distmat    <- attr(distmat, "Labels", exact = TRUE)
#   rn_covariates <- rownames(covariates)
#   rn_intersect  <- intersect(rn_distmat, rn_covariates)
#   
#   if (length(rn_intersect) == 0)
#     stop("No row names in common between the distance matrix and covariates.")
#   
#   if (length(rn_intersect) < 3)
#     stop("Please provide at least 3 samples as input to apcoa().")
#   
#   distmat    <- stats::as.dist(as.matrix(distmat)[rn_intersect, rn_intersect])
#   covariates <- covariates[rn_intersect,,drop=FALSE]
#   
#   formula <- stats::reformulate(response = "distmat", termlabels = names(covariates))
#   
#   lhs          <- distmat
#   formula[[2]] <- NULL
#   rhs.frame    <- stats::model.frame(formula, covariates, drop.unused.levels = TRUE)
#   
#   rhs <- stats::model.matrix(formula, rhs.frame)
#   
#   grps   <- attr(rhs, "assign", exact = TRUE)
#   qrhs   <- qr(rhs)
#   rhs    <- rhs[, qrhs$pivot,  drop = FALSE]
#   rhs    <- rhs[, 1:qrhs$rank, drop = FALSE]
#   grps   <- grps[qrhs$pivot][1:qrhs$rank]
#   u.grps <- unique(grps)
#   nterms <- length(u.grps) - 1
#   if (nterms < 1) 
#     stop("right-hand-side of formula has no usable terms")
#   dmat <- as.matrix(lhs^2)
#   
#   X <- rhs
#   y <- lhs
#   
#   X <- X[rownames(dmat),]
#   X <- as.matrix(X[,-1],nrow=nrow(X))
#   
#   H <- X%*%solve(t(X)%*%X)%*%t(X)
#   
#   A <- -1/2*as.matrix(y)^2
#   J <- diag(nrow(X))-matrix(rep(1/(nrow(X)),length(A)),nrow=nrow(A))
#   E <- (diag(nrow(H))-H)%*%J%*%A%*%J%*%(diag(nrow(H))-H)
#   
#   rownames(E)      <- rownames(covariates)
#   colnames(E)      <- rownames(covariates)
#   eigenE           <- eigen(E)$vectors
#   eigenvalue       <- eigen(E)$values
#   rownames(eigenE) <- rownames(covariates)
#   
#   
#   plotMatrix <- eigenE*matrix(rep(eigenvalue^(1/2),each=nrow(eigenE)),nrow=nrow(eigenE))
#   plotMatrix <- plotMatrix[,!is.na(apply(plotMatrix,2,sum))]
#   plotMatrix
#   
#   
#   set_cache_value(cache_file, plotMatrix)
#   return (plotMatrix)
# }

