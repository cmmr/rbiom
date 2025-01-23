

#' Cluster samples by taxa abundances k-means.
#' 
#' @inherit documentation_clusters
#' @inherit documentation_default
#' 
#' @family taxa_abundance
#' @family clustering
#' 
#' @param ... Passed on to `stats::kmeans()`.
#' 
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     biom$metadata$otu_cluster <- taxa_clusters(biom)
#'     
#'     pull(biom, 'otu_cluster')[1:10]
#'     
#'     bdiv_ord_plot(biom, layers = "p", stat.by = "otu_cluster")

taxa_clusters <- function (biom, rank = ".otu", k = 5, ...) {
  mtx <- taxa_matrix(biom = biom, rank = rank)
  factor(stats::kmeans(t(mtx), centers = k, ...)$cluster)
}


#' Cluster samples by beta diversity k-means.
#' 
#' @inherit documentation_clusters
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' @family clustering
#' 
#' @param ... Passed on to `stats::kmeans()`.
#' 
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     biom$metadata$bray_cluster <- bdiv_clusters(biom)
#'     
#'     pull(biom, 'bray_cluster')[1:10]
#'     
#'     bdiv_ord_plot(biom, stat.by = "bray_cluster")

bdiv_clusters <- function (biom, bdiv = "Bray-Curtis", weighted = TRUE, normalized = TRUE, tree = NULL, k = 5, ...) {
  dm <- bdiv_distmat(biom = biom, bdiv = bdiv, weighted = weighted, normalized = TRUE, tree = tree)
  factor(stats::kmeans(dm, centers = k, ...)$cluster)
}
