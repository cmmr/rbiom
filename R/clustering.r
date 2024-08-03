

#' Define sample kmeans clusters from taxa abundances.
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


#' Define sample PAM clusters from beta diversity.
#' 
#' @inherit documentation_clusters
#' @inherit documentation_default
#' 
#' @family beta_diversity
#' @family clustering
#' 
#' @param ... Passed on to `cluster::pam()`.
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

bdiv_clusters <- function (biom, bdiv = "Bray-Curtis", weighted = TRUE, tree = NULL, k = 5, ...) {
  dm <- bdiv_matrix(biom = biom, bdiv = bdiv, weighted = weighted, tree = tree)
  factor(cluster::pam(dm, k = k, ...)$clustering)
}
