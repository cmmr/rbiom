

#' Cluster samples by taxa abundances (k-means).
#'
#' Runs [stats::kmeans()] on the samples x taxa abundance matrix, with
#' `nstart = 10` random starts so that the partition is stable, and with the
#' random number generator seeded by `seed` so that it is reproducible.
#'
#' @inherit documentation_clusters
#' @inherit documentation_default
#'
#' @family taxa_abundance
#' @family clustering
#'
#' @param seed   Random seed for the k-means starting centers.
#'        Default: `0`
#'
#' @param ... Passed on to [stats::kmeans()].
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

taxa_clusters <- function (biom, rank = ".otu", k = 5, seed = 0, ...) {

  biom <- as_rbiom(biom)
  schema_check('taxa_clusters', environment())
  validate_seed()

  mtx  <- taxa_matrix(biom = biom, rank = rank)
  args <- c(list(x = t(mtx), centers = k), list(...))
  if (is.null(args$nstart)) args$nstart <- 10

  km <- with_seed(seed, do.call(stats::kmeans, args))

  factor(km$cluster)
}


#' Cluster samples by beta diversity (k-medoids).
#'
#' Partitions the samples around `k` medoids with [cluster::pam()], which
#' works directly on the dissimilarity matrix. (k-means is not a valid
#' clustering of a dissimilarity matrix, since it would treat the rows of
#' the matrix as coordinates.)
#'
#' @inherit documentation_clusters
#' @inherit documentation_default
#'
#' @family beta_diversity
#' @family clustering
#'
#' @param seed   Random seed; only used when [cluster::pam()] is asked to
#'        sample medoids (`medoids = "random"`). Default: `0`
#'
#' @param ... Passed on to [cluster::pam()].
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

bdiv_clusters <- function (
    biom, bdiv = "bray", weighted = NULL, normalized = NULL,
    tree = NULL, k = 5, seed = 0, alpha = 0.5, cpus = n_cpus(), ...) {

  biom <- as_rbiom(biom)
  schema_check('bdiv_clusters', environment())
  validate_bdiv()
  validate_seed()

  dm <- bdiv_distmat(
    biom  = biom,
    bdiv  = bdiv,
    tree  = tree,
    alpha = alpha,
    cpus  = cpus )

  validate_finite_dist(dm, bdiv, "clustering")

  pm <- with_seed(seed, cluster::pam(dm, k = k, ...))

  factor(pm$clustering)
}
