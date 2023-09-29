
#' Convert a phyloseq object to an rbiom BIOM object.
#' 
#' @param phy  A phyloseq object.
#' 
#' @return A BIOM object.
#' 
#' @export
#' @examples
#'   \dontrun{
#'     library(rbiom)
#'     
#'     fp   <- system.file("extdata", "rich_sparse_otu_table.biom", package="phyloseq")
#'     phy  <- phyloseq::import_biom(fp)
#'     biom <- convert_from_phyloseq(phy)
#'  }

convert_from_phyloseq <- function (phy) {
  
  stopifnot(is(phy, "phyloseq"))
  
  if (nchar(system.file(package = "phyloseq")) == 0)
    stop("Bioconductor R package 'phyloseq' must be installed to use convert_from_phyloseq().")
  
  structure(
    class = "BIOM", 
    .Data = list(
      counts    = phy %>% phyloseq::otu_table() %>% slam::as.simple_triplet_matrix(), 
      metadata  = phy %>% phyloseq::sample_data(errorIfNULL = FALSE) %>% data.frame(), 
      taxonomy  = phy %>% phyloseq::tax_table(errorIfNULL = FALSE) %>% data.frame() %>% as.matrix(), 
      sequences = phy %>% phyloseq::refseq(errorIfNULL = FALSE) %>% as.vector(), 
      phylogeny = phy %>% phyloseq::phy_tree(errorIfNULL = FALSE), 
      info      = list(
        id        = "Imported PhyloSeq Data",
        shape     = c(phyloseq::ntaxa(phy), phyloseq::nsamples(phy)) )))
  
  biom <- biom_repair(biom, prune = TRUE)
  return (biom)
}
