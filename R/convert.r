
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



#' Convert a BIOM object to a TreeSummarizedExperiment object.
#' 
#' Requires the Bioconductor R package 'TreeSummarizedExperiment' to be installed.
#' 
#' @inherit adiv_boxplot params
#' 
#' @return A TreeSummarizedExperiment object.
#' 
#' @export
#' @examples
#'   \dontrun{
#'     library(rbiom)
#'     
#'     biom <- sample_rarefy(hmp50)
#'     tse  <- convert_to_TreeSummarizedExperiment(biom)
#'  }

convert_to_TreeSummarizedExperiment <- function (biom) {
  
  stopifnot(is(biom, "BIOM"))
  
  if (nchar(system.file(package = "TreeSummarizedExperiment")) == 0)
    stop("Bioconductor R package 'TreeSummarizedExperiment' must be installed to use convert_to_TreeSummarizedExperiment().")
  
  TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays       = list(otu_matrix(biom)),
    rowData      = otu_taxonomy(biom),
    colData      = sample_metadata(biom),
    rowTree      = otu_tree(biom),
    referenceSeq = otu_sequences(biom) )
}



#' Convert a BIOM object to a SummarizedExperiment object.
#' 
#' Requires the Bioconductor R package 'SummarizedExperiment' to be installed.
#' 
#' @inherit adiv_boxplot params
#' 
#' @return A SummarizedExperiment object.
#' 
#' @export
#' @examples
#'   \dontrun{
#'     library(rbiom)
#'     
#'     biom <- sample_rarefy(hmp50)
#'     se   <- convert_to_SummarizedExperiment(biom)
#'  }

convert_to_SummarizedExperiment <- function (biom) {
  
  stopifnot(is(biom, "BIOM"))
  
  if (nchar(system.file(package = "SummarizedExperiment")) == 0)
    stop("Bioconductor R package 'SummarizedExperiment' must be installed to use convert_to_SummarizedExperiment().")
  
  SummarizedExperiment::SummarizedExperiment(
    assays  = list(otu_matrix(biom)),
    rowData = otu_taxonomy(biom),
    colData = sample_metadata(biom) )
}

