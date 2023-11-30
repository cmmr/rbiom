
#' Convert a phyloseq object to an rbiom object.
#' 
#' @family conversion
#' 
#' @param phy  A phyloseq object.
#' 
#' @return An rbiom object.
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
    class = 'rbiom', 
    .Data = list(
      counts    = phy %>% phyloseq::otu_table() %>% slam::as.simple_triplet_matrix(), 
      metadata  = phy %>% phyloseq::sample_data(errorIfNULL = FALSE) %>% data.frame(), 
      taxonomy  = phy %>% phyloseq::tax_table(errorIfNULL = FALSE) %>% data.frame(), 
      sequences = phy %>% phyloseq::refseq(errorIfNULL = FALSE) %>% as.vector(), 
      phylogeny = phy %>% phyloseq::phy_tree(errorIfNULL = FALSE), 
      info      = list(id = "Imported PhyloSeq Data") ))
  
  biom <- biom_repair(biom)
  return (biom)
}



#' Convert an rbiom object to a TreeSummarizedExperiment object.
#' 
#' Requires the Bioconductor R package 'TreeSummarizedExperiment' to be installed.
#' 
#' @inherit documentation_default
#' 
#' @family conversion
#' 
#' @return A TreeSummarizedExperiment object.
#' 
#' @export
#' @examples
#'       
#'     library(rbiom)
#'     
#'     biom <- sample_rarefy(hmp50)
#'     print(biom)
#'       
#'     # Requires 'TreeSummarizedExperiment' Bioconductor R package
#'     if (nzchar(system.file(package = "TreeSummarizedExperiment"))) {
#'       tse <- convert_to_TSE(biom)
#'       print(tse)
#'     }

convert_to_TSE <- function (biom) {
  
  validate_biom(clone = FALSE)
  
  if (nchar(system.file(package = "TreeSummarizedExperiment")) == 0)
    stop("Bioconductor R package 'TreeSummarizedExperiment' must be installed to use convert_to_TreeSummarizedExperiment().")
  
  TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays       = list('OTU table' = otu_matrix(biom)),
    rowData      = otu_taxonomy(biom)    %>% tibble::column_to_rownames(".otu") %>% as.matrix(),
    colData      = sample_metadata(biom) %>% tibble::column_to_rownames(".sample"),
    rowTree      = otu_tree(biom),
    referenceSeq = Biostrings::DNAStringSet(otu_sequences(biom)) )
}



#' Convert an rbiom object to a SummarizedExperiment object.
#' 
#' Requires the Bioconductor R package 'SummarizedExperiment' to be installed.
#' 
#' @inherit documentation_default
#' 
#' @family conversion
#' 
#' @return A SummarizedExperiment object.
#' 
#' @export
#' @examples
#'     
#'     library(rbiom)
#'     
#'     biom <- sample_rarefy(hmp50)
#'     print(biom)
#'       
#'     # Requires 'SummarizedExperiment' Bioconductor R package
#'     if (nzchar(system.file(package = "SummarizedExperiment"))) {
#'       se <- convert_to_SE(biom)
#'       print(se)
#'     }

convert_to_SE <- function (biom) {
  
  validate_biom(clone = FALSE)
  
  if (nchar(system.file(package = "SummarizedExperiment")) == 0)
    stop("Bioconductor R package 'SummarizedExperiment' must be installed to use convert_to_SummarizedExperiment().")
  
  SummarizedExperiment::SummarizedExperiment(
    assays  = list('OTU table' = otu_matrix(biom)),
    rowData = otu_taxonomy(biom)    %>% tibble::column_to_rownames(".otu") %>% as.matrix(),
    colData = sample_metadata(biom) %>% tibble::column_to_rownames(".sample") )
}

