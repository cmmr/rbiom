
#' Convert biom data to an external package class.
#'
#' Requires the relevant Bioconductor R package to be installed:
#' \describe{
#'   \item{`convert_to_animalcules` - }{ [animalcules](https://bioconductor.org/packages/animalcules/) }
#'   \item{`convert_to_biomformat` - }{ [biomformat](https://bioconductor.org/packages/biomformat/) }
#'   \item{`convert_to_phyloseq` - }{ [phyloseq](https://bioconductor.org/packages/phyloseq/) }
#'   \item{`convert_to_SE` - }{ [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/) }
#'   \item{`convert_to_TSE` - }{ [TreeSummarizedExperiment](https://bioconductor.org/packages/TreeSummarizedExperiment/) }
#' }
#'
#' animalcules, SummarizedExperiment, and biomformat objects include counts, metadata, and taxonomy.
#'
#' phyloseq and TreeSummarizedExperiment additionally include the tree and sequences.
#'
#' @name convert_to
#' @inherit documentation_default
#'
#' @param ...  Not Used.
#'
#' @return An animalcules (MultiAssayExperiment), biomformat (biom), phyloseq, SummarizedExperiment, or TreeSummarizedExperiment object.
#'
#' @export
#' @examples
#' \dontrun{
#'     library(rbiom)
#'
#'     print(hmp50)
#'
#'     # Requires 'animalcules', a Bioconductor R package
#'     if (nzchar(system.file(package = "animalcules"))) {
#'       ani <- convert_to_animalcules(hmp50)
#'       print(ani)
#'     }
#'
#'     # Requires 'biomformat', a Bioconductor R package
#'     if (nzchar(system.file(package = "biomformat"))) {
#'       bio <- convert_to_biomformat(hmp50)
#'       print(bio)
#'     }
#'
#'     # Requires 'phyloseq', a Bioconductor R package
#'     if (nzchar(system.file(package = "phyloseq"))) {
#'       phy <- convert_to_phyloseq(hmp50)
#'       print(phy)
#'     }
#'
#'     # Requires 'SummarizedExperiment', a Bioconductor R package
#'     if (nzchar(system.file(package = "SummarizedExperiment"))) {
#'       se <- convert_to_SE(hmp50)
#'       print(se)
#'     }
#'
#'     # Requires 'TreeSummarizedExperiment', a Bioconductor R package
#'     if (nzchar(system.file(package = "TreeSummarizedExperiment"))) {
#'       tse <- convert_to_TSE(hmp50)
#'       print(tse)
#'     }
#' }

convert_to_animalcules <- function (biom, ...) {
  
  require_package('S4Vectors',            'to use convert_to_animalcules()')
  require_package('SummarizedExperiment', 'to use convert_to_animalcules()')
  require_package('MultiAssayExperiment', 'to use convert_to_animalcules()')
  
  dots <- list(...)
  biom <- as_rbiom(biom)
  
  DataFrame            <- getFromNamespace('DataFrame',  'S4Vectors')
  SimpleList           <- getFromNamespace('SimpleList', 'S4Vectors')
  SummarizedExperiment <- getFromNamespace('SummarizedExperiment', 'SummarizedExperiment')
  MultiAssayExperiment <- getFromNamespace('MultiAssayExperiment', 'MultiAssayExperiment')
  
  assays <- SimpleList('MGX' = as.matrix(biom$counts))
  
  colData <- biom$metadata %>%
    tibble::column_to_rownames(".sample") %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    DataFrame()
  
  rowData <- biom$taxonomy %>%
    tibble::column_to_rownames(".otu") %>%
    dplyr::mutate_all(as.character) %>%
    DataFrame()
  
  se <- SummarizedExperiment(
    assays  = assays,
    colData = colData,
    rowData = rowData )
  
  MultiAssayExperiment(
    experiments = SimpleList('MicrobeGenetics' = se),
    colData     = colData )
}



#' @rdname convert_to
#' @export

convert_to_biomformat <- function (biom, ...) {
  
  require_package('biomformat', 'to use convert_to_biomformat()')
  
  dots <- list(...)
  biom <- as_rbiom(biom)
  
  make_biom <- getFromNamespace('make_biom', 'biomformat')
  make_biom(
    data                 = as.matrix(biom$counts),
    sample_metadata      = biom$metadata %>% tibble::column_to_rownames(".sample"),
    observation_metadata = biom$taxonomy %>% tibble::column_to_rownames(".otu"),
    id                   = biom$id )
}



#' @rdname convert_to
#' @export

convert_to_phyloseq <- function (biom, ...) {
  
  require_package('phyloseq', 'to use convert_to_phyloseq()')
  
  dots <- list(...)
  biom <- as_rbiom(biom)
  
  import_biom <- getFromNamespace('import_biom', 'phyloseq')
  
  BIOMfilename   <- write_biom(biom, file = tempfile())
  treefilename   <- if (!is.null(biom$tree))      write_tree(biom,  file = tempfile())
  refseqfilename <- if (!is.null(biom$sequences)) write_fasta(biom, file = tempfile())
  
  on.exit(unlink(c(BIOMfilename, treefilename, refseqfilename)), add = TRUE)
  
  import_biom(
    BIOMfilename   = BIOMfilename, 
    treefilename   = treefilename, 
    refseqfilename = refseqfilename )
}



#' @rdname convert_to
#' @export

convert_to_SE <- function (biom, ...) {
  
  require_package('SummarizedExperiment', 'to use convert_to_SE()')
  
  dots <- list(...)
  biom <- as_rbiom(biom)
  
  SummarizedExperiment <- getFromNamespace('SummarizedExperiment', 'SummarizedExperiment')
  SummarizedExperiment(
    assays  = list('OTU table' = as.matrix(biom$counts)),
    rowData = biom$taxonomy %>% tibble::column_to_rownames(".otu") %>% as.matrix(),
    colData = biom$metadata %>% tibble::column_to_rownames(".sample") )
  
}



#' @rdname convert_to
#' @export

convert_to_TSE <- function (biom, ...) {

  require_package('TreeSummarizedExperiment', 'to use convert_to_TSE()')

  dots <- list(...)
  biom <- as_rbiom(biom)
  
  TreeSummarizedExperiment <- getFromNamespace('TreeSummarizedExperiment', 'TreeSummarizedExperiment')
  DNAStringSet             <- getFromNamespace('DNAStringSet', 'Biostrings')
  
  TreeSummarizedExperiment(
    assays       = list('OTU table' = as.matrix(biom$counts)),
    rowData      = biom$taxonomy %>% tibble::column_to_rownames(".otu") %>% as.matrix(),
    colData      = biom$metadata %>% tibble::column_to_rownames(".sample"),
    rowTree      = biom$tree,
    referenceSeq = DNAStringSet(biom$sequences) )
  
}


