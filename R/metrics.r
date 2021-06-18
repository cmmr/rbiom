
# List all the options for each type of metric.

metrics <- function (mode, query = NULL, biom = NULL) {
  
  mode <- tolower(mode)
  opts <- NULL
  
  if (mode == "alpha") opts <- c('Depth', 'OTUs', 'Shannon', 'Chao1', 'Simpson', 'InvSimpson')
  if (mode == "beta")  opts <- c("Manhattan", "Euclidean", "Bray-Curtis", "Jaccard", "UniFrac")
  if (mode == "ord")   opts <- c("PCoA", "NMDS", "tSNE")
  
  if (mode == "meta" && !is.null(biom))
    opts <- names(biom[['metadata']])
  
  if (mode == "taxa") {
    opts <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
    if (!is.null(biom)) opts <- colnames(biom[['taxonomy']])
    opts <- c(opts, "OTU")
  }
  
  
  if (!is.null(opts) && !is.null(query)) {
    
    if (mode == "meta") {
      opts <- query[sub("^[!=]=", "", query) %in% opts]
    } else {
      opts <- opts[tolower(opts) %in% tolower(query)]
    }
    
    
    if (length(opts) == 0) return (NULL)
  }
  
  return (opts)
}
