
ord_metrics  <- function () c("pcoa", "tsne", "nmds")
adiv_metrics <- function () c("OTUs", "Shannon", "Chao1", "Simpson", "InvSimpson")
bdiv_metrics <- function () c("manhattan", "euclidean", "bray-curtis", "jaccard", "unifrac")

validate_metrics <- function (biom, metrics) {
  
  taxa <- c('Abundance', taxa.ranks(biom), 'OTU')
  adiv <- c('Diversity', adiv_metrics())
  bdiv <- c('Distance',  bdiv_metrics())
  opts <- c(taxa, adiv, bdiv)
  okay <- pmatch(tolower(metrics), tolower(opts))
  vals <- opts[okay]
  
  missing <- which(is.na(okay))
  if (length(missing) > 0)
    stop("Invalid metric(s): ", paste(collapse = ", ", metrics[missing]))
  
  if (length(vals) > 1 && !all(vals %in% adiv))
    stop("Invalid combination of metrics: ", paste(collapse = ", ", vals))
  
  if ("Diversity" %in% vals) vals <- adiv_metrics()
  if ("Abundance" %in% vals) vals <- tail(c('OTU', taxa.ranks(biom)), 1)
  if ("Distance"  %in% vals) vals <- ifelse(is(phylogeny(biom), 'phylo'), 'unifrac', 'bray-curtis')
  
  if (all(vals %in% taxa)) attr(vals, 'mode') <- "taxa"
  if (all(vals %in% adiv)) attr(vals, 'mode') <- "adiv"
  if (all(vals %in% bdiv)) attr(vals, 'mode') <- "bdiv"
  
  return (vals)
}


validate_md_cols <- function (biom, md_cols) {
  
  ignore  <- c(".", "pcoa", "tsne", "nmds")
  md_cols <- md_cols[!md_cols %in% ignore]
  
  missing <- setdiff(sub("^[!=]=", "", md_cols), names(metadata(biom)))
  if (length(missing) > 0)
    stop("Invalid metadata columns: ", paste(collapse = ", ", missing))
  
  return (md_cols)
}

