#' Convert BIOM data to an rbiom object.
#' 
#' @param biom   Object which can be coerced to an \code{rbiom}-class object.
#'        For example:
#'        \itemize{
#'          \item{\emph{file} - }{ Filepath or URL to a biom file. }
#'          \item{\emph{matrix} - }{ An abundance matrix with OTUs in rows and samples in columns. }
#'          \item{\code{phyloseq}-class object - }{ From the phyloseq Bioconductor R package. }
#'        }
#' 
#' @return An \code{rbiom}-class object.
#' 
#' @seealso [biom_build()], [read_biom()]
#' 
#' @export
#' 
as_rbiom <- function (biom) {
  
  
  #________________________________________________________
  # Already validated. Don't process further.
  #________________________________________________________
  if (is(biom, 'rbiom') && is(biom, 'environment'))
    return (biom)
  
  
  #________________________________________________________
  # Make sure this object hasn't been manually edited.
  #________________________________________________________
  if (is(biom, 'rbiom'))
    return (biom_repair(biom = biom))
  
  
  
  #________________________________________________________
  # Read new rbiom object from filename / URL / JSON.
  #________________________________________________________
  if (is_scalar_character(biom))
    return (read_biom(src = biom))
  
  
  
  #________________________________________________________
  # Create an rbiom object from just count data.
  #________________________________________________________
  if (is(biom, "matrix") || is(biom, "simple_triplet_matrix"))
    return (biom_build(counts = biom))
  
  
  
  #________________________________________________________
  # Allow passing phyloseq objects to rbiom functions.
  #________________________________________________________
  if (is(biom, "phyloseq"))
    return (convert_from_phyloseq(phy = biom))
  
  
  
  stop("Cannot parse `biom` argument with class '", paste(class(biom), collapse = "; "), "'.")
  
}




#' Run after manually editing an rbiom object's content.
#'
#' @param biom  The \code{rbiom} object to repair.
#'
#' @return An \code{rbiom} object.
#' @export
#' 

biom_repair <- function (biom) {
  
  
  #________________________________________________________
  # Sanity Pre-checks
  #________________________________________________________
  stopifnot(is(biom, 'rbiom'))
  stopifnot(is(biom, 'list') || is(biom, 'environment'))
  stopifnot(hasName(biom, 'counts'))
  
  
  #________________________________________________________
  # Check / sanitize contents of the rbiom object.
  #________________________________________________________
  if (!hasName(biom, 'info'))     biom[['info']]     <- list()
  if (!hasName(biom, 'metadata')) biom[['metadata']] <- tibble()
  if (!hasName(biom, 'taxonomy')) biom[['taxonomy']] <- tibble()
  
  
  
  
  #________________________________________________________
  # Check the counts matrix. Get sample and otu names.
  #________________________________________________________
  
  if (!is.simple_triplet_matrix(biom[['counts']]))
    biom[['counts']] <- tryCatch(
      expr  = as.simple_triplet_matrix(biom[['counts']]),
      error = function (e) stop ("Can't convert `counts` to matrix.\n", e) )
  
  if (!is.numeric(biom[['counts']][['v']]))     stop("The `counts` matrix must be numeric.")
  if (!all(is.finite(biom[['counts']][['v']]))) stop("Non-finite values in `counts` matrix.")
  
  sns <- colnames(biom[['counts']])
  ons <- rownames(biom[['counts']])
  
  if (is.null(sns)) stop ("`counts` must have sample names as column names.")
  if (is.null(ons)) stop ("`counts` must have OTU names as row names.")
  
  
  
  
  #________________________________________________________
  # Minimal metadata and taxonomy objects.
  #________________________________________________________
  
  if (!is_tibble(biom[['metadata']]))  biom[['metadata']] %<>% as_tibble(rownames = ".sample")
  if (!is_tibble(biom[['taxonomy']]))  biom[['taxonomy']] %<>% as_tibble(rownames = ".otu")
  if (empty(biom[['metadata']])) biom[['metadata']] <- tibble(.sample = sns)
  if (empty(biom[['taxonomy']])) biom[['taxonomy']] <- tibble(.otu    = ons)
  
  colnames(biom[['metadata']]) %<>% trimws()
  colnames(biom[['taxonomy']]) %<>% trimws()
  biom[['metadata']] <- biom[['metadata']][,!duplicated(colnames(biom[['metadata']]))]
  biom[['taxonomy']] <- biom[['taxonomy']][,!duplicated(colnames(biom[['taxonomy']]))]
  
  stopifnot(hasName(biom[['metadata']], '.sample'))
  stopifnot(hasName(biom[['taxonomy']], '.otu'))
  stopifnot(!duplicated(biom[['metadata']][['.sample']]))
  stopifnot(!duplicated(biom[['taxonomy']][['.otu']]))
  
  
  
  #________________________________________________________
  # Ensure appropriate column types and names.
  #________________________________________________________
  
  # Convert all character cols to factor
  biom[['metadata']] %<>% relocate(.sample) %>%
    mutate(across(where(is.character), as.factor))
  
  # Ensure all columns are factors
  biom[['taxonomy']] %<>% relocate(.otu) %>%
    mutate(across(!where(is.factor), as.factor))
  
  # No column names can start with a '.' (except id column)
  stopifnot(!startsWith(setdiff(colnames(biom[['metadata']]), '.sample'), '.'))
  stopifnot(!startsWith(setdiff(colnames(biom[['taxonomy']]), '.otu'),    '.'))
  
  sns %<>% intersect(as.character(biom[['metadata']][['.sample']]))
  ons %<>% intersect(as.character(biom[['taxonomy']][['.otu']]))
  
  
  
  
  #________________________________________________________
  # Drop OTUs that are absent in optional tree / seqs.
  #________________________________________________________
  
  if (!is_null(biom[['phylogeny']])) {
    stopifnot(is(biom[['phylogeny']], 'phylo'))
    ons %<>% intersect(tree_tips(biom[['phylogeny']]))
  }
  
  if (!is_null(biom[['sequences']])) {
    stopifnot(is_character(biom[['sequences']]))
    stopifnot(!is.null(names(biom[['sequences']])))
    ons %<>% intersect(names(biom[['sequences']]))
  }
  
  
  
  #________________________________________________________
  # Only keep samples / otus named in all components.
  #________________________________________________________
  
  if (!eq(ons, rownames(biom[['counts']])))
    biom[['counts']] <- biom[['counts']][ons,]
  
  if (!eq(sns, colnames(biom[['counts']])))
    biom[['counts']] <- biom[['counts']][,sns]
  
  
  
  #________________________________________________________
  # Drop taxa/samples with zero observations
  #________________________________________________________
  
  if (any(row_sums(biom[['counts']]) == 0)) {
    biom[['counts']] <- biom[['counts']][row_sums(biom[['counts']]) > 0,]
    ons <- rownames(biom[['counts']])
  }
  
  if (any(col_sums(biom[['counts']]) == 0)) {
    biom[['counts']] <- biom[['counts']][,col_sums(biom[['counts']]) > 0]
    sns <- colnames(biom[['counts']])
  }
  
  
  
  #________________________________________________________
  # Sanity Post-checks
  #________________________________________________________
  
  if (length(ons) == 0) stop ("All OTUs have been dropped.")
  if (length(sns) == 0) stop ("All samples have been dropped.")
  
  
  
  #________________________________________________________
  # Subset and order all components to match.
  #________________________________________________________
  if (!eq(sns, as.character(biom[['metadata']][['.sample']])))
    biom[['metadata']] <- left_join(
      x  = tibble(.sample = sns), 
      y  = biom[['metadata']], 
      by = ".sample" )
  
  if (!is.factor(biom[['metadata']][['.sample']]))
    biom[['metadata']][['.sample']] %<>% as.factor()
  
  if (!eq(ons, as.character(biom[['taxonomy']][['.otu']])))
    biom[['taxonomy']] <- left_join(
      x  = tibble(.otu = ons), 
      y  = biom[['taxonomy']], 
      by = ".otu" )
  
  if (!is.factor(biom[['taxonomy']][['.otu']]))
    biom[['taxonomy']][['.otu']] %<>% as.factor()
  
  if (!is.null(biom[['phylogeny']]) && !all(biom[['phylogeny']][['tip.label']] %in% ons))
    biom[['phylogeny']] %<>% tree_subset(tips = ons)
  
  if (!is.null(biom[['sequences']]) && !eq(names(biom[['sequences']]), ons))
    biom[['sequences']] <- biom[['sequences']][ons]
  
  
  
  #________________________________________________________
  # Drop missing factor levels from metadata
  #________________________________________________________
  for (i in seq_along(biom[['metadata']]))
    if (is.factor(biom[['metadata']][[i]]))
      if (!all(levels(biom[['metadata']][[i]]) %in% biom[['metadata']][[i]]))
        biom[['metadata']][[i]] %<>% {factor(., levels = intersect(levels(.), .))}
  remove("i")
  
  
  
  #________________________________________________________
  # Sanitize `info` values.
  #________________________________________________________
  
  biom[['info']] %<>% within({
    
    if (!exists('id'))      id      <- ""
    if (!exists('comment')) comment <- ""
    
    if (!is_scalar_character(id)      || is_na(id))      id      <- ""
    if (!is_scalar_character(comment) || is_na(comment)) comment <- ""
    
  })
  
  
  
  return (biom)
}

