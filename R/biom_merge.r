
#' Combine several rbiom objects into one.
#' 
#' WARNING: It is generally ill-advised to merge BIOM datasets, as OTUs
#' mappings are dependent on upstream clustering and are not equivalent
#' between BIOM files.
#' 
#' @inherit documentation_return.biom return
#' 
#' @family biom
#' 
#' @param ...  Any number of rbiom objects (e.g. from [as_rbiom()]), lists of
#'        rbiom objects, or valid arguments to the \code{biom} parameter of 
#'        [as_rbiom()] (for instance file names).
#'        
#' @param metadata,taxonomy,tree,sequences,id,comment Replace the corresponding 
#'        data in the merged rbiom object with these values. Set to `NULL` to
#'        not inherit a particular component. The default, `NA`, will attempt
#'        to create the component based on `...` values. The merged 
#'        phylogenetic tree cannot be inferred.
#'        
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     b1 <- as_rbiom(hmp50$counts[,1:4])
#'     b2 <- as_rbiom(hmp50$counts[,5:8])
#'     
#'     biom <- biom_merge(b1, b2)
#'     print(biom)
#'     
#'     biom$tree     <- hmp50$tree
#'     biom$metadata <- hmp50$metadata
#'     print(biom)

biom_merge <- function (..., metadata = NA, taxonomy = NA, tree = NULL, sequences = NA, id = NA, comment = NA) { 
  
  # --- 1. Flatten and Validate Inputs ---
  dots      <- list(...)
  biom_list <- lapply(dots, function (dot) {
    if (inherits(dot, 'rbiom'))          return (dot)
    if (length(dot) > 1)                 return (do.call(biom_merge, dot))
    if (length(dot) < 1)                 return (NULL)
    if (inherits(dot[[1]], 'rbiom'))     return (dot[[1]])
    if (rlang::is_scalar_character(dot)) return (read_biom(src = dot))
    cli::cli_abort("Invalid argument to biom_merge(): {.typeof {dot}}")
  })
  
  biom_list <- biom_list[sapply(biom_list, inherits, 'rbiom')]
  
  if (length(biom_list) == 0) cli::cli_abort("No BIOM datasets provided to biom_merge().")
  if (length(biom_list) == 1) return (biom_list[[1]])
  
  # --- 2. Determine Global Dimensions ---
  # Extract all count matrices
  counts_list <- lapply(biom_list, `[[`, 'counts')
  
  # Get all unique OTU names (rows) and sample names (columns)
  all_otus    <- unique(unlist(lapply(counts_list, rownames), use.names = FALSE))
  all_samples <- unlist(lapply(counts_list, colnames), use.names = FALSE)
  
  # Validation: Check for duplicate samples
  if (any(duplicated(all_samples))) {
    dups <- unique(all_samples[duplicated(all_samples)])
    cli::cli_abort("Sample names are not unique among BIOM datasets: {.val {dups}}")
  }
  
  # Validation: Warn if no OTU overlap
  total_n_otus <- sum(vapply(counts_list, nrow, integer(1)))
  if (length(all_otus) == total_n_otus) {
    cli::cli_warn("No overlapping OTU names. Likely incompatible datasets.")
  }
  
  
  # --- 3. Efficient Matrix Merge ---
  # Calculate column offsets for each matrix (0, ncol(m1), ncol(m1)+ncol(m2), ...)
  # This allows us to shift 'j' indices without string matching.
  nsamples_per_mx <- vapply(counts_list, ncol, integer(1))
  col_offsets     <- c(0, cumsum(nsamples_per_mx)[-length(nsamples_per_mx)])
  
  # Extract triplets (i, j, x) and remap indices
  triplets <- lapply(seq_along(counts_list), function(k) {
    mtx <- counts_list[[k]]
    sm  <- Matrix::summary(mtx) # returns 1-based i, j, x
    
    # Map local row indices to global row indices
    # map[local_idx] -> global_idx
    row_map <- match(rownames(mtx), all_otus)
    new_i   <- row_map[sm$i]
    
    # Shift column indices based on previous matrices' widths
    new_j   <- sm$j + col_offsets[k]
    
    list(i = new_i, j = new_j, x = sm$x)
  })
  
  # Reconstruct the global sparse matrix
  counts <- Matrix::sparseMatrix(
    i        = unlist(lapply(triplets, `[[`, "i"), use.names = FALSE),
    j        = unlist(lapply(triplets, `[[`, "j"), use.names = FALSE),
    x        = unlist(lapply(triplets, `[[`, "x"), use.names = FALSE),
    dims     = c(length(all_otus), length(all_samples)),
    dimnames = list(all_otus, all_samples)
  )
  
  
  # --- 4. Merge Metadata and Attributes ---
  if (is.na(metadata))  metadata  <- dplyr::bind_rows(lapply(biom_list, `[[`, 'metadata'))
  if (is.na(taxonomy))  taxonomy  <- dplyr::bind_rows(lapply(biom_list, `[[`, 'taxonomy'))
  if (is.na(sequences)) sequences <- do.call(c,       lapply(biom_list, `[[`, 'sequences'))
  
  if (is.na(id)) {
    id <- unique(unlist(lapply(biom_list, `[[`, 'id')))
    id <- setdiff(id, 'Untitled Dataset')
    if (length(id) != 1) id <- "Merged BIOM"
  }
  
  if (is.na(comment)) {
    comment <- unique(unlist(lapply(biom_list, `[[`, 'comment')))
    comment <- setdiff(comment, '')
    if (length(comment) != 1) comment <- NULL
  }
  
  rbiom$new(
    id        = id,
    comment   = comment,
    counts    = counts, 
    metadata  = metadata, 
    taxonomy  = taxonomy,
    sequences = sequences,
    tree      = tree 
  )
}
