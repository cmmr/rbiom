
#' @rdname write_biom
#' @export
write_xlsx <- function (biom, file, depth = 0.1, n = NULL, seed = 0, unc = 'singly') {
  
  require_package('openxlsx', 'to save as an Excel file')
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  unc <- match.arg(tolower(unc), c("asis", "singly", "grouped", "drop"))
  
  stopifnot(is_scalar_character(file) && !is_na(file))
  file <- normalizePath(file, winslash = "/", mustWork = FALSE)
  if (!dir.exists(dirname(file))) dir.create(dirname(file), recursive = TRUE)
  
  
  
  #________________________________________________________
  # We'll output worksheets for both full and rarefied.
  #________________________________________________________
  full <- as_rbiom(biom)
  wb   <- openxlsx::createWorkbook(creator="rbiom", title=full$id)
  
  
  #________________________________________________________
  # Bypass some sheets for non-integer data.
  #________________________________________________________
  if (all(full$counts$v %% 1 == 0)) {
    
    rare <- full$clone()
    
    if (!isTRUE(depth == 0))
      rarefy(biom = rare, depth = depth, seed = seed, clone = FALSE)
    
  } else {
    rare <- NULL
  }
  
  
  
  #________________________________________________________
  # Track "Reads Per Sample" counts before/after rarefying.
  #________________________________________________________
  if (!is.null(rare))
    local ({
      
      if ('Reads Per Step' %in% names(attributes(full))) {
        rps <- attr(full, 'Reads Per Step', exact = TRUE)
      } else {
        rps <- data.frame(Mapped=col_sums(full$counts))
      }
      
      rps <- data.frame(rps[order(rownames(rps)),,drop=FALSE])
      rps[['Rarefied']] <- col_sums(rare$counts)[rownames(rps)]
      rps[which(is.na(rps[['Rarefied']])), 'Rarefied'] <- 0
      
      df <- tibble::as_tibble(rps, rownames = '.sample')
      
      
      openxlsx::addWorksheet(wb, 'Reads Per Sample')
      openxlsx::writeData(wb, 'Reads Per Sample', df)
  })
  
  
  
  #________________________________________________________
  # "Mapped OTU Counts" and "Rarefied OTU Counts".
  #________________________________________________________
  if (!is.null(rare))
    local ({
      
      tbl <- as.matrix(full$counts) %>% tibble::as_tibble(rownames = ".otu")
      
      if (ncol(full$taxonomy) > 1) {
        lin <- taxa_map(full, rank = 0, unc = unc, lineage = TRUE)
        tbl %<>% mutate(.lineage = unname(lin[tbl$.otu]), .after = ".otu")
      }
      openxlsx::addWorksheet(wb, 'Mapped OTU Counts')
      openxlsx::writeData(wb, 'Mapped OTU Counts', tbl)
      
      
      
      tbl <- as.matrix(rare$counts) %>% tibble::as_tibble(rownames = ".otu")
      
      if (ncol(rare$taxonomy) > 1) {
        lin <- taxa_map(rare, rank = 0, unc = unc, lineage = TRUE)
        tbl %<>% mutate(.lineage = unname(lin[tbl$.otu]), .after = ".otu")
      }
      openxlsx::addWorksheet(wb, 'Rarefied OTU Counts')
      openxlsx::writeData(wb, 'Rarefied OTU Counts', tbl)
      
    })
  
  
  #________________________________________________________
  # "Sample Metadata": .sample, <metadata fields>
  #________________________________________________________
  openxlsx::addWorksheet(wb, 'Sample Metadata')
  openxlsx::writeData(wb, 'Sample Metadata', full$metadata)
  
  
  
  
  #________________________________________________________
  # "Alpha Diversity": .otu, Depth, OTUs, Shannon, ...
  #________________________________________________________
  openxlsx::addWorksheet(wb, 'Alpha Diversity')
  if (!is.null(rare)) {
    openxlsx::writeData(wb, 'Alpha Diversity', {
      adiv_matrix(rare) %>% tibble::as_tibble(rownames = ".otu")
    })
  }
  
  
  
  #________________________________________________________
  # "OTU Taxonomy": .otu, <ranks>, .sequence
  #________________________________________________________
  openxlsx::addWorksheet(wb, 'OTU Taxonomy')
  openxlsx::writeData(wb, 'OTU Taxonomy', local({
    
    tbl <- biom$taxonomy %>%
      relocate('.otu')
    
    if (!is.null(seqs <- full$sequences))
      tbl %<>% mutate(.sequence = unname(seqs[tbl$.otu]))
    
    return (tbl)
  }))
  
  
  
  #________________________________________________________
  # Create worksheets for any other attributes on biom.
  #________________________________________________________
  
  for (key in names(attributes(full))) {
    if (key %in% c("class", "Reads Per Step")) next
    
    val <- attr(full, key, exact = TRUE)
    
    if (identical(class(val),  "dist"))   val <- data.frame(as.matrix(val))
    if (identical(class(val),  "matrix")) val <- data.frame(val)
    if (!identical(class(val), "data.frame")) next
    
    rn <- !identical(rownames(val), as.character(1:nrow(val)))
    openxlsx::addWorksheet(wb, key)
    openxlsx::writeData(wb, key, val, rowNames = rn)
    remove("val", "rn")
  }
  
  remove(list = intersect(ls(), "key"))
  
  
  
  #________________________________________________________
  # Roll up abundances at each taxonomic level.
  #________________________________________________________
  if (is.null(rare))
    rare <- full
  
  for (rank in rare$ranks[-1]) {
    openxlsx::addWorksheet(wb, rank)
    taxa_matrix(rare, rank, lineage = TRUE) %>%
      tibble::as_tibble(rownames = paste0('.', tolower(rank))) %>%
      openxlsx::writeData(wb = wb, sheet = rank)
  }
  
  remove(list = intersect(ls(), c("rank")))
  
  
  
  #________________________________________________________
  # Write everything to the specified output file
  #________________________________________________________
  
  openxlsx::saveWorkbook(wb, file=file, overwrite=TRUE)
  
  
    
  
  return (invisible(file))
}
