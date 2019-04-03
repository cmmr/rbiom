#' Write data and summary information to a Microsoft Excel-compatible workbook.
#'
#' @param biom     The BIOM object to save to the file.
#' @param outfile  Path to the output xlsx file.
#' @param depth    Depth to rarefy to. See \code{rarefy} function for details.
#' @param seed     Random seed to use in rarefying. See \code{rarefy} function
#'                   for details.
#' @param rps    A \code{data frame} of additional information to include in
#'                 the 'Reads per Sample' tab. The rownames of this 
#'                 \code{data frame} must match the IDs in the BIOM object.
#' @return On success, returns \code{NULL} invisibly.
#' @export


write.xlsx <- function (biom, outfile, depth=NULL, seed=0, rps=NULL) {
  
  #========================================================
  # Sanity Checks
  #========================================================
  
  if (!is(biom, "BIOM"))
    stop(simpleError("Invalid BIOM object."))
  
  
  
  #========================================================
  # Define the initial structure of our workbook
  #========================================================
  
  wb <- openxlsx::createWorkbook(creator="CMMR", title=biom$info$id)
  
  openxlsx::addWorksheet(wb, 'Reads Per Sample')
  openxlsx::addWorksheet(wb, "Mapped OTU Counts")
  openxlsx::addWorksheet(wb, "Rarefied OTU Counts")
  openxlsx::addWorksheet(wb, 'Alpha Diversity')
  openxlsx::addWorksheet(wb, 'Centroid Taxonomy')
  
  
  
  #========================================================
  # Output counts and taxonomy, & seqs prior to rarefying
  #========================================================
  
  RawCounts <- data.frame(rbiom::counts(biom),   check.names=FALSE)
  Taxonomy  <- data.frame(rbiom::taxonomy(biom), check.names=FALSE)
  
  # Set the full taxonomy string as the first column of RawCounts
  TaxaStrings <- data.frame(Taxonomy=apply(biom$taxonomy, 1L, paste, collapse="; "))
  RawCounts   <- cbind(TaxaStrings, RawCounts[rownames(TaxaStrings),,F])
  
  # Set the centroid sequence as the first column of Taxonomy
  if (is(biom[['sequences']], "character"))
    Taxonomy <- cbind(Sequence=biom[['sequences']], Taxonomy)
  
  openxlsx::writeData(wb, "Mapped OTU Counts", RawCounts, rowNames=TRUE)
  openxlsx::writeData(wb, 'Centroid Taxonomy', Taxonomy,  rowNames=TRUE)
  
  
  
  #========================================================
  # Rarefy, then output counts again and alpha.div
  #========================================================
  
  rare <- rbiom::rarefy(biom, depth, seed)
  
  RareCounts <- data.frame(rbiom::counts(rare), check.names=FALSE)
  AlphaDiv   <- rbiom::alpha.div(rare)
  
  # Set the full taxonomy string as the first column of RareCounts
  TaxaStrings <- data.frame(Taxonomy=apply(rare$taxonomy, 1L, paste, collapse="; "))
  RareCounts  <- cbind(TaxaStrings, RareCounts[rownames(TaxaStrings),,F])
  
  openxlsx::writeData(wb, "Rarefied OTU Counts", RareCounts, rowNames=TRUE)
  openxlsx::writeData(wb, 'Alpha Diversity',     AlphaDiv,   rowNames=FALSE)
  
  
  
  #========================================================
  # Track the each sample's read counts before and after
  #========================================================
  
  if (is.null(rps))
    rps <- data.frame(Mapped=slam::col_sums(biom$counts))
  
  rps <- data.frame(rps[order(rownames(rps)),,drop=FALSE])
  rps[['Rarefied']] <- slam::col_sums(rare$counts)[rownames(rps)]
  rps[which(is.na(rps[['Rarefied']])), 'Rarefied'] <- 0
  openxlsx::writeData(wb, 'Reads Per Sample', rps, rowNames=TRUE)
  
  biom <- rare
  
  remove("rare")
  
  
  
  #========================================================
  # Roll up abundances at each taxonomic level
  #========================================================
  
  ranks <- c('OTU', rbiom::taxa.ranks(biom))
  
  for (rank in ranks) {
    
    df <- data.frame(rbiom::taxa.rollup(biom, rank, lineage = TRUE), check.names=FALSE)
    df <- df[sort(rownames(df)), sort(colnames(df)), drop=FALSE]
    
    openxlsx::addWorksheet(wb, rank)
    openxlsx::writeData(wb, rank, df, rowNames=TRUE)
  }
  
  remove("df", "rank", "ranks")
  
  
  
  #========================================================
  # Write everything to the specified output file
  #========================================================
  
  openxlsx::saveWorkbook(wb, file=outfile, overwrite=TRUE)
  
  
  return(invisible(NULL))
}
