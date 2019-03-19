#' Write data and summary information to a Microsoft Excel-compatible workbook.
#'
#' @param biom      The BIOM object to save to the file.
#' @param outfile   Path to the output xlsx file.
#' @param depth     Depth to rarefy to. See \code{rarefy} function for details.
#' @param seed      Random seed to use in rarefying. See \code{rarefy} function
#'                    for details.
#' @return On success, returns \code{NULL} invisibly.
#' @export


write.xlsx <- function (biom, outfile, depth=NULL, seed=0) {
  
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
  openxlsx::addWorksheet(wb, "Raw OTU Counts")
  openxlsx::addWorksheet(wb, "Rarefied OTU Counts")
  openxlsx::addWorksheet(wb, 'Alpha Diversity')
  openxlsx::addWorksheet(wb, 'Taxonomy')
  
  
  
  #========================================================
  # Output counts and taxonomy, & seqs prior to rarefying
  #========================================================
  
  RawCounts <- data.frame(rbiom::counts(biom),   check.names=FALSE)
  Taxonomy  <- data.frame(rbiom::taxonomy(biom), check.names=FALSE)
  
  if (is(biom[['sequences']], "character")) {
    Taxonomy[['Sequence']] <- unname(biom[['sequences']][rownames(Taxonomy)])
    if (ncol(Taxonomy) > 1)
      Taxonomy <- Taxonomy[,c(ncol(Taxonomy), 1:(ncol(Taxonomy) - 1)),drop=FALSE]
  }
  
  openxlsx::writeData(wb, "Raw OTU Counts", RawCounts, rowNames=TRUE)
  openxlsx::writeData(wb, 'Taxonomy',       Taxonomy,  rowNames=TRUE)
  
  
  
  #========================================================
  # Rarefy, then output counts again and alpha.div
  #========================================================
  
  rare <- rbiom::rarefy(biom, depth, seed)
  
  openxlsx::writeData(wb, "Rarefied OTU Counts", data.frame(rbiom::counts(rare),  check.names=FALSE), rowNames=TRUE)
  openxlsx::writeData(wb, 'Alpha Diversity',     rbiom::alpha.div(rare),                              rowNames=TRUE)
  
  
  
  #========================================================
  # Track the each sample's read counts before and after
  #========================================================
  
  df <- merge(
    data.frame(slam::col_sums(biom$counts)),
    data.frame(slam::col_sums(rare$counts)),
    by="row.names", all=TRUE)
  names(df) <- c("Sample", "Raw", "Rarefied")
  df <- df[order(df$Sample),]
  df[which(is.na(df[['Rarefied']])), 'Rarefied'] <- 0
  openxlsx::writeData(wb, 'Reads Per Sample', df, rowNames=FALSE)
  
  biom <- rare
  
  remove("df", "rare")
  
  
  
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
  # Skip Distance Matrices and Ordinations unless >1 sample
  #========================================================
  if (rbiom::nsamples(biom) >= 2) {
    
    
    #========================================================
    # Distance Matrices
    #========================================================
    
    if (is(biom[['phylogeny']], 'character')) {
      distance_matrices <- list( U.UniFrac    = rbiom::beta.div(biom, 'unifrac',     FALSE),
                                 W.UniFrac    = rbiom::beta.div(biom, 'unifrac',     TRUE) )
    } else {
      distance_matrices <- list( W.BrayCurtis = rbiom::beta.div(biom, 'bray-curtis', TRUE),
                                 W.Jaccard    = rbiom::beta.div(biom, 'jaccard',     TRUE) )
    }
  
    for (sheetName in names(distance_matrices)) {
      df <- data.frame(as.matrix(distance_matrices[[sheetName]]), check.names=FALSE)
      df <- df[sort(rownames(df)), sort(colnames(df)), drop=FALSE]
      sheetName <- paste(sheetName, " DM")
      openxlsx::addWorksheet(wb, sheetName)
      openxlsx::writeData(wb, sheetName, df, rowNames=TRUE)
    }
    
    remove("df", "sheetName", "distance_matrices")
    
  }
  
  
  openxlsx::saveWorkbook(wb, file=outfile, overwrite=TRUE)
  
  
  return(invisible(NULL))
}
