#' Write data and summary information to a Microsoft Excel-compatible workbook.
#'
#' @param biom     The BIOM object to save to the file.
#' @param outfile  Path to the output xlsx file.
#' @param depth    Depth to rarefy to. See \code{rarefy} function for details.
#'                 \code{depth = NULL} auto-selects a rarefaction level.
#'                 \code{depth = 0} disables rarefaction.
#'                   Only use \code{depth} with \code{BIOM} files of type
#'                   'OTU table' and integer count values.
#' @param seed     Random seed to use in rarefying. See \code{rarefy} function
#'                   for details.
#' @return On success, returns \code{NULL} invisibly.
#' 
#' @section Note:
#' Any \code{data frame}, \code{matrix}, or \code{dist} attributes on 
#' \code{biom} will be included as separate worksheets. An attribute named 
#' 'Reads Per Step' is treated specially and merged with the usual 'Reads Per 
#' Sample' tab.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- select(hmp50, 1:10) %>% rarefy()
#'     
#'     attr(biom, "Weighted UniFrac")   <- unifrac(biom)
#'     attr(biom, "Unweighted Jaccard") <- bdiv_distmat(biom, 'jaccard', weighted=FALSE)
#'     
#'     outfile <- tempfile(fileext = ".xlsx")
#'     write_xlsx(biom, outfile)
#'


write_xlsx <- function (biom, outfile, depth=NULL, seed=0) {
  
  #========================================================
  # Sanity Checks
  #========================================================
  
  if (!is(biom, "BIOM"))
    stop(simpleError("Invalid BIOM object."))
  
  
  
  #========================================================
  # Biom is an 'OTU table' and all counts are int
  #========================================================
  
  if (identical(tolower(biom$info$type), 'otu table') && !any(biom$counts$v %% 1 > 0)) {
    
    #--------------------------------------------------------
    # Define the initial structure of our workbook
    #--------------------------------------------------------
    
    wb <- openxlsx::createWorkbook(creator="CMMR", title=biom$info$id)
    
    openxlsx::addWorksheet(wb, 'Reads Per Sample')
    openxlsx::addWorksheet(wb, 'Mapped OTU Counts')
    openxlsx::addWorksheet(wb, 'Rarefied OTU Counts')
    openxlsx::addWorksheet(wb, 'Alpha Diversity')
    openxlsx::addWorksheet(wb, 'Centroid Taxonomy')
    
    
    
    #--------------------------------------------------------
    # Output counts and taxonomy, & seqs prior to rarefying
    #--------------------------------------------------------
    
    RawCounts <- data.frame(counts(biom),   check.names=FALSE)
    Taxonomy  <- data.frame(taxonomy(biom), check.names=FALSE)
    
    # Set the full taxonomy string as the first column of RawCounts
    TaxaStrings <- data.frame(Taxonomy=apply(biom$taxonomy, 1L, paste, collapse="; "))
    RawCounts   <- cbind(TaxaStrings, RawCounts[rownames(TaxaStrings),,F])
    
    # Set the centroid sequence as the first column of Taxonomy
    if (is(biom[['sequences']], "character"))
      Taxonomy <- cbind(Sequence=biom[['sequences']], Taxonomy)
    
    openxlsx::writeData(wb, 'Mapped OTU Counts', RawCounts, rowNames=TRUE)
    openxlsx::writeData(wb, 'Centroid Taxonomy', Taxonomy,  rowNames=TRUE)
    
    
    
    #--------------------------------------------------------
    # Rarefy, then output counts again and adiv_table
    #--------------------------------------------------------
    
    # Allow the user to override rarefaction by setting depth = 0.
    if (identical(depth, 0) || isTRUE(is_rarefied(biom))) {
      rare <- biom
    } else {
      rare <- rarefy(biom, depth, seed)
    }
    
    RareCounts <- data.frame(counts(rare), check.names=FALSE)
    AlphaDiv   <- adiv_table(rare)
    
    # Set the full taxonomy string as the first column of RareCounts
    TaxaStrings <- data.frame(Taxonomy=apply(rare$taxonomy, 1L, paste, collapse="; "))
    RareCounts  <- cbind(TaxaStrings, RareCounts[rownames(TaxaStrings),,F])
    
    openxlsx::writeData(wb, 'Rarefied OTU Counts', RareCounts, rowNames=TRUE)
    openxlsx::writeData(wb, 'Alpha Diversity',     AlphaDiv,   rowNames=FALSE)
    
    
    
    #--------------------------------------------------------
    # Track each sample's read counts before and after
    #--------------------------------------------------------
    
    if ('Reads Per Step' %in% names(attributes(biom))) {
      rps <- attr(biom, 'Reads Per Step', exact = TRUE)
    } else {
      rps <- data.frame(Mapped=slam::col_sums(biom$counts))
    }
    
    rps <- data.frame(rps[order(rownames(rps)),,drop=FALSE])
    rps[['Rarefied']] <- slam::col_sums(rare$counts)[rownames(rps)]
    rps[which(is.na(rps[['Rarefied']])), 'Rarefied'] <- 0
    openxlsx::writeData(wb, 'Reads Per Sample', rps, rowNames=TRUE)
    
    biom <- rare
    
    remove(list = intersect(ls(), "rare"))
    
    
    
    #--------------------------------------------------------
    # Create worksheets for any other attributes on biom
    #--------------------------------------------------------
    
    for (key in names(attributes(biom))) {
      if (key %in% c("names", "class", "Reads Per Step"))   next
      
      val <- attr(biom, key, exact = TRUE)
      
      if (identical(class(val), "dist"))   val <- data.frame(as.matrix(val))
      if (identical(class(val), "matrix")) val <- data.frame(val)
      if (!identical(class(val), "data.frame")) next
      
      rn <- !identical(rownames(val), as.character(1:nrow(val)))
      openxlsx::addWorksheet(wb, key)
      openxlsx::writeData(wb, key, val, rowNames = rn)
      remove("val", "rn")
    }
    
    remove(list = intersect(ls(), "key"))
    
    
    
    #--------------------------------------------------------
    # Roll up abundances at each taxonomic level
    #--------------------------------------------------------
    
    ranks <- unique(c('OTU', taxa_ranks(biom)))
    
    for (rank in ranks) {
      
      df <- t(taxa_rollup(biom, rank, lineage = TRUE, safe=NA))
      df <- df[sort(rownames(df)), sort(colnames(df)), drop=FALSE]
      
      openxlsx::addWorksheet(wb, rank)
      openxlsx::writeData(wb, rank, df, rowNames=TRUE)
    }
    
    remove(list = intersect(ls(), c("df", "rank", "ranks")))
    
    
    
    #--------------------------------------------------------
    # Write everything to the specified output file
    #--------------------------------------------------------
    
    openxlsx::saveWorkbook(wb, file=outfile, overwrite=TRUE)
    
    
    
    
  #========================================================
  # Biom isn't an 'OTU table' or not all counts are int
  #========================================================
  
  } else {
    
    type <- sub(" table", "", biom$info$type, ignore.case=TRUE)
    
    if (!is.null(depth))
      return(simpleError("Please use write_xlsx(depth=NULL) for non-OTU or non-Integer data."))
    
    
    #--------------------------------------------------------
    # Define the initial structure of our workbook
    #--------------------------------------------------------
    
    wb <- openxlsx::createWorkbook(creator="CMMR", title=biom$info$id)
    
    ranks <- c(paste0(type, "s"), taxa_ranks(biom))
    for (rank in ranks)
      openxlsx::addWorksheet(wb, rank)
    
    openxlsx::addWorksheet(wb, 'Alpha Diversity')
    openxlsx::addWorksheet(wb, paste(type, 'List'))
    
    
    
    #--------------------------------------------------------
    # Output totals, taxonomy/seqs, & adiv_table
    #--------------------------------------------------------
    
    AlphaDiv  <- data.frame(adiv_table(biom), check.names=FALSE)
    Taxonomy  <- data.frame(taxonomy(biom),  check.names=FALSE)
    
    # Set the centroid sequence as the first column of Taxonomy
    if (is(biom[['sequences']], "character"))
      Taxonomy <- cbind(Sequence=biom[['sequences']], Taxonomy)
    
    # Rename a few of the Alpha Div columns
    names(AlphaDiv) <- sub("Depth", "Total",   names(AlphaDiv))
    names(AlphaDiv) <- sub("OTUs", ranks[[1]], names(AlphaDiv))
    
    openxlsx::writeData(wb, 'Alpha Diversity',   AlphaDiv, rowNames=FALSE)
    openxlsx::writeData(wb, paste(type, 'List'), Taxonomy, rowNames=TRUE)
    
    
    
    #--------------------------------------------------------
    # Create worksheets for any other attributes on biom
    #--------------------------------------------------------
    
    for (key in names(attributes(biom))) {
      if (key %in% c("names", "class")) next
      
      val <- attr(biom, key, exact = TRUE)
      
      if (identical(class(val), "dist"))   val <- data.frame(as.matrix(val))
      if (identical(class(val), "matrix")) val <- data.frame(val)
      if (!identical(class(val), "data.frame")) next
      
      rn <- !identical(rownames(val), as.character(1:nrow(val)))
      openxlsx::addWorksheet(wb, key)
      openxlsx::writeData(wb, key, val, rowNames = rn)
      remove("val", "rn")
    }
    
    remove(list = intersect(ls(), "key"))
    
    
    
    #--------------------------------------------------------
    # Roll up abundances at each taxonomic level
    #--------------------------------------------------------
      
    for (i in seq_along(ranks)) {
      
      df <- t(taxa_rollup(biom, i - 1, lineage=TRUE, safe=NA))
      df <- df[sort(rownames(df)), sort(colnames(df)), drop=FALSE]
      
      # Set the full taxonomy string as the first column
      if (identical(i, 1L)) {
        Description <- data.frame(Description=apply(biom$taxonomy, 1L, paste, collapse="; "))
        df          <- cbind(Description, df[rownames(Description),,F])
      }
      
      openxlsx::writeData(wb, ranks[[i]], df, rowNames=TRUE)
    }
    
    remove(list = intersect(ls(), c("df", "i")))
    
    
    
    #--------------------------------------------------------
    # Write everything to the specified output file
    #--------------------------------------------------------
    
    openxlsx::saveWorkbook(wb, file=outfile, overwrite=TRUE)
    
  }
  
  return(invisible(NULL))
}