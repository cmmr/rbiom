#' Map OTUs names to taxa names at a given rank.
#' 
#' @inherit documentation_default
#' 
#' @family taxonomy
#' 
#' @param rank  When `NULL`, the entire biom$taxonomy data.frame is returned, 
#'        transformed as per `unc`. Alternatively, a single taxonomic rank 
#'        (rank name or integer position in `biom$ranks`) which returns a named
#'        character vector for mapping OTUs to taxa names.
#' 
#' @return A tibble data.frame when `rank=NULL`, or a character vector named 
#'         with the OTU names.
#' 
#' @seealso `pull.rbiom()`
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     library(dplyr, warn.conflicts = FALSE)
#'     
#'     # In $taxonomy, .otu is the first column (like a row identifier)  -----
#'     hmp50$taxonomy %>% head(4)
#'     
#'     # In taxa_map, .otu is the last column (most precise rank)  -----------
#'     taxa_map(hmp50) %>% head(4)
#'     
#'     # Generate an OTU to Genus mapping  -----------------------------------
#'     taxa_map(hmp50, "Genus") %>% head(4)
#'     
#'     # Sometimes taxonomic names are incomplete ----------------------------
#'     otus <- c('GemAsacc', 'GcbBacte', 'Unc58411')
#'     taxa_map(hmp50) %>% filter(.otu %in% otus) %>% select(Phylum:.otu)
#'     
#'     # rbiom can replace them with unique placeholders ---------------------
#'     taxa_map(hmp50, unc = "singly") %>% filter(.otu %in% otus) %>% select(Class:.otu)
#'     
#'     # Or collapse them into groups ----------------------------------------
#'     taxa_map(hmp50, unc = "grouped") %>% filter(.otu %in% otus) %>% select(Class:Genus)
#'

taxa_map <- function (biom, rank = NULL, unc = "singly", lineage = FALSE) {
  
  biom <- as_rbiom(biom)
  
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- eval_envir(environment())
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  
  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  validate_rank(null_ok = TRUE)
  validate_unc()
  validate_bool("lineage")
  
  
  #________________________________________________________
  # Move '.otu' to last column of taxonomy map.
  #________________________________________________________
  tbl <- relocate(biom$taxonomy, .otu, .after = last_col())
  if (!is.null(rank)) rank <- which(colnames(tbl) == rank)
  
  
  
  #________________________________________________________
  # Transform the taxa names.
  #________________________________________________________
  if (unc != "asis" && ncol(tbl) > 1)
    tbl <- tryCatch(
      expr = local({
        
        mtx  <- as.matrix(tbl)
        otus <- mtx[,'.otu']
        
      
        # Discard technical prefixes/suffixes.
        #________________________________________________________
        mtx <- sub("^.__", "", mtx) # Remove leading p__ c__ etc
        mtx <- sub("^_+",  "", mtx) # Remove leading underscores
        mtx <- sub(";$",   "", mtx) # Remove trailing semicolons
        
        
        # "g" => NA; "Unknown Order" => NA
        #________________________________________________________
        regex <- ".*(unknown|uncultured|unclassified|unidentified|incertae.sedis).*"
        mtx[which(nchar(mtx) < 2)] <- NA
        mtx[grep(regex, mtx, ignore.case = TRUE)] <- NA
        
        
        # "R_7_group" => "R_7"
        #________________________________________________________
        mtx[] <- sub("_group$", "", mtx, ignore.case = TRUE)
        
        
        # "Family XIII" => "Clostridiales XIII"
        #________________________________________________________
        mtx <- t(apply(mtx, 1L, function (x) {
          regex <- "^family\\s"
          if (!is.null(prefixed <- grep(regex, x, ignore.case = TRUE)))
            for (i in prefixed)
              if (!is.null(ideal <- grep("^[A-Z][a-z]+$", x[seq_len(i - 1)], value = TRUE)))
                x[i] <- sub(regex, paste0(tail(ideal, 1), " "), x[i], ignore.case = TRUE)
          
          return (x)
        }))
        
        
        # Replace NA with "Unc. <OTU ID>".
        #________________________________________________________
        if (eq(unc, "singly")) {
          x <- which(is.na(mtx))
          mtx[x] <- paste("Unc.", mtx[row(mtx)[x], '.otu'])
        }
        
        
        # Replace NA with "Unc. <Higher Rank>".
        #________________________________________________________
        if (eq(unc, "grouped"))
          for (i in which(!complete.cases(mtx)))
            for (j in rev(which(is.na(mtx[i,]))))
              if (!is.null(x <- na.omit(c("N/A", mtx[i,seq_len(j - 1)]))))
                mtx[i,j] <- paste("Unc.", tail(x, 1))
        
        
        # Drop any row with an NA in a higher-order rank column.
        #________________________________________________________
        if (eq(unc, "drop"))
          mtx <- mtx[complete.cases(mtx[,1:if.null(rank, ncol(mtx)),drop=FALSE]),,drop=FALSE]
        
        
        
        tbl <- as_tibble(mtx) %>%
          relocate(.otu) %>%
          mutate(across(everything(), as.factor))
        
        
        # Always keep .otu column unmodified.
        #________________________________________________________
        tbl[['.otu']] <- as.character(otus)
        
        
        return (tbl)
      }), 
      
      error = function (e)
        stop("Error in renaming taxa: ", e) )
  
  
  
  tbl %<>% relocate(.otu, .after = last_col())
  if (!is.null(rank)) {
    
    if (isTRUE(lineage)) {
      tbl <- setNames(
        object = plyr::splat(paste)(as.list(tbl[,1:rank]), sep = "; "),
        nm     = as.character(tbl[['.otu']]) )
      
    } else {
      tbl <- setNames(tbl[[rank]], as.character(tbl[['.otu']]))
    }
  }
  
  
  set_cache_value(cache_file, tbl)
  return (tbl)
}
