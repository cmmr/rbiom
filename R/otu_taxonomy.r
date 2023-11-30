#' Get or set the taxonomy table.
#' 
#' @inherit documentation_default
#' 
#' @family taxonomy
#' 
#' 
#' @param rank  A single taxonomic rank to return as a named vector, or 
#'        \code{NULL} to return the complete table. Default: \code{NULL}
#' 
#' @param value  An object coercible with [tibble::as_tibble()]. Must either 
#'        have rownames or an '.otu' column with OTU names matching 
#'        [otu_names(x)].
#' 
#' @return Depending on \code{rank}, a named character vector or a tibble data 
#'         frame with '.otu' as the first column name.
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     # Display the full taxonomic data for each OTU ------------------------
#'     otu_taxonomy(hmp50)
#'     
#'     # Only show the OTU -> Genus mapping  ---------------------------------
#'     otu_taxonomy(hmp50, "Genus") %>% head()
#'     
#'     # Sometimes taxonomic names are incomplete ----------------------------
#'     otu_taxonomy(hmp50)[,4:7] %>%
#'       dplyr::filter(.otu %in% c('GemAsacc', 'GcbBacte', 'Unc58411'))
#'     
#'     # rbiom can insert more descriptive placeholders ----------------------
#'     otu_taxonomy(hmp50, unc = "singly")[,4:7] %>%
#'       dplyr::filter(.otu %in% c('GemAsacc', 'GcbBacte', 'Unc58411'))
#'     
#'     # Or collapse them into groups ----------------------------------------
#'     otu_taxonomy(hmp50, unc = "grouped")[,4:7] %>%
#'       dplyr::filter(.otu %in% c('GemAsacc', 'GcbBacte', 'Unc58411'))
#'

otu_taxonomy <- function (biom, rank = NULL, unc = "asis", lineage = FALSE) {
  
  validate_biom(clone = FALSE)
  
  
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
  tbl <- relocate(biom[['taxonomy']], .otu, .after = last_col())
  if (!is.null(rank)) rank <- which(colnames(tbl) == rank)
  
  
  
  #________________________________________________________
  # Transform the taxa names.
  #________________________________________________________
  if (unc != "asis")
    tbl <- tryCatch(
      expr = local({
        
        mtx <- tbl %>% 
          as.matrix()
        
      
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


#' @rdname otu_taxonomy
#' @export

`otu_taxonomy<-` <- function (biom, value) {
  
  validate_biom(clone = TRUE)
  
  
  #________________________________________________________
  # Parse a file/URL.
  #________________________________________________________
  if (is_scalar_character(value)) {
    value <- import_table(value)
    colnames(value) <- c('.otu', default_taxa_ranks(ncol(value) - 1))
  }
  
  otus <- otu_names(biom)
  
  
  
  #________________________________________________________
  # Convert data.frame/matrix to tibble with .otu column.
  #________________________________________________________
  value <- as_tibble(value, rownames = NA)
  stopifnot(xor(hasName(value, '.otu'), has_rownames(value)))
  
  if (has_rownames(value))
    value %<>% rownames_to_column(var = ".otu")
  
  value %<>% relocate(.otu)
  value %<>% mutate(across(everything() & !.otu, as.factor))
  
  
  
  #________________________________________________________
  # Ignore OTUs that aren't currently in the rbiom object.
  #________________________________________________________
  ignored <- setdiff(as.character(value[['.otu']]), otus)
  n       <- length(ignored)
  if (n > 0) {
    if (n > 4) ignored <- c(head(ignored, 4), "...")
    ignored <- paste(collapse = ", ", ignored)
    msg <- "Ignoring %i extra OTU ID%s: %s"
    message(sprintf(msg, n, ifelse(n == 1, "", "s"), ignored))
  }
  
  
  
  #________________________________________________________
  # Drop any OTUs that aren't in the new taxonomy map.
  #________________________________________________________
  drop <- setdiff(otus, as.character(value[['.otu']]))
  n    <- length(drop)
  if (n > 0) {
    if (n > 4) drop <- c(head(drop, 4), "...")
    drop <- paste(collapse = ", ", drop)
    msg <- "Dropping %i OTU%s from the otu matrix: %s"
    message(sprintf(msg, n, ifelse(n == 1, "", "s"), drop))
  }
  
  
  
  #________________________________________________________
  # Reorder and subset to match incoming OTU ids.
  #________________________________________________________
  otus <- intersect(as.character(value[['.otu']]), otus)
  biom[['counts']]   <- biom[['counts']][otus,]
  biom[['taxonomy']] <- left_join(
    x  = tibble(.otu = otus),
    y  = value, 
    by = ".otu" )
  
  
  
  biom <- biom_repair(biom)
  
  invalidate_biom()
  return (biom) 
}


