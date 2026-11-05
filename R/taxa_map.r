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
#' @param taxa   Which taxa to keep in the mapping. An integer value keeps the
#'        top n most abundant taxa, and a value 0 <= n < 1 keeps taxa with
#'        that mean relative abundance or greater, both ranked exactly as in
#'        [taxa_matrix()]. A character vector of taxa names keeps only those
#'        taxa. OTUs belonging to other taxa are dropped from the mapping, or
#'        mapped to `other` when it is given. Default: `Inf` (all taxa)
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
#'     # Keep the five most abundant OTUs; map the rest to "Other" -----------
#'     taxa_map(hmp50, ".otu", taxa = 5, other = TRUE) %>% table()
#'
#'     # Sometimes taxonomic names are incomplete ----------------------------
#'     otus <- c('GemAsacc', 'GcbBacte', 'Unc58411')
#'     taxa_map(hmp50, unc = "asis") %>% filter(.otu %in% otus) %>% select(Phylum:.otu)
#'
#'     # rbiom can replace them with unique placeholders ---------------------
#'     taxa_map(hmp50, unc = "singly") %>% filter(.otu %in% otus) %>% select(Class:.otu)
#'
#'     # Or collapse them into groups ----------------------------------------
#'     taxa_map(hmp50, unc = "grouped") %>% filter(.otu %in% otus) %>% select(Class:Genus)
#'

taxa_map <- function (biom, rank = NULL, taxa = Inf, unc = "singly", lineage = FALSE, other = FALSE) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  cmd    <- current_cmd('taxa_map', params)


  #________________________________________________________
  # Validate user's arguments.
  #________________________________________________________
  validate_rank(null_ok = TRUE)
  validate_taxa()
  validate_unc()
  validate_bool("lineage")
  validate_string('other', if_true = 'Other', if_false = NULL, null_ok = TRUE)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('taxa_map', params)
  tbl        <- get_cache_value(cache_file)

  if (is.null(tbl)) {

    #________________________________________________________
    # Move '.otu' to last column of taxonomy map. `rank` stays a rank
    # name; `rank_idx` is its column position in this table.
    #________________________________________________________
    tbl      <- relocate(biom$taxonomy, '.otu', .after = dplyr::last_col())
    rank_idx <- if (is.null(rank)) ncol(tbl) else which(colnames(tbl) == rank)


    #________________________________________________________
    # Transform the taxa names.
    #________________________________________________________
    if (unc != "asis" && ncol(tbl) > 1)
      tbl <- tryCatch(
        error = function (e) cli_abort(c("Error in renaming taxa.", 'x' = conditionMessage(e))),
        expr  = local({

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
            for (i in grep(regex, x, ignore.case = TRUE)) {
              ideal <- grep("^[A-Z][a-z]+$", x[seq_len(i - 1)], value = TRUE)
              if (length(ideal) > 0)
                x[i] <- sub(regex, paste0(tail(ideal, 1), " "), x[i], ignore.case = TRUE)
            }
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
            for (i in which(!stats::complete.cases(mtx)))
              for (j in rev(which(is.na(mtx[i,])))) {
                higher   <- na.omit(c("N/A", mtx[i, seq_len(j - 1)]))
                mtx[i,j] <- paste("Unc.", tail(higher, 1))
              }


          # Drop any row with an NA in a higher-order rank column.
          #________________________________________________________
          if (eq(unc, "drop")) {
            keep <- stats::complete.cases(mtx[,seq_len(rank_idx),drop=FALSE])
            otus <- otus[keep]
            mtx  <- mtx[keep,,drop=FALSE]
          }


          tbl <- as_tibble(mtx) %>%
            relocate('.otu') %>%
            mutate(across(everything(), as.factor))


          # Always keep .otu column unmodified.
          #________________________________________________________
          tbl[['.otu']] <- as.character(otus)


          return (tbl)
        }))


    tbl %<>% relocate('.otu', .after = dplyr::last_col())

    if (is.null(rank)) {
      tbl %<>% as_rbiom_tbl()

    } else {

      if (isTRUE(lineage)) {
        tbl <- setNames(
          object = do.call(paste, c(as.list(tbl[, seq_len(rank_idx)]), sep = "; ")),
          nm     = as.character(tbl[['.otu']]) )

      } else {
        tbl <- setNames(tbl[[rank]], as.character(tbl[['.otu']]))
      }


      #________________________________________________________
      # Drop or rename taxa to 'Other'. Top taxa are chosen exactly as
      # taxa_matrix() does: by mean relative abundance.
      #________________________________________________________
      if (!all(is.infinite(taxa))) {

        if (is.numeric(taxa)) {
          mtx  <- taxa_matrix(biom, rank = rank, lineage = lineage, unc = unc, sparse = TRUE)
          taxa <- top_taxa(mtx, taxa)
        }

        tbl <- factor(tbl, levels = taxa)

        if (is.null(other)) {
          tbl <- tbl[!is.na(tbl)]
        } else {
          tbl <- addNA(tbl)
          levels(tbl) <- ifelse(is.na(levels(tbl)), other, levels(tbl))
        }
      }
    }

    set_cache_value(cache_file, tbl)
  }


  return (prov_result(tbl, cmd, biom))
}
