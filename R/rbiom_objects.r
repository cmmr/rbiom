#' Working with rbiom Objects.
#' 
#' \cr
#' Rbiom objects make it easy to access and manipulate your BIOM data, ensuring 
#' all the disparate components remain in sync. These objects behave largely 
#' like lists, in that you can access and assign to them using the `$` 
#' operator. The sections below list all the fields which can be read and/or 
#' written, and the helper functions for common tasks like rarefying and 
#' subsetting. To create an rbiom object, see [as_rbiom()].
#' \cr\cr
#' Use `$clone()` to create a copy of an rbiom object. This is necessary 
#' because rbiom objects are **passed by reference**. The usual `<-` assignment 
#' operator will simply create a second reference to the same object - it will 
#' not create a second object. See [speed ups][speed] for more details.
#' \cr\cr
#' 
#' 
#' @name rbiom_objects
#' @keywords internal
#' 
#' 
#' 
#' @section Readable Fields:
#' 
#' Reading from fields will not change the rbiom object.
#' 
#' | **Accessor**             | **Content**                                              |
#' | ------------------------ | -------------------------------------------------------- |
#' | `$counts`                | Abundance of each OTU in each sample.                    |
#' | `$metadata`              | Sample mappings to metadata (treatment, patient, etc).   |
#' | `$taxonomy`              | OTU mappings to taxonomic ranks (genus, phylum, etc).    |
#' | `$otus`, `$n_otus`       | OTU names.                                               |
#' | `$samples`, `$n_samples` | Sample names.                                            |
#' | `$fields`, `$n_fields`   | Metadata field names.                                    |
#' | `$ranks`, `$n_ranks`     | Taxonomic rank names.                                    |
#' | `$tree`, `$sequences`    | Phylogenetic tree / sequences for the OTUs, or `NULL`.   |
#' | `$id`, `$comment`        | Arbitrary strings for describing the dataset, or `NULL`. |
#' | `$depth`                 | Rarefaction depth, or `NULL` if unrarefied.              |
#' | `$date`                  | Date from BIOM file.                                     |
#' 
#' \cr
#' 
#' 
#' @section Writable Fields:
#' 
#' Assigning new values to these components will trigger 
#' validation checks and inter-component synchronization.
#' 
#' | **Component**     | **What can be assigned.**                               |
#' | ----------------- | ------------------------------------------------------- |
#' | `$counts`         | matrix of abundances; OTUs (rows) by samples (columns)  |
#' | `$metadata`       | data.frame with `'.sample'` column, or a file name      |
#' | `$taxonomy`       | data.frame with `'.otu'` as the first column            |
#' | `$otus`           | character vector with new names for the OTUs            |
#' | `$samples`        | character vector with new names for the samples         |
#' | `$tree`           | phylo object with the phylogenetic tree for the OTUs    |
#' | `$sequences`      | character vector of OTU reference sequences             |
#' | `$id`, `$comment` | string with a title or comment for the dataset          |
#' | `$date`           | date-like object, or `"%Y-%m-%dT%H:%M:%SZ"` string      |
#' 
#' \cr
#' 
#' 
#' @section Transformations:
#' 
#' All functions return an rbiom object.
#' 
#' | **Function**                            | **Transformation**                               |
#' | --------------------------------------- | ------------------------------------------------ |
#' | \code{\emph{<rbiom>}$clone()}           | Safely duplicate an rbiom object.                |
#' | \code{\link[=subset]{\emph{<rbiom>}[]}} | Subset to a specific set of sample names.        |
#' | [subset()]                              | Subset samples according to metadata properties. |
#' | [slice()]                               | Subset to a specific number of samples.          |
#' | [mutate()]                              | Create, modify, and delete metadata fields.      |
#' | [rarefy()]                              | Sub-sample OTU counts to an even sampling depth. |
#' 
#' \cr
#' 
#' 
#' @examples
#'     library(rbiom)
#'     
#'     # Duplicate the HMP50 example dataset.
#'     biom <- hmp50$clone()
#'     
#'     
#'     # Display an overall summary of the rbiom object.
#'     biom
#'     
#'     
#'     # Demonstrate a few accessors.
#'     biom$n_samples
#'     biom$fields
#'     biom$metadata
#'     
#'     
#'     # Edit the metadata table.
#'     biom$metadata$rand <- sample(1:50)
#'     biom %<>% mutate(Obese = BMI >= 30, Sex = NULL)
#'     biom %<>% rename('Years Old' = "Age")
#'     biom$metadata
#'     
#'     
#'     # Subset the rbiom object
#'     biom %<>% subset(`Body Site` == "Saliva" & !Obese)
#'     biom$metadata
#'     
#'     
#'     # Rarefy to an even sampling depth
#'     sample_sums(biom)
#'     
#'     biom %<>% rarefy()
#'     sample_sums(biom)
#' 
NULL

rbiom <- R6::R6Class(
  
  classname  = "rbiom", 
  parent_env = rlang::ns_env('rbiom'), 
  lock_class = TRUE,
  
  private = list(
    
    .counts       = NULL,
    .metadata     = NULL,
    .taxonomy     = NULL,
    .tree         = NULL,
    .sequences    = NULL,
    .depth        = NULL,
    .id           = "", 
    .comment      = "", 
    .date         = "",
    .generated_by = "",
    
    .pkg_version  = packageVersion("rbiom"),
    .hashes       = list(
      .counts       = rlang::hash(NULL),
      .metadata     = rlang::hash(NULL),
      .taxonomy     = rlang::hash(NULL),
      .tree         = rlang::hash(NULL),
      .sequences    = rlang::hash(NULL) ),
    
    
    
    
    #________________________________________________________
    # Intersect OTUs/Samples across private variables.
    # 
    # OTUs are ordered by overall abundance.
    # Sample order is taken from metadata.
    #________________________________________________________
    sync = function (sids = NULL) {
      
      #________________________________________________________
      # Sync sample names between counts and metadata.
      #________________________________________________________
      if (is.null(sids))
        sids <- private$.metadata[['.sample']]
      
      sids <- sids %>%
        intersect(private$.metadata[['.sample']]) %>%
        intersect(colnames(private$.counts))
      
      private$.counts              <- private$.counts[,sids]
      private$.hashes[['.counts']] <- rlang::hash(private$.counts)
      
      if (!identical(sids, private$.metadata[['.sample']]))
        private$.metadata <- private$.metadata[match(sids, private$.metadata[['.sample']]),]
      
      # Drop missing factor levels
      for (i in which(sapply(private$.metadata, is.factor)))
        if (!all(levels(private$.metadata[[i]]) %in% private$.metadata[[i]]))
          private$.metadata[[i]] %<>% {factor(., levels = intersect(levels(.), .))}
      
      private$.hashes[['.metadata']] <- rlang::hash(private$.metadata)
      
      
      #________________________________________________________
      # Sync OTU names between counts, taxonomy, and seqs.
      #________________________________________________________
      otus <- rev(sort(slam::row_sums(private$.counts)))
      otus <- names(otus[otus > 0])
      
      if (!identical(otus, rownames(private$.counts))) {
        private$.counts              <- private$.counts[otus,]
        private$.hashes[['.counts']] <- rlang::hash(private$.counts)
      }
      
      if (!identical(otus, private$.taxonomy[['.otu']])) {
        private$.taxonomy <- private$.taxonomy[match(otus, private$.taxonomy[['.otu']]),]
        private$.hashes[['.taxonomy']] <- rlang::hash(private$.taxonomy)
      }
      
      if (!is.null(private$.sequences))
        if (!identical(names(private$.sequences), otus)) {
          private$.sequences              <- private$.sequences[otus]
          private$.hashes[['.sequences']] <- rlang::hash(private$.sequences)
        }
      
      
      return (invisible(private))
    }
    
    
  ), # /private
  
  
  public = list(
    
    #________________________________________________________
    # Instantiate a new object via `rbiom$new(...)`
    #________________________________________________________
    initialize = function(
    counts, metadata = NULL, taxonomy = NULL, tree = NULL, sequences = NULL, id = "", comment = "", ...) {
      
      if (!is.null(private$.counts))
        cli_abort(c('x' = "Cannot re-initialize an rbiom object."))
      
      dots <- list(...)
      
      
      #________________________________________________________
      # Use active bindings to load components.
      #________________________________________________________
      self$counts <- counts
      
      if (!is.null(metadata))      self$metadata  <- metadata
      if (!is.null(taxonomy))      self$taxonomy  <- taxonomy
      if (!is.null(tree))          self$tree      <- tree
      if (!is.null(sequences))     self$sequences <- sequences
      if (!is.null(id))            self$id        <- id
      if (isTRUE(nzchar(comment))) self$comment   <- comment
      
      
      #________________________________________________________
      # Date found by read_biom(), or current date/time.
      #________________________________________________________
      private$.date <- strftime(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC")
      if (!is.null(date <- dots[['date']]))
        if (is_scalar_character(date) && !is.na(date))
          private$.date <- date
      
      #________________________________________________________
      # Value found by read_biom(). Read-only hereafter.
      #________________________________________________________
      private$.generated_by <- paste("rbiom", packageVersion("rbiom"))
      if (!is.null(generated_by <- dots[['generated_by']]))
        if (is_scalar_character(generated_by) && !is.na(generated_by))
          private$.generated_by <- generated_by
      
    },
    
    
    
    #________________________________________________________
    # Display a summary of an rbiom dataset.
    #________________________________________________________
    print = function (...) {
      
      sids   <- colnames(private$.counts)
      otus   <- rownames(private$.counts)
      ranks  <- colnames(private$.taxonomy)
      fields <- colnames(private$.metadata)
      
      cat(ifelse(nzchar(private$.id), private$.id, "An rbiom object"))
      cat(sprintf(" (%s)\n", substr(private$.date, 1, 10)))
      
      if (nzchar(private$.comment)) {
        width         <- min(80, floor(getOption("width") * 0.75))
        manual_breaks <- strsplit(babies$comment, "\n")[[1]]
        cat(paste(collapse = "\n", strwrap(manual_breaks, width)), "\n")
        cat("-----------\n")
      }
      
      width <- min(80, floor(getOption("width") * 0.75)) - 20
      cat(sprintf("%7.0f Samples:  %s\n", length(sids),   vw(sids,   width)))
      cat(sprintf("%7.0f OTUs:     %s\n", length(otus),   vw(otus,   width)))
      cat(sprintf("%7.0f Ranks:    %s\n", length(ranks),  vw(ranks,  width)))
      cat(sprintf("%7.0f Metadata: %s\n", length(fields), vw(fields, width)))
      cat(sprintf("        Tree:     <%s>\n", ifelse(is.null(private$.tree), "absent", "present")))
      
    }
    
  ), # /public
  
  
  active = list(
    
    
    #________________________________________________________
    # Get or set an rbiom object's arbitrary ID.
    #________________________________________________________
    id = function (value) {
      
      if (missing(value)) return (private$.id)
      
      if (is.null(value) || is.na(value)) value <- ""
      if (!is_scalar_character(value)) cli_abort("Invalid `id`: {.val {value}}")
      private$.id <- value
      
      return (self)
    },
    
    
    
    #________________________________________________________
    # Get or set an rbiom object's arbitrary comment.
    #________________________________________________________
    comment = function (value) {
      
      if (missing(value)) return (private$.comment)
      
      if (is.null(value)) value <- ""
      stopifnot(is_scalar_character(value))
      private$.comment <- value
      
      return (self)
    },
    
    
    
    #________________________________________________________
    # Get a 128-bit hash of this object (string).
    #________________________________________________________
    hash = function () {
      
      rlang::hash(list(
        private$.hashes$.counts,
        private$.hashes$.metadata,
        private$.hashes$.taxonomy,
        private$.hashes$.tree,
        private$.hashes$.sequences,
        private$.depth, 
        private$.id, 
        private$.comment, 
        private$.date, 
        private$.generated_by ))
    },
    
    
    
    #________________________________________________________
    # Get or set the nucleotide sequences for each OTU. 
    # Format: `c('OTU_NAME' = "ATCTAGTA", ...)` or `NULL`.
    #________________________________________________________
    sequences = function (value) {
      
      if (missing(value)) return (private$.sequences)
      
      
      if (!is_null(value) && !(is_character(value) && !is.null(names(value))))
        cli_abort(c(x = "`value` must be NULL or a named character vector, not {.type {value}}."))
      
      
      if (is_character(value)) {
        
        rownames(private$.counts) %<>% trimws()
        otus <- rownames(private$.counts)
        
        if (length(missing <- setdiff(otus, names(value))) != 0)
          cli_abort(c(x = "No sequences given for OTUs: {missing}."))
        
        value <- value[otus]
      }
      
      
      private$.sequences              <- value
      private$.hashes[['.sequences']] <- rlang::hash(private$.sequences)
      
      return (self)
    },
    
    
    
    #________________________________________________________
    # Get or set the phylogenetic tree.
    #________________________________________________________
    tree = function (value) {
      
      if (missing(value)) {
        
        # Subsetting a tree is slow. Only do it when required.
        if (!is.null(private$.tree)) {
          otus <- private$.counts$dimnames[[1]]
          tips <- private$.tree$tip.label
          if (!all(tips %in% otus))
            private$.tree <- tree_subset(private$.tree, otus)
        }
        
        return (private$.tree)
      }
      
      
      if (!is.null(value)) {
        
        validate_tree("value")
        
        otus <- private$.counts$dimnames[[1]]
        tips <- value$tip.label <- trimws(value$tip.label)
        if (length(missing <- setdiff(otus, tips)) > 0)
          cli_abort(('x' = sprintf("OTUs missing from tree: {missing}")))
      }
      
      private$.tree              <- value
      private$.hashes[['.tree']] <- rlang::hash(private$.tree)
      
      return (self)
    },
    
    
    
    #________________________________________________________
    # Get or set the OTU taxonomy.
    #________________________________________________________
    taxonomy = function (value) {
      
      if (missing(value)) return (private$.taxonomy)
      
      
      #________________________________________________________
      # Convert data.frame/matrix to tibble.
      #________________________________________________________
      if (is.null(value)) value <- tibble(.otu = rownames(private$.counts))
      value <- tibble::as_tibble(value, rownames = NA, .name_repair = trimws)
      
      
      #________________________________________________________
      # Make '.otu' the first column.
      #________________________________________________________
      if (!has_rownames(value) && !hasName(value, '.otu'))
        cli_abort(c(x = "Taxonomy table must have row names or an '.otu' column."))
      
      if (has_rownames(value) && hasName(value, '.otu'))
        cli_abort(c(x = "Row names are not allowed when an '.otu' column is present."))
      
      if (has_rownames(value)) { value %<>% rownames_to_column('.otu')
      } else                   { value %<>% relocate(.otu) }
      
      
      #________________________________________________________
      # Disallow column names starting with a "."
      #________________________________________________________
      if (length(badnames <- grep("^\\.", colnames(value)[-1], value = TRUE)) != 0)
        cli_abort(c(x = "Taxonomy rank{?s} can't start with a '.': {badnames}."))
      
      
      #________________________________________________________
      # Convert all cols (except .otu) to factor.
      #________________________________________________________
      for (i in seq_len(ncol(value))[-1]) {
        
        if (!is.factor(value[[i]])) value[[i]] %<>% as.factor()
        
        lvls <- levels(value[[i]])
        if (length(i <- which(lvls != trimws(lvls))) > 0) {
          cli_warn("Whitespace trimmed from {.field {i}} column value{?s} {.val {lvls[i]}}.")
          levels(value[[i]]) %<>% trimws()
        }
      }
      value[['.otu']] %<>% as.character() %>% trimws()
      
      
      #________________________________________________________
      # Match OTU ordering in counts.
      #________________________________________________________
      otus <- rownames(private$.counts)
      if (length(missing <- setdiff(otus, value[['.otu']])) != 0)
        cli_abort(c(x = "OTUs missing from new taxonomy table: {missing}."))
      value <- value[match(otus, value[['.otu']]),]
      
      
      private$.taxonomy              <- value
      private$.hashes[['.taxonomy']] <- rlang::hash(private$.taxonomy)
      
      return (self)
    },
    
    
    
    #________________________________________________________
    # Get or set the sample metadata.
    #________________________________________________________
    metadata = function (value) {
      
      if (missing(value)) return (private$.metadata)
      
      
      #________________________________________________________
      # User wants to erase all the metadata.
      #________________________________________________________
      if (is.null(value)) value <- tibble(.sample = colnames(private$.counts))
      
      #________________________________________________________
      # They're handing us a filename.
      #________________________________________________________
      if (is.character(value) && length(value) == 1)
        value <- local({
        
          fp <- normalizePath(value, mustWork = FALSE, winslash = "/")
          if (!file.exists(fp)) cli_abort("File doesn't exist: {fp}")
          
          
          if (grepl("\\.(xls|xlsx)$", tolower(fp))) {
            
            value <- openxlsx::read.xlsx(fp, sheet = 1)
            
          } else {
            
            #________________________________________________________
            # Auto-detect the field separator.
            #________________________________________________________
            ch  <- c(strsplit(readLines(con = fp, n = 1), '')[[1]], "\t")
            sep <- ch[[head(which(ch %in% c("\t", ",", ";")), 1)]]
            
            value <- read.delim(
              file        = value, 
              sep         = sep, 
              check.names = FALSE,
              strip.white = TRUE )
          }
          
          
          #________________________________________________________
          # Auto-detect the column with sample names.
          #________________________________________________________
          if (!any(colnames(value) == ".sample"))
            for (i in seq_len(ncol(value)))
              if (!anyDuplicated(vals <- value[[i]]))
                if (any(as.character(vals) %in% self$samples)) {
                  rownames(value) <- as.character(vals)
                  value <- value[,-i,drop=FALSE]
                  break
                }
          
          return (value)
        })
      
      
      #________________________________________________________
      # Convert data.frame/matrix to tibble.
      #________________________________________________________
      value <- tibble::as_tibble(value, rownames = NA, .name_repair = trimws)
      
      
      #________________________________________________________
      # Make '.sample' the first column.
      #________________________________________________________
      if (!has_rownames(value) && !hasName(value, '.sample'))
        cli_abort(c(x = "Metadata table must have row names or a '.sample' column."))
      
      if (has_rownames(value) && hasName(value, '.sample'))
        cli_abort(c(x = "Row names are not allowed when a '.sample' column is present."))
      
      if (has_rownames(value)) { value %<>% rownames_to_column('.sample')
      } else                   { value %<>% relocate(.sample) }
      
      
      #________________________________________________________
      # Disallow column names starting with a "."
      #________________________________________________________
      if (length(badnames <- grep("^\\.", colnames(value)[-1], value = TRUE)) != 0)
        cli_abort(c(x = "Metadata field{?s} can't start with a '.': {badnames}."))
      
      
      #________________________________________________________
      # Convert all character cols (except .sample) to factor.
      #________________________________________________________
      for (i in colnames(value)[-1]) {
        
        if (is.character(value[[i]])) value[[i]] %<>% as.factor()
        
        if (is.factor(value[[i]])) {
          lvls <- levels(value[[i]])
          if (length(i <- which(lvls != trimws(lvls))) > 0) {
            cli_warn("Whitespace trimmed from {.field {i}} column value{?s} {.val {lvls[i]}}.")
            levels(value[[i]]) %<>% trimws()
          }
        }
      }
      
      value[['.sample']] %<>% as.character() %>% trimws()
      
      
      
      #________________________________________________________
      # Enforce unnamed vectors.
      #________________________________________________________
      for (i in colnames(value))
        if (!is.null(names(value[[i]])))
          value[[i]] %<>% unname()
      
      
      
      #________________________________________________________
      # Prevent dropping all samples.
      #________________________________________________________
      if (!any(colnames(private$.counts) %in% value[['.sample']]))
        cli_abort(c(
          'i' = "No matching sample names.", 
          '*' = "Expected: {vw(colnames(private$.counts))}", 
          '*' = "Provided: {vw(value[['.sample']])}", 
          'x' = "Can't subset to zero samples." ))
      
      
      private$.metadata              <- value
      private$.hashes[['.metadata']] <- rlang::hash(private$.metadata)
      
      private$sync()
      
      
      return (self)
    },
    
    
    
    #________________________________________________________
    # Get or set the abundance counts.
    #________________________________________________________
    counts = function (value) {
      
      if (missing(value)) return (private$.counts)
      
      
      #________________________________________________________
      # Coerce value to slam matrix.
      #________________________________________________________
      value <- slam::as.simple_triplet_matrix(value)
      stopifnot(slam::is.simple_triplet_matrix(value))
      
      
      #________________________________________________________
      # Sanity check the new matrix.
      #________________________________________________________
      if (!is.numeric(value)) cli_abort("Matrix must be numeric, not {.type {value$v}}.")
      if (any(value$v < 0))   cli_abort("Matrix can't have negative values.")
      if (sum(value) == 0)    cli_abort("Matrix is all zeroes.")
      
      sids <- colnames(value)
      otus <- rownames(value)
      if (!is.character(sids)) cli_abort("Column names must be sample names, not {.type {sids}}.")
      if (!is.character(otus)) cli_abort("Row names must be OTU names, not {.type {otus}}.")
      
      sids <- trimws(sids)
      otus <- trimws(otus)
      if (length(i <- which(sids != colnames(value))) > 0)
        cli_warn("Trimming whitespace from sample name{?s}: {.val {colnames(value)[i]}}.")
      if (length(i <- which(otus != rownames(value))) > 0)
        cli_warn("Trimming whitespace from OTU name{?s}: {.val {rownames(value)[i]}}.")
      
      if (length(dups <- sids[duplicated(sids)]) > 0) cli_abort("Duplicated sample name{?s}: {dups}.")
      if (length(dups <- otus[duplicated(otus)]) > 0) cli_abort("Duplicated OTU name{?s}: {dups}.")
      colnames(value) <- sids
      rownames(value) <- otus
      
      
      #________________________________________________________
      # First time loading counts.
      #________________________________________________________
      if (is.null(private$.counts)) {
        
        private$.counts   <- value
        private$.metadata <- tibble::tibble(.sample = colnames(value))
        private$.taxonomy <- tibble::tibble(.otu    = rownames(value))
        
        private$.hashes[['.counts']]    <- rlang::hash(private$.counts)
        private$.hashes[['.metadata']]  <- rlang::hash(private$.metadata)
        private$.hashes[['.taxonomy']]  <- rlang::hash(private$.taxonomy)
      }
      
      
      #________________________________________________________
      # Intersect with current row/col names.
      #________________________________________________________
      else {
        sids <- intersect(colnames(value), colnames(private$.counts))
        otus <- intersect(rownames(value), rownames(private$.counts))
        
        if (length(sids) == 0)
          cli_abort(c(
            'i' = "Column names don't match any sample names.", 
            '*' = "Expected: {vw(colnames(private$.counts))}", 
            '*' = "Provided: {vw(colnames(value))}", 
            'x' = "Can't subset to zero samples." ))
        
        if (length(otus) == 0)
          cli_abort(c(
            'i' = "Row names don't match any OTU names.", 
            '*' = "Expected: {vw(rownames(private$.counts))}", 
            '*' = "Provided: {vw(rownames(value))}", 
            'x' = "Can't subset to zero OTUs." ))
        
        value <- value[otus,sids]
      }
      
      
      #________________________________________________________
      # Drop zero abundance samples/otus.
      #________________________________________________________
      sids <- names(which(slam::col_sums(value) > 0))
      otus <- names(which(slam::row_sums(value) > 0))
      
      if (length(sids) == 0)
        cli_abort(c(
          'i' = "Matrix is all zeros.", 
          'x' = "Can't subset to zero samples/OTUs." ))
      
      private$.counts              <- value[otus, sids]
      private$.hashes[['.counts']] <- rlang::hash(private$.counts)
      
      
      #________________________________________________________
      # Note the rarefaction depth, or NULL if unrarefied.
      #________________________________________________________
      ss <- unique(round(col_sums(private$.counts), 2))
      private$.depth <- if (length(ss) == 1) ss else NULL
      
      
      private$sync(sids = sids)
      
      return (self)
    },
    
    
    
    #________________________________________________________
    # Get or set the date/time.
    #________________________________________________________
    date = function (value) {
      
      if (missing(value)) return (private$.date)
      
      if (length(value) != 1) cli_abort("Length must be 1.")
      
      fmt <- "%Y-%m-%dT%H:%M:%SZ"
      
      posix <- try(silent = TRUE, {
        if (is.character(value)) { strptime(value, fmt, tz = "UTC") 
        } else                   { as.POSIXlt(value, tz = "UTC")    } })
      
      if (!inherits(posix, 'POSIXt') || is.na(posix))
        cli_abort(c(
          x = "Can't parse date given as: {value}.",
          i = "Expected POSIXt object or string in {.code {fmt}} format."))
      
      
      private$.date <- strftime(posix, fmt, tz="UTC")
      
      return (self)
    },
    
    
    
    #________________________________________________________
    # Get or set sample names.
    #________________________________________________________
    samples = function (value) {
      
      if (missing(value)) return (private$.counts$dimnames[[2]])
      
      stopifnot(is.character(value))
      stopifnot(!anyDuplicated(value))
      stopifnot(!anyNA(value))
      stopifnot(is.null(names(value)))
      stopifnot(length(value) == self$n_samples)
      
      prev <- private$.counts$dimnames[[2]]
      map  <- function (x) value[match(x, prev)]
      
      private$.counts$dimnames[[2]]  %<>% map()
      private$.metadata[['.sample']] %<>% map()
      
      private$.hashes[['.counts']]   <- rlang::hash(private$.counts)
      private$.hashes[['.metadata']] <- rlang::hash(private$.metadata)
      
      return (self)
    },
    
    
    
    #________________________________________________________
    # Get or set OTU names.
    #________________________________________________________
    otus = function (value) {
      
      if (missing(value)) return (private$.counts$dimnames[[1]])
      
      stopifnot(is.character(value))
      stopifnot(!anyDuplicated(value))
      stopifnot(!anyNA(value))
      stopifnot(is.null(names(value)))
      stopifnot(length(value) == self$n_otus)
      
      prev <- private$.counts$dimnames[[1]]
      map  <- function (x) value[match(x, prev)]
      
      private$.counts$dimnames[[1]] %<>% map()
      private$.taxonomy[['.otu']]   %<>% map()
      if (!is.null(private$.sequences)) names(private$.sequences) %<>% map()
      if (!is.null(private$.tree))      private$.tree$tip.label   %<>% map()
      
      private$.hashes[['.counts']]    <- rlang::hash(private$.counts)
      private$.hashes[['.taxonomy']]  <- rlang::hash(private$.taxonomy)
      private$.hashes[['.sequences']] <- rlang::hash(private$.sequences)
      private$.hashes[['.tree']]      <- rlang::hash(private$.tree)
      
      return (self)
    },
    
    
    
    #________________________________________________________
    # Read-only accessors.
    #________________________________________________________
    fields       = function () { attr(private$.metadata, 'names')         },
    ranks        = function () { attr(private$.taxonomy, 'names')         },
    n_otus       = function () { private$.counts$nrow                     },
    n_samples    = function () { private$.counts$ncol                     },
    n_fields     = function () { length(attr(private$.metadata, 'names')) },
    n_ranks      = function () { length(attr(private$.taxonomy, 'names')) },
    generated_by = function () { private$.generated_by                    },
    pkg_version  = function () { private$.pkg_version                     },
    depth        = function () { private$.depth                           }
    
  ) # /active
)

