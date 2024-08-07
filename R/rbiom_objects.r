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
#' | **Accessor**             | **Content**                                            |
#' | ------------------------ | ------------------------------------------------------ |
#' | `$counts`                | Abundance of each OTU in each sample.                  |
#' | `$metadata`              | Sample mappings to metadata (treatment, patient, etc). |
#' | `$taxonomy`              | OTU mappings to taxonomic ranks (genus, phylum, etc).  |
#' | `$otus`, `$n_otus`       | OTU names.                                             |
#' | `$samples`, `$n_samples` | Sample names.                                          |
#' | `$fields`, `$n_fields`   | Metadata field names.                                  |
#' | `$ranks`, `$n_ranks`     | Taxonomic rank names.                                  |
#' | `$tree`, `$sequences`    | Phylogenetic tree / sequences for the OTUs, or `NULL`. |
#' | `$id`, `$comment`        | Arbitrary strings for describing the dataset.          |
#' | `$depth`                 | Rarefaction depth, or `NULL` if unrarefied.            |
#' | `$date`                  | Date from BIOM file.                                   |
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
#' | `$counts`         | Matrix of abundances; OTUs (rows) by samples (columns). |
#' | `$metadata`       | Data.frame with `'.sample'` column, or a file name.     |
#' | `$taxonomy`       | Data.frame with `'.otu'` as the first column.           |
#' | `$otus`           | Character vector with new names for the OTUs.           |
#' | `$samples`        | Character vector with new names for the samples.        |
#' | `$tree`           | Phylo object with the phylogenetic tree for the OTUs.   |
#' | `$sequences`      | Named character vector of OTU reference sequences.      |
#' | `$id`, `$comment` | String with dataset's title or comment.                 |
#' | `$date`           | Date-like object, or `"%Y-%m-%dT%H:%M:%SZ"` string.     |
#' 
#' \cr
#' 
#' 
#' @section Transformations:
#' 
#' All functions return an rbiom object.
#' 
#' | **Function**                     | **Transformation**                               |
#' | -------------------------------- | ------------------------------------------------ |
#' | \code{<rbiom>$clone()}           | Safely duplicate an rbiom object.                |
#' | \code{\link[=subset]{<rbiom>[}}  | Subset to a specific set of sample names.        |
#' | [subset()]                       | Subset samples according to metadata properties. |
#' | [slice()]                        | Subset to a specific number of samples.          |
#' | [mutate()]                       | Create, modify, and delete metadata fields.      |
#' | [rarefy()]                       | Sub-sample OTU counts to an even sampling depth. |
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
#'     # Markdown syntax for comments is recommended.
#'     biom$comment %>% cli::cli_text()
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
    .generated_by = paste("rbiom", utils::packageVersion("rbiom")),
    
    .pkg_version  = utils::packageVersion("rbiom"),
    .hash         = NULL,
    .tree_hash    = NULL,
    
    
    
    
    #________________________________________________________
    # Intersect OTUs/Samples across private variables.
    # 
    # OTUs are ordered by overall abundance.
    # Sample order is taken from metadata.
    #________________________________________________________
    sync = function (sids = NULL, otus = NULL) {
      
      mtx <- private$.counts
      
      #________________________________________________________
      # Ordering of sids/otus when not explicitly provided.
      #________________________________________________________
      if (is.null(sids)) sids <- intersect(colnames(mtx), private$.metadata[['.sample']])
      if (is.null(otus)) otus <- intersect(rownames(mtx), private$.taxonomy[['.otu']])
      
      
      #________________________________________________________
      # Maintain the $counts matrix.
      #________________________________________________________
      if (!identical(sids, colnames(mtx)) || !identical(otus, rownames(mtx)))
        mtx <- mtx[otus,sids]
      
      
      #________________________________________________________
      # Drop zero abundance otus/samples.
      #________________________________________________________
      i <- sort(unique(mtx$i))
      j <- sort(unique(mtx$j))
      if (!length(c(i, j)) == nrow(mtx) + ncol(mtx)) {
        mtx  <- mtx[i,j]
        otus <- otus[i]
        sids <- sids[j]
      }
      
      
      #________________________________________________________
      # Update any components that are out of sync.
      #________________________________________________________
      if (!identical(sids, private$.metadata[['.sample']]))
        private$.metadata <- private$.metadata[match(sids, private$.metadata[['.sample']]),]
      
      if (!identical(otus, private$.taxonomy[['.otu']]))
        private$.taxonomy <- private$.taxonomy[match(otus, private$.taxonomy[['.otu']]),]
      
      if (!is.null(private$.sequences))
        if (!identical(names(private$.sequences), otus))
          private$.sequences <- private$.sequences[otus]
      
      
      private$.hash     <- NULL
      private$.counts   <- mtx
      private$.metadata <- relevel(private$.metadata)
      private$.taxonomy <- relevel(private$.taxonomy)
      
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
      # Read-only hereafter.
      #________________________________________________________
      for (i in c('generated_by', 'pkg_version'))
        if (is_scalar_character(dots[[i]]) && !is.na(dots[[i]]))
          private[[paste0('.', i)]] <- dots[[i]]
      
    },
    
    
    
    #________________________________________________________
    # Display a summary of an rbiom dataset.
    #________________________________________________________
    print = function (...) {
      
      sids   <- colnames(private$.counts)
      otus   <- rownames(private$.counts)
      ranks  <- colnames(private$.taxonomy)
      fields <- colnames(private$.metadata)
      
      width <- min(80L, as.integer(getOption("width") * 0.75))
      rlang::local_options(cli.width = width)
      
      
      div <- cli::cli_div(theme = list(
        'rule' = list(color = "grey50", 'line-type' = "double"),
        'span' = list(color = "white") ))
      
      cli::cli_text("\n")
      title <- ifelse(nzchar(private$.id), private$.id, "An rbiom object")
      cli::cli_rule(paste0("{.span ", title, "}"), .envir = stop())
      
      cli::cli_end(div)
      cli::cli_text("\n")
      
      
      div <- cli::cli_div(theme = list(
        'span.url' = list(color = "grey50"),
        'rule'     = list(color = "grey50"),
        'span.q'   = list(fmt = function (x) gsub(" ", "\u00a0", x)) ))
      
      if (nzchar(private$.comment)) {
        private$.comment %>%    # Make clickable hyperlinks.
          gsub('(\\[.+?\\]\\(.+?\\))', '{.href \\1}', .) %>%
          cli::cli_text(.envir = stop())
        cli::cli_text("\n")
      }
      
      cli::cli_text("{.q {sprintf('%7.0f Samples: %s\n', length(sids),   vw(sids,   width - 20))}}")
      cli::cli_text("{.q {sprintf('%7.0f OTUs:    %s\n', length(otus),   vw(otus,   width - 20))}}")
      cli::cli_text("{.q {sprintf('%7.0f Ranks:   %s\n', length(ranks),  vw(ranks,  width - 20))}}")
      cli::cli_text("{.q {sprintf('%7.0f Fields:  %s\n', length(fields), vw(fields, width - 20))}}")
      cli::cli_text("{.q {'        Tree:    '}}{.emph {ifelse(is.null(private$.tree), '<absent>', '<present>')}}")
      cli::cli_text("\n")
      
      if (is.null(depth <- private$.depth))
        depth <- range(col_sums(private$.counts)) %>%
          {scales::label_number(scale_cut=scales::cut_si(''))}() %>% 
          sub(' ', '', .) %>% 
          paste(collapse = " - ")
      
      cli::cli_rule(
        left   = paste(depth, 'reads/sample'), 
        right  = substr(private$.date, 1, 10) )
      
      cli::cli_end(div)
      cli::cli_text("\n")
    }
    
  ), # /public
  
  
  active = list(
    
    
    #________________________________________________________
    # Get or set an rbiom object's arbitrary ID.
    #________________________________________________________
    id = function (value) {
      
      if (missing(value)) return (private$.id)
      
      if (is.null(value) || is.na(value)) value <- "Untitled Dataset"
      if (!is_scalar_character(value)) cli_abort("Invalid `id`: {.type {value}}")
      value <- trimws(value, whitespace = "[\\h\\v]")
      if (nchar(value) < 1)   value <- "Untitled Dataset"
      if (nchar(value) > 100) {
        if (interactive()) cli_warn("Truncating `id` to 100 characters.")
        value <- substr(value, 0, 100)
      }
      
      private$.hash <- NULL
      private$.id   <- value
      
      return (self)
    },
    
    
    
    #________________________________________________________
    # Get or set an rbiom object's arbitrary comment.
    #________________________________________________________
    comment = function (value) {
      
      if (missing(value)) return (private$.comment)
      
      if (is.null(value) || is.na(value)) value <- ""
      if (!is_scalar_character(value)) cli_abort("Invalid `comment`: {.type {value}}")
      value <- trimws(value, whitespace = "[\\h\\v]")
      if (nchar(value) > 5000) {
        if (interactive()) cli_warn("Truncating `comment` to 5000 characters.")
        value <- substr(value, 0, 5000)
      }
      
      private$.hash    <- NULL
      private$.comment <- value
      
      return (self)
    },
    
    
    
    # #________________________________________________________
    # # Get or set an rbiom object's comment in HTML format.
    # #________________________________________________________
    # comment_html = function (value) {
    #   
    #   if (!nzchar(system.file(package = "pandoc")))
    #     cli_abort("`$comment_html` requires the {.pkg pandoc} R package.")
    #   
    #   # Return existing markdown comment, converted to html.
    #   if (missing(value)) {
    #     md   <- private$.comment
    #     html <- pandoc::pandoc_convert(text = md, from = 'markdown', to = 'html')
    #     html <- paste0(collapse = '\n', html)
    #     return (html)
    #   }
    #   
    #   # Sanity check incoming html string for type and length.
    #   if (is.null(value) || is.na(value)) value <- ""
    #   if (!is_scalar_character(value)) cli_abort("Invalid `comment_html`: {.type {value}}")
    #   if (nchar(value) > 5000) cli_abort("Invalid `comment_html`: maximum length is 5000 characters.")
    #   
    #   # Save the provided html comment in markdown format.
    #   html <- value
    #   md   <- pandoc::pandoc_convert(text = html, from = 'html', to = 'markdown')
    #   md   <- paste0(collapse = '\n', md)
    #   private$.comment <- md
    #   
    #   return (self)
    # },
    
    
    
    #________________________________________________________
    # Get a 128-bit hash of this object (string).
    #________________________________________________________
    hash = function () {
      
      if (is.null(private$.hash)) {
        
        private$.hash <- rlang::hash(list(
          private$.counts,
          private$.metadata,
          private$.taxonomy,
          private$.sequences,
          private$.tree_hash, 
          private$.depth, 
          private$.id, 
          private$.comment, 
          private$.date, 
          private$.generated_by ))
      }
      
      return (private$.hash)
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
        
        names(value) %<>% trimws()
        otus <- rownames(private$.counts)
        
        if (length(missing <- setdiff(otus, names(value))) != 0)
          cli_abort(c(x = "No sequences given for OTUs: {missing}."))
        
        value <- value[otus]
      }
      
      
      private$.hash      <-  NULL
      private$.sequences <- value
      
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
      
      private$.hash      <- NULL
      private$.tree      <- value
      private$.tree_hash <- rlang::hash(value)
      
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
      if (!tibble::has_rownames(value) && !hasName(value, '.otu'))
        cli_abort(c(x = "Taxonomy table must have row names or an '.otu' column."))
      
      if (tibble::has_rownames(value) && hasName(value, '.otu'))
        cli_abort(c(x = "Row names are not allowed when an '.otu' column is present."))
      
      if (tibble::has_rownames(value)) { value %<>% tibble::rownames_to_column('.otu')
      } else                           { value %<>% relocate(.otu) }
      
      
      #________________________________________________________
      # Disallow column names starting with a "."
      #________________________________________________________
      if (length(i <- grep("^\\.", colnames(value)[-1], value = TRUE)) != 0)
        cli_abort(c(x = "Taxonomy rank{?s} can't start with a '.': {.val {i}}."))
      
      
      #________________________________________________________
      # Convert all cols (except .otu) to factor.
      #________________________________________________________
      for (i in seq_len(ncol(value))[-1])
        if (!is.factor(value[[i]])) value[[i]] %<>% as.factor()
      
      
      #________________________________________________________
      # Intersect with OTU names from $counts.
      #________________________________________________________
      value[['.otu']] %<>% as.character() %>% trimws()
      otus <- intersect(value[['.otu']], rownames(private$.counts))
      
      if (length(otus) == 0) cli_abort(c(
        'i' = "No matching OTU names.", 
        '*' = "Expected: {vw(rownames(private$.counts))}", 
        '*' = "Provided: {vw(value[['.otu']])}", 
        'x' = "Can't subset to zero OTUs." ))
      
      if (length(i <- setdiff(rownames(private$.counts), otus)) != 0)
        cli_warn(c('i' = paste(
          "Dropping {length(i)} OTU{?s} from biom object",
          "since they are not in the new taxonomy: {.val {i}}." )))
      
      if (length(i <- setdiff(value[['.otu']], otus)) != 0)
        cli_warn(c('i' = paste(
          "Ignoring taxonomy for {length(i)} OTU{?s}", 
          "not currently in biom object: {.val {i}}." )))
      
      
      private$.taxonomy <- value
      private$sync(otus = otus)
      
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
            
            if (nchar(system.file(package = "openxlsx")) == 0)
              stop("CRAN R package 'openxlsx' must be installed to import Excel files.")
            
            value <- openxlsx::read.xlsx(fp, sheet = 1)
            
          } else {
            
            #________________________________________________________
            # Auto-detect the field separator.
            #________________________________________________________
            ch  <- c(strsplit(readLines(con = fp, n = 1), '')[[1]], "\t")
            sep <- ch[[head(which(ch %in% c("\t", ",", ";")), 1)]]
            
            value <- utils::read.delim(
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
      if (!tibble::has_rownames(value) && !hasName(value, '.sample'))
        cli_abort(c(x = "Metadata table must have row names or a '.sample' column."))
      
      if (tibble::has_rownames(value) && hasName(value, '.sample'))
        cli_abort(c(x = "Row names are not allowed when a '.sample' column is present."))
      
      if (tibble::has_rownames(value)) { value %<>% tibble::rownames_to_column('.sample')
      } else                           { value %<>% relocate(.sample) }
      
      
      #________________________________________________________
      # Disallow column names starting with a "."
      #________________________________________________________
      if (length(i <- grep("^\\.", colnames(value)[-1], value = TRUE)) != 0)
        cli_abort(c(x = "Metadata field{?s} can't start with a '.': {.val {i}}."))
      
      
      #________________________________________________________
      # Convert all character cols (except .sample) to factor.
      #________________________________________________________
      for (i in colnames(value)[-1])
        if (is.character(value[[i]])) value[[i]] %<>% as.factor()
      
      
      #________________________________________________________
      # Intersect with sample names from $counts.
      #________________________________________________________
      value[['.sample']] %<>% as.character() %>% trimws()
      sids <- intersect(value[['.sample']], colnames(private$.counts))
      
      if (length(sids) == 0) cli_abort(c(
        'i' = "No matching sample names.", 
        '*' = "Expected: {vw(colnames(private$.counts))}", 
        '*' = "Provided: {vw(value[['.sample']])}", 
        'x' = "Can't subset to zero samples." ))
      
      if (length(i <- setdiff(colnames(private$.counts), sids)) != 0) {
        cli_warn(c('i' = paste(
          "Dropping {length(i)} sample{?s} from biom object",
          "since they are not in the new metadata: {.val {i}}." )))
      }
      
      if (length(i <- setdiff(value[['.sample']], sids)) != 0)
        cli_warn(c('i' = paste(
          "Ignoring metadata for {length(i)} sample{?s}", 
          "not currently in biom object: {.val {i}}." )))
      
      
      private$.metadata <- value
      private$sync(sids = sids)
      
      
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
      mtx <- as.simple_triplet_matrix(value)
      stopifnot(is.simple_triplet_matrix(mtx))
      
      
      #________________________________________________________
      # Sanity check the new matrix.
      #________________________________________________________
      if (length(mtx$v) == 0) cli_abort("Matrix is all zeroes.")
      if (!is.numeric(mtx$v)) cli_abort("Matrix must be numeric, not {.type {mtx$v}}.")
      if (any(mtx$v < 0))     cli_abort("Matrix can't have negative values.")
      
      sids <- colnames(mtx)
      otus <- rownames(mtx)
      if (!is.character(sids)) cli_abort("Column names must be sample names, not {.type {sids}}.")
      if (!is.character(otus)) cli_abort("Row names must be OTU names, not {.type {otus}}.")
      
      sids <- trimws(sids)
      otus <- trimws(otus)
      if (length(i <- which(sids != colnames(mtx))) > 0)
        cli_warn("Trimming whitespace from sample name{?s}: {.val {colnames(mtx)[i]}}.")
      if (length(i <- which(otus != rownames(mtx))) > 0)
        cli_warn("Trimming whitespace from OTU name{?s}: {.val {rownames(mtx)[i]}}.")
      
      if (length(dups <- sids[duplicated(sids)]) > 0) cli_abort("Duplicated sample name{?s}: {dups}.")
      if (length(dups <- otus[duplicated(otus)]) > 0) cli_abort("Duplicated OTU name{?s}: {dups}.")
      colnames(mtx) <- sids
      rownames(mtx) <- otus
      
      
      #________________________________________________________
      # First time loading counts.
      #________________________________________________________
      if (is.null(private$.counts)) {
        private$.metadata <- tibble::tibble(.sample = sids)
        private$.taxonomy <- tibble::tibble(.otu    = otus)
      }
      
      #________________________________________________________
      # Intersect with current row/col names.
      #________________________________________________________
      else {
        sids <- intersect(sids, colnames(private$.counts))
        otus <- intersect(otus, rownames(private$.counts))
        
        if (length(sids) == 0)
          cli_abort(c(
            'i' = "Column names don't match any sample names.", 
            '*' = "Expected: {vw(colnames(private$.counts))}", 
            '*' = "Provided: {vw(sids)}", 
            'x' = "Can't subset to zero samples." ))
        
        if (length(otus) == 0)
          cli_abort(c(
            'i' = "Row names don't match any OTU names.", 
            '*' = "Expected: {vw(rownames(private$.counts))}", 
            '*' = "Provided: {vw(otus)}", 
            'x' = "Can't subset to zero OTUs." ))
        
        mtx <- mtx[otus, sids]
      }
      
      
      #________________________________________________________
      # Note the rarefaction depth, or NULL if unrarefied.
      #________________________________________________________
      ss <- unique(round(col_sums(mtx), 2))
      private$.depth <- if (length(ss) == 1) ss else NULL
      
      
      private$.counts <- mtx
      private$sync()
      
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
      
      
      private$.hash <- NULL
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
      
      private$.hash                   <-  NULL
      private$.counts$dimnames[[2]]  %<>% map()
      private$.metadata[['.sample']] %<>% map()
      
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
      
      private$.hash                  <-  NULL
      private$.counts$dimnames[[1]] %<>% map()
      private$.taxonomy[['.otu']]   %<>% map()
      if (!is.null(private$.sequences)) names(private$.sequences) %<>% map()
      if (!is.null(private$.tree))      private$.tree$tip.label   %<>% map()
      
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

