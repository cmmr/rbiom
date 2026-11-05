#' Working with rbiom Objects.
#' 
#' @description
#' Rbiom objects make it easy to access and manipulate your BIOM data, ensuring 
#' all the disparate components remain in sync. These objects behave largely 
#' like lists, in that you can access and assign to them using the `$` 
#' operator. The sections below list all the fields which can be read and/or 
#' written, and the helper functions for common tasks like rarefying and 
#' subsetting. To create an rbiom object, see [as_rbiom()].
#' 
#' Use `$clone()` to create a copy of an rbiom object. This is necessary 
#' because rbiom objects are **passed by reference**. The usual `<-` assignment 
#' operator will simply create a second reference to the same object - it will 
#' not create a second object. See [speed ups][speed] for more details.
#' 
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
#' | `$relative`              | `TRUE` when counts are relative abundances (sum to 1). |
#' | `$date`                  | Date from BIOM file.                                   |
#' | `$history`               | List of the calls that produced this object.           |
#' | `$code`                  | R code that reproduces this object (see below).        |
#'
#'
#'
#' @section Provenance:
#'
#' Every rbiom object records how it was made. The constructor call
#' ([read_biom()], [as_rbiom()]) and each subsequent transformation
#' ([subset()], [rarefy()], [mutate()], etc.) are appended to `$history`, and
#' `$code` prints them as a runnable script. Tables, statistics, and plots
#' derived from the object inherit this history, so their `$code` is a
#' complete, self-contained recipe.
#'
#' Assigning directly to a field (e.g. `biom$metadata <- df`) is also
#' recorded, but as a comment, since the assigned value only exists in memory.
#' Prefer the transformation functions, or build the object from files, when
#' a fully reproducible script is wanted.
#'
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
    .generated_by = paste("rbiom", packageVersion("rbiom")),
    
    .pkg_version  = packageVersion("rbiom"),
    .hash         = NULL,
    .tree_hash    = NULL,
    .history      = list(),
    .relative     = FALSE,
    underscores   = FALSE,




    #________________________________________________________
    # Intersect OTUs/Samples across private variables.
    #
    # When `sids`/`otus` are given (e.g. after assigning new metadata or
    # taxonomy), that ordering wins. Otherwise the order of the counts
    # matrix is kept.
    #________________________________________________________
    sync = function (sids = NULL, otus = NULL) rb__sync(self, private, sids, otus)
    
    
  ), # /private
  
  
  public = list(
    
    #________________________________________________________
    # Instantiate a new object via `rbiom$new(...)`
    #________________________________________________________
    initialize = function (
      counts, metadata = NULL, taxonomy = NULL, tree = NULL, 
      sequences = NULL, id = "", comment = "", underscores = FALSE, ...) {
      
        rb_init(
        self, private,
        counts, metadata, taxonomy, tree, 
        sequences, id, comment, underscores, ...)
      },
    
    
    #________________________________________________________
    # Display a summary of an rbiom dataset.
    #________________________________________________________
    print = function (...) rb_print(self, private)
    
  ), # /public
  
  
  active = list(
    
    
    #________________________________________________________
    # Get or set an rbiom object's arbitrary ID.
    #________________________________________________________
    id = function (value) rb_id(self, private, value),
    
    
    
    #________________________________________________________
    # Get or set an rbiom object's arbitrary comment.
    #________________________________________________________
    comment = function (value) rb_comment(self, private, value),
    
    
    
    #________________________________________________________
    # Get a 128-bit hash of this object (string).
    #________________________________________________________
    hash = function () rb_hash(self, private),
    
    
    
    #________________________________________________________
    # Get or set the nucleotide sequences for each OTU. 
    # Format: `c('OTU_NAME' = "ATCTAGTA", ...)` or `NULL`.
    #________________________________________________________
    sequences = function (value) rb_sequences(self, private, value),
    
    
    
    #________________________________________________________
    # Get or set the phylogenetic tree.
    #________________________________________________________
    tree = function (value) rb_tree(self, private, value),
    
    
    
    #________________________________________________________
    # Get or set the OTU taxonomy.
    #________________________________________________________
    taxonomy = function (value) rb_taxonomy(self, private, value),
    
    
    
    #________________________________________________________
    # Get or set the sample metadata.
    #________________________________________________________
    metadata = function (value) rb_metadata(self, private, value),
    
    
    
    #________________________________________________________
    # Get or set the abundance counts.
    #________________________________________________________
    counts = function (value) rb_counts(self, private, value),
    
    
    
    #________________________________________________________
    # Get or set the date/time.
    #________________________________________________________
    date = function (value) rb_date(self, private, value),
    
    
    
    #________________________________________________________
    # Get or set sample names.
    #________________________________________________________
    samples = function (value) rb_samples(self, private, value),
    
    
    
    #________________________________________________________
    # Get or set OTU names.
    #________________________________________________________
    otus = function (value) rb_otus(self, private, value),
    
    
    
    #________________________________________________________
    # Read-only accessors.
    #________________________________________________________
    fields       = function () attr(private$.metadata, 'names'),
    ranks        = function () attr(private$.taxonomy, 'names'),
    n_otus       = function () nrow(private$.counts),
    n_samples    = function () ncol(private$.counts),
    n_fields     = function () length(attr(private$.metadata, 'names')),
    n_ranks      = function () length(attr(private$.taxonomy, 'names')),
    generated_by = function () private$.generated_by,
    pkg_version  = function () private$.pkg_version,
    depth        = function () private$.depth,
    relative     = function () isTRUE(private$.relative),


    #________________________________________________________
    # Provenance: the calls that produced this object.
    #________________________________________________________
    history      = function () lapply(private$.history, `[[`, 'call'),
    code         = function () add_class(paste(collapse = "\n", c(
                     sprintf("# Generated by rbiom %s", packageVersion('rbiom')),
                     "library(rbiom)", "", prov_lines(self) )), 'rbiom_code')

  ) # /active
)


rb__sync <- function (self, private, sids, otus) {
  
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
    mtx <- mtx[otus, sids, drop = FALSE]


  #________________________________________________________
  # Drop zero abundance otus/samples.
  #________________________________________________________
  i <- which(rowSums(mtx) > 0)
  j <- which(colSums(mtx) > 0)
  if (length(i) < nrow(mtx) || length(j) < ncol(mtx)) {
    mtx  <- mtx[i, j, drop = FALSE]
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
  
  return (invisible(NULL))
}


rb_init <- function(
    self, private,
    counts, metadata, taxonomy, tree, 
    sequences, id, comment, underscores, ...) {
  
  if (!is.null(private$.counts))
    cli_abort(c('x' = "Cannot re-initialize an rbiom object."))
  
  dots <- list(...)
  
  validate_bool('underscores')
  private$underscores <- underscores

  #________________________________________________________
  # Use active bindings to load components.
  #________________________________________________________
  prov_quiet({

    self$counts <- counts

    if (!is.null(metadata))      self$metadata  <- metadata
    if (!is.null(taxonomy))      self$taxonomy  <- taxonomy
    if (!is.null(tree))          self$tree      <- tree
    if (!is.null(sequences))     self$sequences <- sequences
    if (!is.null(id))            self$id        <- id
    if (isTRUE(nzchar(comment))) self$comment   <- comment
  })

  #________________________________________________________
  # Date found by read_biom(), or current date/time.
  #________________________________________________________
  private$.date <- strftime(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC")
  if (!is.null(date <- dots[['date']]))
    if (is_scalar_character(date) && !anyNA(date))
      private$.date <- date

  #________________________________________________________
  # Read-only hereafter.
  #________________________________________________________
  for (i in c('generated_by', 'pkg_version'))
    if (is_scalar_character(dots[[i]]) && !anyNA(dots[[i]]))
      private[[paste0('.', i)]] <- dots[[i]]

  #________________________________________________________
  # Provenance. read_biom() and as_rbiom() replace this generic
  # constructor entry with the actual call.
  #________________________________________________________
  prov_init(self, prov_call('as_rbiom', list(counts)))

  return (invisible(self))
}


rb_print <- function (self, private) {
  
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
  
  if (isTRUE(private$.relative)) {
    depth <- "relative abundances"

  } else {
    if (is.null(depth <- private$.depth))
      depth <- range(colSums(private$.counts)) %>%
      {scales::label_number(scale_cut=scales::cut_si(''))}() %>%
      sub(' ', '', .) %>%
      paste(collapse = " - ")
    depth <- paste(depth, 'reads/sample')
  }

  cli::cli_rule(
    left   = depth,
    right  = substr(private$.date, 1, 10) )
  
  cli::cli_end(div)
  cli::cli_text("\n")
}



rb_id <- function (self, private, value) {
  
  if (missing(value)) return (private$.id)
  
  if (is.null(value) || anyNA(value)) value <- "Untitled Dataset"
  if (!is_scalar_character(value)) cli_abort("Invalid `id`: {.type {value}}")
  value <- trimws(value, whitespace = "[\\h\\v]")
  if (nchar(value) < 1)   value <- "Untitled Dataset"
  if (nchar(value) > 100) {
    cli_warn("Truncating `id` to 100 characters.")
    value <- substr(value, 0, 100)
  }
  
  private$.hash <- NULL
  private$.id   <- value
  prov_setter(self, 'id', value)
}


rb_comment <- function (self, private, value) {
  
  if (missing(value)) return (private$.comment)
  
  if (is.null(value) || anyNA(value)) value <- ""
  if (!is_scalar_character(value)) cli_abort("Invalid `comment`: {.type {value}}")
  value <- trimws(value, whitespace = "[\\h\\v]")
  if (nchar(value) > 5000) {
    cli_warn("Truncating `comment` to 5000 characters.")
    value <- substr(value, 0, 5000)
  }
  
  private$.hash    <- NULL
  private$.comment <- value
  prov_setter(self, 'comment', value)
}


rb_hash <- function (self, private) {
  
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
}


rb_sequences <- function (self, private, value) {
  
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
  prov_setter(self, 'sequences', value)
}


rb_tree <- function (self, private, value) {
  
  if (missing(value)) {
    
    # Subsetting a tree is slow. Only do it when required.
    if (!is.null(private$.tree)) {
      otus <- private$.counts@Dimnames[[1]]
      tips <- private$.tree$tip.label
      if (!all(tips %in% otus))
        private$.tree <- tree_subset(private$.tree, otus)
    }
    
    return (private$.tree)
  }
  
  
  if (!is.null(value)) {
    
    validate_tree("value", underscores = private$underscores)
    
    otus <- private$.counts@Dimnames[[1]]
    tips <- value$tip.label <- trimws(value$tip.label)
    if (length(missing <- setdiff(otus, tips)) > 0)
      cli_abort(c('x' = sprintf("OTUs missing from tree: {missing}")))
  }
  
  private$.hash      <- NULL
  private$.tree      <- value
  private$.tree_hash <- rlang::hash(value)
  prov_setter(self, 'tree', value)
}


rb_taxonomy <- function (self, private, value) {

  if (missing(value)) return (private$.taxonomy)

  private$.taxonomy <- import_taxonomy(value, otus = rownames(private$.counts))
  private$sync(otus = private$.taxonomy[['.otu']])
  prov_setter(self, 'taxonomy', value)
}


rb_metadata <- function (self, private, value) {

  if (missing(value)) return (private$.metadata)

  private$.metadata <- import_metadata(value, sids = colnames(private$.counts))
  private$sync(sids = private$.metadata[['.sample']])
  prov_setter(self, 'metadata', value)
}


rb_counts <- function (self, private, value) {

  if (missing(value)) return (private$.counts)


  #________________________________________________________
  # Coerce value to a sparse numeric matrix (dgCMatrix).
  #________________________________________________________
  if (is.data.frame(value)) value <- as.matrix(value)

  if (!(is.matrix(value) || inherits(value, 'Matrix')))
    cli_abort("`counts` must be a matrix, not {.type {value}}.")

  if (is.matrix(value) && !(is.numeric(value) || is.logical(value)))
    cli_abort("`counts` must be a numeric matrix, not a {typeof(value)} matrix.")

  # generalMatrix: Matrix returns dtCMatrix for triangular/diagonal/1x1 input.
  mtx <- tryCatch(
    expr  = as(as(as(value, "CsparseMatrix"), "generalMatrix"), "dMatrix"),
    error = function (e) cli_abort(c(
      'x' = "Can't convert `counts` to a numeric matrix.",
      'i' = conditionMessage(e) )))

  stopifnot(inherits(mtx, "dgCMatrix"))


  #________________________________________________________
  # Sanity check the new matrix.
  #________________________________________________________
  if (anyNA(mtx@x))       cli_abort("Matrix can't contain NA values.")
  mtx <- Matrix::drop0(mtx)
  if (length(mtx@x) == 0) cli_abort("Matrix is all zeroes.")
  if (any(mtx@x < 0))     cli_abort("Matrix can't have negative values.")
  
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
    
    mtx <- mtx[otus, sids, drop = FALSE]
  }


  #________________________________________________________
  # Note the rarefaction depth, or NULL if unrarefied. Counts that all sum
  # to 1 and aren't integers are relative abundances, not rarefied counts.
  #________________________________________________________
  ss <- unique(round(colSums(mtx), 2))
  private$.relative <- length(ss) == 1 && isTRUE(all.equal(ss, 1)) && any(mtx@x %% 1 != 0)
  private$.depth    <- if (length(ss) == 1 && !private$.relative) ss else NULL


  private$.counts <- mtx
  private$sync()
  prov_setter(self, 'counts', value)
}


rb_date <- function (self, private, value) {
  
  if (missing(value)) return (private$.date)
  
  if (length(value) != 1) cli_abort("Length must be 1.")
  
  fmt <- "%Y-%m-%dT%H:%M:%SZ"
  
  posix <- try(silent = TRUE, {
    if (is.character(value)) { strptime(value, fmt, tz = "UTC") 
    } else                   { as.POSIXlt(value, tz = "UTC")    } })
  
  if (!inherits(posix, 'POSIXt') || anyNA(posix))
    cli_abort(c(
      x = "Can't parse date given as: {value}.",
      i = "Expected POSIXt object or string in {.code {fmt}} format."))
  
  
  private$.hash <- NULL
  private$.date <- strftime(posix, fmt, tz="UTC")
  prov_setter(self, 'date', private$.date)
}


rb_samples <- function (self, private, value) {
  
  if (missing(value)) return (private$.counts@Dimnames[[2]])

  value <- validate_names(value, self$n_samples, "sample")
  
  prev <- private$.counts@Dimnames[[2]]
  map  <- function (x) value[match(x, prev)]
  
  private$.hash                   <-  NULL
  private$.counts@Dimnames[[2]]  %<>% map()
  private$.metadata[['.sample']] %<>% map()
  prov_setter(self, 'samples', value)
}


rb_otus <- function (self, private, value) {
  
  if (missing(value)) return (private$.counts@Dimnames[[1]])

  value <- validate_names(value, self$n_otus, "OTU")
  
  prev <- private$.counts@Dimnames[[1]]
  map  <- function (x) value[match(x, prev)]
  
  private$.hash                  <-  NULL
  private$.counts@Dimnames[[1]] %<>% map()
  private$.taxonomy[['.otu']]   %<>% map()
  if (!is.null(private$.sequences)) names(private$.sequences) %<>% map()
  if (!is.null(private$.tree))      private$.tree$tip.label   %<>% map()
  prov_setter(self, 'otus', value)
}


#________________________________________________________
# Validate new sample or OTU names. Leading/trailing whitespace is trimmed
# (with a warning), matching how names are handled when counts and 
# metadata are assigned.
#________________________________________________________
validate_names <- function (value, n, what) {
  
  if (!is.character(value))
    cli_abort("New {what} names must be character, not {.type {value}}.")
  if (length(value) != n)
    cli_abort("Expected {n} {what} name{?s}, got {length(value)}.")
  if (anyNA(value))
    cli_abort("{what} names can't be NA.")
  
  names(value) <- NULL
  trimmed      <- trimws(value)

  if (!all(nzchar(trimmed)))
    cli_abort("{what} names can't be empty.")
  if (length(dups <- unique(trimmed[duplicated(trimmed)])) > 0)
    cli_abort("Duplicated {what} name{?s}: {.val {dups}}.")

  if (length(i <- which(trimmed != value)) > 0)
    cli_warn("Trimming whitespace from {what} name{?s}: {.val {value[i]}}.")

  return (trimmed)
}
