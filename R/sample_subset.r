#' Subset samples using the BIOM's metadata or taxonomy.
#' 
#' @section Taxonomic abundance filtering:
#' 
#' For taxonomic subsetting, several functions are added or overridden to
#' behave as expected within the subsetting expression. They are:
#' 
#' \code{mean()}, \code{median()}, \code{min()}, \code{max()}, \code{n()}, 
#' \code{count()}, \code{percent()}, \code{rank()}, and \code{apply()}.
#' 
#' Therefore you can write 
#' \code{sample_subset(hmp50, mean(Genus) >= 0.1)} and the returned BIOM object will 
#' contain only the genera that average at least 10% relative abundance across
#' all the samples.
#' 
#' If you want only orders that are present in three or more samples, you can 
#' do: \code{sample_subset(hmp50, count(Order) >= 3)}. To require presence in 25% of 
#' samples, you'd use: \code{sample_subset(hmp50, percent(Order) >= 0.25)}.
#' 
#' Both \code{count()} and \code{percent()} have default arguments of 
#' \code{gt=0, le=1, ge=NULL, lt=NULL}, which can be overridden to find, e.g., 
#' which genera comprise at least 2% of the community in 10% or more of the 
#' samples: \code{sample_subset(hmp50, percent(Genus, ge=0.02) >= 0.10)}. 
#' 
#' \emph{gt = greater than, ge = greater than or equal to. lt/le similarly with 
#' 'less than'.}
#' 
#' To keep only the top 5 most abundant genera (based on mean), run: 
#' \code{sample_subset(hmp50, rank(Genus) <= 5)}.
#' 
#' \code{apply()} allows you to run any function on a per-taxon basis. For 
#' example, to filter genera by the root mean square of their relative abundances:
#' \code{rms <- function(x) sqrt(mean(x^2)); sample_subset(hmp50, apply(Genus, rms) >= 0.1)}.
#' 
#' If you prefer to work on the raw values (e.g. read counts) instead of 
#' relative abundances, set \code{apply(..., raw = TRUE)}. For instance: 
#' \code{sample_subset(hmp50, apply(Genus, mean, raw=TRUE) >= 100)}.
#' 
#' 
#' @param biom,x   A BIOM object, as returned from [read_biom()].
#' 
#' @param expr   Logical expression to run on the metadata or taxonomy (not both)
#'        to identify samples or taxa to retain.
#'        
#' @param env   The environment to search for variables used in \code{expr}.
#'        Default: \code{parent.frame()}
#'        
#' @param drop.na   When \code{expr} is e.g. \code{Age > 30}, should 
#'        \code{!is.na(Age)} be automatically applied too?
#'        Default: \code{TRUE}
#'        
#' @param refactor  When \code{expr} is e.g. 
#'        \code{`Body Site` \%in\% c("Stool", "Saliva")}, should
#'        \code{`Body Site`} be redefined as 
#'        \code{factor(`Body Site`, levels=c('Stool', 'Saliva'))}?
#'        Applies to categorical metadata only.
#'        Default: \code{TRUE}
#'
#' @param fast  Should subsetting the phylogenetic tree and sequences be 
#'        skipped? These slow steps are often not necessary.
#'        Default: \code{FALSE}
#'        
#' @return A \code{BIOM} object.
#' 
#' @seealso [sample_select()]
#' 
#' @export
#' @examples
#'   \dontrun{
#'     library(rbiom) 
#'     
#'     sample_subset(hmp50, `Body Site` %in% c("Saliva", "Stool"))
#'     sample_subset(hmp50, Age < 25 & BMI > 22)
#'     sample_subset(hmp50, Phylum %in% c("Firmicutes", "Actinobacteria"))
#'     sample_subset(hmp50, mean(Genus) > 0.1)
#'     sample_subset(hmp50, rank(Genus) <= 5)
#'     sample_subset(hmp50, a == b, list(a = as.name("Body Site"), b ="Saliva"))
#'  }
#'
sample_subset <- function (biom, expr, env = parent.frame(), drop.na = TRUE, refactor = TRUE, fast = FALSE) {
  
  stopifnot(is(biom, 'BIOM'))
  
  
  #________________________________________________________
  # Replace variable names with their values
  #________________________________________________________
  expr <- do.call(substitute, list(substitute(expr), as.list(env)))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  
  cache_file <- get_cache_file("sample_subset", params = list(
    biom     = eval(biom),
    expr     = rlang::enexpr(expr),
    drop.na  = eval(drop.na),
    refactor = eval(refactor),
    fast     = eval(refactor)
  ))
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Don't record sub-calls in the biom's history
  #________________________________________________________
  hist <- attr(biom, 'history', exact = TRUE)
  
  
  # #________________________________________________________
  # # Convert an expression to call. expression(Age > 9) => Age > 9
  # #________________________________________________________
  # 
  # if (mode(substitute(expr)) == "name")
  #   if (identical(expr[[1]], as.name("expression")))
  #     return (do.call(sample_subset, list(biom, expr[[1]], env, drop.na, refactor, fast)))
  # 
  # if (identical(substitute(expr)[[1]], as.name("expression")))
  #   return (do.call(sample_subset, list(biom, expr[[1]], env, drop.na, refactor, fast)))
  #
  # Just use do.call instead to make sure you pass a call in. E.g.:
  #   expr <- expression(Age < 30)
  #   biom <- do.call(subset, list(biom, expr[[1]]))
  #
  
  
  #________________________________________________________
  # Metadata or Taxonomy based subsetting?
  #________________________________________________________
  vars   <- all.vars(expr)
  ranks  <- taxa_ranks(biom)
  mdcols <- colnames(sample_metadata(biom))
  if (any(vars %in% ranks) && any(vars %in% mdcols))
    stop("subset expression must be either all metadata or all taxonomy.")
  mode  <- ifelse(any(vars %in% mdcols), "metadata", "taxonomy")
  remove("vars", "ranks", "mdcols")
  
  
  #________________________________________________________
  # Taxonomy vs metadata subsetting are done differently
  #________________________________________________________
  if (mode == "taxonomy") {
    envir <- c(
      as.list(as.data.frame(otu_taxonomy(biom))),
      list(OTU = otu_names(biom)) )
    
    
    #________________________________________________________
    # These functions get mapped to new apply; remember originals
    #________________________________________________________
    funcs <- c(
      'mean'  = mean, 'median' = median, 'n'     = length,
      'min'   = min,  'max'    = max,    'apply' = apply,
      'count' = function (x, gt=0, le=1, ge=NULL, lt=NULL) {
        lower_pass <- if (is_null(ge)) x >  gt else x >= ge
        upper_pass <- if (is_null(lt)) x <= le else x >  lt
        sum(lower_pass & upper_pass) },
      'percent' = function (x, gt=0, le=1, ge=NULL, lt=NULL) {
        lower_pass <- if (is_null(ge)) x >  gt else x >= ge
        upper_pass <- if (is_null(lt)) x <= le else x >  lt
        sum(lower_pass & upper_pass) / length(x) })
    
    #________________________________________________________
    # Take any function and run it on each taxon's rel. abundances.
    # `f` should take a vector of numbers and return one value.
    #________________________________________________________
    envir[['apply']] <- function (x, f = NULL, raw = FALSE, ...) {
      for (i in names(funcs)) assign(i, funcs[[i]])
      
      if (missing(f)) f <- as.character(match.call()[[1]])
      fn  <- ifelse(is.character(f), f, capture.output(substitute(f)))
      fun <- funcs[[fn]] %||% ifelse(is.function(f), f, get(fn))
      if (!is.function(fun)) stop("unknown function: ", fn)
      
      rank <- capture.output(substitute(x))
      if (!rank %in% taxa_ranks(biom))
        return (fun(x, ...))
        
      mtx  <- taxa_matrix(biom, rank)
      mtx  <- if (isTRUE(raw)) mtx else mtx / rowSums(mtx)
      res  <- apply(mtx, 2L, fun, ...)
      res[otu_taxonomy(biom, rank)[,1]]
    }
    for (i in names(funcs)) envir[[i]] <- envir[['apply']]
    
    
    #________________________________________________________
    # Rank calls apply then does additional post-processing.
    #________________________________________________________
    funcs %<>% c('rank' = rank)
    envir[['rank']] <- function (...) {
      for (i in names(funcs)) assign(i, funcs[[i]])
      means    <- envir[['apply']](f = mean, ...)
      rankings <- rank(-means, ties.method = "min")
      gapfill  <- sort(unique(rankings))
      gapfill  <- setNames(rank(gapfill), as.character(gapfill))
      gapfill[as.character(rankings)]
    }
    
    
    #________________________________________________________
    # Evaluate expr to determine which taxa to keep
    #________________________________________________________
    keep <- try(eval(expr, envir = envir), silent=TRUE)
    if (is(keep, "try-error") || is(keep, 'error'))
      stop(sprintf("Subset failed: %s", keep))
    
    
    biom$taxonomy <- biom$taxonomy[keep,,drop=F]
    biom <- biom_repair(biom, fast=fast)
    
  } else {
    envir <- sample_metadata(biom)
    
    #________________________________________________________
    # Evaluate expr to determine which samples to keep
    #________________________________________________________
    keep <- try(eval(expr, envir = envir), silent=TRUE)
    if (is(keep, "try-error") || is(keep, 'error'))
      stop(sprintf("Subset failed: %s", keep))
    
    #________________________________________________________
    # Keep only the indicated samples
    #________________________________________________________
    biom <- sample_select(biom, keep, fast=fast)
    
    
    #________________________________________________________
    # Remove NAs in subsetted columns
    #________________________________________________________
    if (isTRUE(drop.na))
      for (i in all.vars(expr))
        if (i %in% colnames(biom$metadata))
          biom <- sample_select(biom, !is.na(biom$metadata[[i]]))
    
    
    
    #________________________________________________________
    # Recursively traverse subsetting expression to find all
    # instances of, e.g., `Body Site` %in% c('Stool', 'Saliva')
    #________________________________________________________
    in_lists <- function (expr) {
      result <- list()
      
      if (isTRUE(expr[[1]] == "%in%")) {
        
        key <- as.character(expr[[2]])
        val <- expr[[3]]
        
        if (is(val, "call") && val[[1]] == "c")
          val <- sapply(val[-1], unlist)
        
        result[[key]] <- as.character(val)
        
      } else {
        
        for (i in seq_len(length(expr))) {
          if (length(expr[[i]]) > 1) {
            result <- c(result, in_lists(expr[[i]]))
          }
        }
      }
      
      result <- result[!duplicated(names(result))]
      
      return (result)
    }
    
    
    #________________________________________________________
    # Convert/update factors according to subset specs
    #________________________________________________________
    if (isTRUE(refactor)) {
      specs <- in_lists(expr)
      for (key in names(specs)) {
        val <- specs[[key]]
        col <- biom$metadata[[key]]
        if (all(unique(col) %in% val))
          biom$metadata[[key]] <- factor(as.character(col), levels = val)
      }
    }
    
  }
  
  
  #________________________________________________________
  # Attach sample_subset() call to provenance tracking
  #________________________________________________________
  
  cl <- match.call()
  cl <- cl[names(cl) != "env"]
  cl[[1]] <- as.name("sample_subset")
  cl[[2]] <- as.name("biom")
  cl[[3]] <- expr
  for (i in seq_along(names(cl))[-(1:3)]) {
    cl[i] <- list(eval.parent(cl[[i]]))
  }
  names(cl)[[2]] <- ""
  names(cl)[[3]] <- ""
  
  attr(biom, 'history') <- paste0(collapse = "\n", c(
    hist,
    sprintf("biom <- %s", deparse1(cl)) ))
  
  
  set_cache_value(cache_file, biom)
  return (biom)
}


