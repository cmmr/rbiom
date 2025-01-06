
#' Calculate the alpha diversity of each sample.
#' 
#' @inherit documentation_default
#' 
#' @family alpha_diversity
#'        
#' @return A data frame of alpha diversity values.
#'         Each combination of sample/depth/\code{adiv} has its own row.
#'         Column names are \bold{.sample}, \bold{.depth}, \bold{.adiv}, 
#'         and \bold{.diversity}, followed by any metadata fields requested by 
#'         \code{md}.
#' 
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     # Subset to 10 samples.
#'     biom <- slice(hmp50, 1:10)
#'     adiv_table(biom)
#'     
#'     biom <- rarefy(biom)
#'     adiv_table(biom, adiv = ".all", md = NULL)

adiv_table <- function (
    biom, adiv = "Shannon", md = ".all", transform = "none" ) {
  
  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  cmd    <- sprintf("adiv_table(%s)", as.args(params, fun = adiv_table))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('adiv_table', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Check for valid arguments.
  #________________________________________________________
  validate_adiv(max = Inf)
  validate_biom_field('md', null_ok = TRUE, max = Inf)
  
  
  
  
  #________________________________________________________
  # Compute adiv values
  #________________________________________________________
  mtx <- adiv_matrix(biom = biom, transform = transform)
  mtx <- mtx[, unique(c('Depth', adiv))]
  tbl <- tibble(
    '.sample'    = rownames(mtx)[row(mtx)] %>% factor(levels = rownames(mtx)),
    '.depth'     = unname(mtx[,'Depth'][row(mtx)]),
    '.adiv'      = colnames(mtx)[col(mtx)] %>% factor(levels = adiv),
    '.diversity' = as.numeric(mtx) ) %>%
    dplyr::filter(!is.na(.data$.adiv))
  
  
  
  #________________________________________________________
  # Add Metadata
  #________________________________________________________
  if (!is.null(md))
    tbl %<>% left_join(
      y  = biom$metadata[,unique(c('.sample', md))], 
      by = '.sample' )
  
  
  
  #________________________________________________________
  # Descriptive label for y-axis.
  #________________________________________________________
  if (length(adiv) == 1) {
    
    resp_label <- switch(
      EXPR    = adiv,
      'OTUs'  = "Observed OTUs",
      'Depth' = "Sequencing Depth",
      paste(adiv, "Diversity") )
    
  } else {
    
    resp_label <- "Alpha Diversity"
  }
  
  
  
  
  tbl %<>% as_rbiom_tbl()
  attr(tbl, 'cmd')        <- cmd
  attr(tbl, 'response')   <- ".diversity"
  attr(tbl, 'resp_label') <- resp_label
  
  
  
  set_cache_value(cache_file, tbl)
  
  
  return (tbl)
}



#' Create a matrix of samples x alpha diversity metrics.
#' 
#' @inherit documentation_default
#'        
#' @return A numeric matrix with samples as rows and columns named 
#'         \bold{Depth}, \bold{OTUs}, \bold{Shannon}, \bold{Chao1}, 
#'         \bold{Simpson}, and \bold{InvSimpson}.
#'     
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     biom <- slice_head(hmp50, n = 5)
#'     
#'     adiv_matrix(biom)

adiv_matrix <- function (biom, transform = "none") {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment())
  cmd    <- sprintf("adiv_matrix(%s)", as.args(params, fun = adiv_matrix))
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('adiv_matrix', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  remove("params")
  
  
  #________________________________________________________
  # Check for valid arguments.
  #________________________________________________________
  validate_var_choices('transform', c("none", "rank", "log", "log1p", "sqrt"))
  
  
  
  #________________________________________________________
  # We want a numeric matrix of samples x adiv metrics
  #________________________________________________________
  mtx <- biom$counts
  df  <- rcpp_alpha_div(mtx)
  
  
  
  #________________________________________________________
  # Optionally transform the computed diversity values.
  #________________________________________________________
  if (transform != "none") # "rank", "log", "log1p", "sqrt", "percent"
    for (i in c('OTUs', 'Shannon', 'Chao1', 'Simpson', 'InvSimpson')) {
      if (transform == "percent") {
        df[[i]] <- tryCatch(df[[i]] / sum(df[[i]]), error = function (e) { warning(e); df[[i]] })
      } else {
        df[[i]] <- do.call(`::`, list('base', transform))(df[[i]])
      }
    }
  
  
  
  #________________________________________________________
  # Move sample names from df$.sample to rownames(df).
  #________________________________________________________
  rownames(df) <- df[['.sample']]
  mtx <- as.matrix(df[,-1,drop=FALSE])
  
  
  
  attr(mtx, 'cmd') <- cmd
  set_cache_value(cache_file, mtx)
  
  return (mtx)
}



#' Summarize the taxa observations in each sample.
#' 
#' @inherit documentation_default
#' 
#' @family samples
#' @family rarefaction
#' @family taxa_abundance
#' 
#' @param sort  Sort the result. Options: `NULL` - don't sort; `"asc"` - in 
#'        ascending order (smallest to largest); `"desc"` - in descending order 
#'        (largest to smallest). Ignored when the result is not a simple 
#'        numeric vector. Default: `NULL`
#' 
#' @param FUN  The function to apply to each column of `taxa_matrix()`.
#' 
#' @param ...  Optional arguments to `FUN`.
#' 
#' @return For `sample_sums`, A named numeric vector of the number of 
#'         observations in each sample. For `sample_apply`, a named vector or 
#'         list with the results of `FUN`. The names are the taxa IDs.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     library(ggplot2)
#'     
#'     sample_sums(hmp50, sort = 'asc') %>% head()
#'     
#'     # Unique OTUs and "cultured" classes per sample
#'     nnz <- function (x) sum(x > 0) # number of non-zeroes
#'     sample_apply(hmp50, nnz, 'otu') %>% head()
#'     sample_apply(hmp50, nnz, 'class', unc = 'drop') %>% head()
#'     
#'     # Number of reads in each sample's most abundant family
#'     sample_apply(hmp50, base::max, 'f', sort = 'desc') %>% head()
#'     
#'     ggplot() + geom_histogram(aes(x=sample_sums(hmp50)), bins = 20)

sample_sums <- function (biom, rank = -1, sort = NULL, unc = "singly") {
  sample_apply(biom, FUN = base::sum, rank, sort, unc)
}



#' @rdname sample_sums
#' @export

sample_apply <- function (biom, FUN, rank = -1, sort = NULL, unc = "singly", ...) {
  
  stopifnot(is.function(FUN))
  validate_var_choices('sort', choices = c('asc', 'desc'), null_ok = TRUE)
  
  mtx <- taxa_matrix(biom = biom, rank = rank, sparse = TRUE, unc = unc)
  
  res <- if (identical(FUN, base::sum))  { slam::col_sums(mtx)
  } else if (identical(FUN, base::mean)) { slam::col_means(mtx)
  } else { slam::colapply_simple_triplet_matrix(mtx, FUN, ...) }
  
  if (is.numeric(res) && !is.null(sort))
    res <- base::sort(res, decreasing = (sort == 'desc'))
  
  return (res)
}


