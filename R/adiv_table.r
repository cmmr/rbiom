
#' Calculate the alpha diversity of each sample.
#' 
#' @inherit documentation_default
#' 
#' @family alpha_diversity
#'        
#' @return A data frame of alpha diversity values. \cr
#'         Each combination of sample/depth/\code{adiv} has its own row. \cr
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
    biom, adiv = "Shannon", md = ".all", trans = "none" ) {
  
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
  mtx <- adiv_matrix(biom = biom, trans = trans)
  mtx <- mtx[, unique(c('Depth', adiv))]
  tbl <- tibble(
    '.sample'    = rownames(mtx)[row(mtx)] %>% factor(levels = rownames(mtx)),
    '.depth'     = mtx[,'Depth'][row(mtx)],
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
    
    resp_label <- "\u03B1 Diversity" %>% 
      aa(display = '"\\u03B1 Diversity"')
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

adiv_matrix <- function (biom, trans = "none") {
  
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
  validate_var_choices('trans', c("none", "rank", "log", "log1p", "sqrt"))
  
  
  
  #________________________________________________________
  # We want a numeric matrix of samples x adiv metrics
  #________________________________________________________
  mtx <- biom$counts
  df  <- rcpp_alpha_div(mtx)
  
  
  
  #________________________________________________________
  # Optionally transform the computed diversity values.
  #________________________________________________________
  if (trans != "none") # "rank", "log", "log1p", "sqrt"
    for (i in c('OTUs', 'Shannon', 'Chao1', 'Simpson', 'InvSimpson'))
      df[[i]] <- do.call(`::`, list('base', trans))(df[[i]])
  
  
  
  #________________________________________________________
  # Move sample names from df$.sample to rownames(df).
  #________________________________________________________
  rownames(df) <- df[['.sample']]
  mtx <- as.matrix(df[,-1,drop=FALSE])
  
  
  
  attr(mtx, 'cmd') <- cmd
  set_cache_value(cache_file, mtx)
  
  return (mtx)
}



#' Sum the observations in each sample.
#' 
#' @inherit documentation_default
#' 
#' @family samples
#' @family rarefaction
#' 
#' @return A named numeric vector of the number of observations in each 
#'         sample. The names are the sample IDs.
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     library(ggplot2)
#'     
#'     sample_sums(hmp50) %>% sort() %>% head()
#'     
#'     ggplot() + geom_histogram(aes(x=sample_sums(hmp50)), bins = 20)

sample_sums <- function (biom) {
  biom <- as_rbiom(biom)
  col_sums(biom$counts)
}


