
#' Calculate the alpha diversity of each sample.
#' 
#' @inherit documentation_default
#' 
#' @family alpha_diversity
#' 
#' @param rarefy  Control how/whether rarefactions are done prior to alpha 
#'        diversity computations. Options are: 
#'        \itemize{
#'          \item{\code{FALSE} - }{
#'            Use each sample's current set of observations without applying
#'            any rarefaction.}
#'          \item{\code{TRUE} - }{
#'            Automatically select and apply a single rarefaction. }
#'          \item{\code{"multi"}, \code{"multi_log"}, \code{"multi_even"} - }{
#'            Automatically select and apply multiple rarefactions.
#'            \code{"multi"} provides \code{"multi_log"} at the low end 
#'            and \code{"multi_even"} at the high end. }
#'          \item{\emph{integer vector} - }{
#'            Rarefy at the specified depth(s). }
#'        }
#'        Default: \code{FALSE}
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
#'     biom <- sample_select(hmp50, 1:6)
#'     
#'     adiv_table(biom)
#'     
#'     adiv_table(sample_rarefy(biom), adiv = "all", md = FALSE)

adiv_table <- function (
    biom, rarefy = FALSE, adiv = "Shannon", md = TRUE ) {
  
  validate_biom()
  
  params  <- eval_envir(environment())
  history <- append_history('tbl ', params)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  
  #________________________________________________________
  # Check for valid arguments.
  #________________________________________________________
  validate_adiv(max = Inf)
  validate_meta('md', null_ok = TRUE, max = Inf)
  
  
  #________________________________________________________
  # Define the rarefaction depths to sample at
  #________________________________________________________
  rLvls  <- NA
  
  if (eq(rarefy, TRUE)) {
    rLvls <- rare_suggest(biom)
    
  } else if (is.numeric(rarefy)) {
    rLvls <- sort(rarefy)
    
  } else if (is.character(rarefy)) {
    upper <- fivenum(sample_sums(biom))[[4]]
    
    if (eq(rarefy, "multi")) {
      
      # Log intervals until rLvl/2, then even intervals until rLvl*2
      rLvl  <- rare_suggest(biom)
      rLvls <- 10 ** (c(1,3,5,7,8,9) * (log10(rLvl / 2) / 10))
      rLvls <- c(rLvls, seq(from = rLvl / 2, to = rLvl * 2, length.out = 5))
      rLvls <- floor(rLvls)
      remove("rLvl")
      
    } else if (eq(rarefy, "multi_log")) {
      rLvls <- floor(10 ** (c(1:10) * (log10(upper) / 10)))
      
    } else if (eq(rarefy, "multi_even")) {
      rLvls <- floor(seq(from = 5, to = upper, length.out = 10))
    }
    remove("upper")
  }
  
  
  
  #________________________________________________________
  # Compute adiv values at each rLvl.
  #________________________________________________________
  tbl    <- NULL
  counts <- otu_matrix(biom, sparse = TRUE)
  
  # If non-integer counts, scale each sample to 10k "reads".
  if (!all(counts$v %% 1 == 0))
    counts <- sample_rarefy(biom) %>% otu_matrix(sparse = TRUE)
  
  for (rLvl in rLvls) {
    
    otus <- if (is.na(rLvl)) counts else rcpp_rarefy(counts, rLvl, 0)
    df   <- rcpp_alpha_div(otus)
    mtx  <- as.matrix(df[,adiv,drop=FALSE])
    
    
    #________________________________________________________
    # Pivot Longer
    #________________________________________________________
    df <- tibble(
      '.sample'    = df[['.sample']][row(mtx)],
      '.depth'     = df[['Depth']][row(mtx)],
      '.adiv'      = colnames(mtx)[col(mtx)],
      '.diversity' = as.numeric(mtx) )
    
    
    tbl <- rbind(tbl, df)
  }
  
  tbl[['.adiv']] %<>% factor(levels = adiv)
  
  
  
  #________________________________________________________
  # Add Metadata
  #________________________________________________________
  if (is.null(md)) {
    tbl[['.sample']] %<>% as.factor()
    
  } else {
    tbl %<>% left_join(
      y  = sample_metadata(biom)[,unique(c('.sample', md))], 
      by = '.sample' )
  }
  
  
  #________________________________________________________
  # Command history tracking and other attributes.
  #________________________________________________________
  attr(tbl, 'response') <- ".diversity"
  
  
  attr(tbl, 'history') <- history
  set_cache_value(cache_file, tbl)
  
  return (tbl)
}



#' Create a matrix of samples x alpha diversity metrics.
#' 
#' @inherit documentation_default
#'     
#' @param rarefy  Control how/whether rarefactions are done prior to alpha 
#'        diversity computations. Options are: 
#'        \itemize{
#'          \item{\code{FALSE} - }{
#'            Use each sample's current set of observations without applying
#'            any rarefaction.}
#'          \item{\code{TRUE} - }{
#'            Automatically select and apply a single rarefaction. }
#'          \item{\emph{integer} - }{
#'            Rarefy to the specified depth. }
#'        }
#'        Default: \code{FALSE}
#'        
#' @return A numeric matrix with samples as rows and columns named 
#'         \bold{Depth}, \bold{OTUs}, \bold{Shannon}, \bold{Chao1}, 
#'         \bold{Simpson}, and \bold{InvSimpson}.
#'     
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     biom <- sample_select(hmp50, 1:6)
#'     adiv_matrix(biom)

adiv_matrix <- function (biom, rarefy=FALSE) {
  
  validate_biom()
  
  params  <- eval_envir(environment())
  history <- append_history('mtx ', params)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  remove("params")
  
  
  #________________________________________________________
  # Rarefy unless `rarefy`=FALSE.
  #________________________________________________________
  if (!isFALSE(rarefy)) {
    
    if (!(is_scalar_integerish(rarefy) || is_na(rarefy)))
      rarefy <- NULL
    
    biom <- sample_rarefy(biom, depth = rarefy)
  }
  
  
  #________________________________________________________
  # We want a numeric matrix of samples x adiv metrics
  #________________________________________________________
  df <- rcpp_alpha_div(biom[['counts']])
  rownames(df) <- df[['.sample']]
  mtx <- as.matrix(df[,-1,drop=FALSE])
  
  
  
  attr(mtx, 'history') <- history
  set_cache_value(cache_file, mtx)
  
  return (mtx)
}
