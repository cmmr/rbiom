
#' Create a matrix of samples x alpha diversity metrics.
#' 
#' @inherit documentation_default
#'     
#' @param rarefy  Control how/whether rarefactions are done prior to alpha 
#'        diversity computations. Options are: 
#'        \itemize{
#'          \item{`FALSE` - }{
#'            Use each sample's current set of observations without applying
#'            any rarefaction.}
#'          \item{`TRUE` - }{
#'            Automatically select and apply a single rarefaction. }
#'          \item{\emph{integer} - }{
#'            Rarefy to the specified depth. }
#'        }
#'        Default: `FALSE`
#'        
#' @return A numeric matrix with samples as rows and columns named 
#'         \bold{Depth}, \bold{OTUs}, \bold{Shannon}, \bold{Chao1}, 
#'         \bold{Simpson}, and \bold{InvSimpson}.
#'     
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     biom <- hmp50$clone()
#'     biom$metadata %<>% head()
#'     
#'     adiv_matrix(biom)

adiv_matrix <- function (biom, rarefy=FALSE) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment())
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file()
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  remove("params")
  
  
  mtx <- biom$counts
  
  #________________________________________________________
  # Rarefy unless `rarefy`=FALSE.
  #________________________________________________________
  if (!isFALSE(rarefy)) {
    
    if (!is_scalar_integerish(rarefy) || is_na(rarefy))
      rarefy <- 'auto'
    
    mtx %<>% rarefy_cols(mtx, depth = rarefy)
  }
  
  
  #________________________________________________________
  # We want a numeric matrix of samples x adiv metrics
  #________________________________________________________
  df <- rcpp_alpha_div(mtx)
  rownames(df) <- df[['.sample']]
  mtx <- as.matrix(df[,-1,drop=FALSE])
  
  
  
  set_cache_value(cache_file, mtx)
  
  return (mtx)
}



#' Calculate the alpha diversity of each sample.
#' 
#' @inherit documentation_default
#' 
#' @family alpha_diversity
#' 
#' @param rarefy  Control how/whether rarefactions are done prior to alpha 
#'        diversity computations. Options are: 
#'        \itemize{
#'          \item{`FALSE` - }{
#'            Use each sample's current set of observations without applying
#'            any rarefaction.}
#'          \item{`TRUE` - }{
#'            Automatically select and apply a single rarefaction. }
#'          \item{\code{"multi"}, \code{"multi_log"}, \code{"multi_even"} - }{
#'            Automatically select and apply multiple rarefactions.
#'            \code{"multi"} provides \code{"multi_log"} at the low end 
#'            and \code{"multi_even"} at the high end. }
#'          \item{\emph{integer vector} - }{
#'            Rarefy at the specified depth(s). }
#'        }
#'        Default: `FALSE`
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
    biom, rarefy = FALSE, adiv = "Shannon", md = ".all" ) {
  
  biom <- as_rbiom(biom)
  
  params <- eval_envir(environment())
  
  
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
    rLvls <- rare_suggest(biom$counts)
    
  } else if (is.numeric(rarefy)) {
    rLvls <- sort(rarefy)
    
  } else if (is.character(rarefy)) {
    upper <- fivenum(sample_sums(biom))[[4]]
    
    if (eq(rarefy, "multi")) {
      
      # Log intervals until rLvl/2, then even intervals until rLvl*2
      rLvl  <- rare_suggest(biom$counts)
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
  counts <- biom$counts
  
  # If non-integer counts, scale each sample to 10k "reads".
  if (!all(counts$v %% 1 == 0))
    counts %<>% rarefy_cols()
  
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
      y  = biom$metadata[,unique(c('.sample', md))], 
      by = '.sample' )
  }
  
  
  attr(tbl, 'response') <- ".diversity"
  
  
  set_cache_value(cache_file, tbl)
  
  return (tbl)
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
#'     
#'     sample_sums(hmp50) %>% sort() %>% head()
#'     
#'     hist(sample_sums(hmp50))

sample_sums <- function (biom) {
  biom <- as_rbiom(biom)
  col_sums(biom$counts)
}


