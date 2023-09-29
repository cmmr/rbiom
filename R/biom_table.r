#' Convert to a data.frame with metadata and diversity/abundance values.
#' 
#' @noRd
#' 
#' @param biom   A BIOM object, as returned from [read_biom()].
#' 
#' @param metric   The diversity/abundance values of interest. Options are:
#'        \itemize{
#'          \item{\bold{Alpha Diversity Metrics - }}{
#'            \code{"OTUs"}, \code{"Shannon"}, \code{"Chao1"}, 
#'            \code{"Simpson"}, and/or \code{"InvSimpson"}.
#'          }
#'          \item{\bold{Beta Diversity Metrics - }}{
#'            \code{"Manhattan"}, \code{"Euclidean"}, \code{"Bray-Curtis"}, 
#'            \code{"Jaccard"}, or \code{"UniFrac"}. \cr
#'            Use in combination with the \code{weighted} parameter. \cr
#'            Metadata column names can be prefixed with \bold{==} or \bold{!=} 
#'            to limit distance calculations to \emph{within} or 
#'            \emph{between}, respectively, those categories. \cr
#'            See examples below.
#'          }
#'          \item{\bold{Taxa Abundances - }}{
#'            \code{"Kingdom"}, \code{"Phylum"}, \code{"Class"}, 
#'            \code{"Order"}, \code{"Family"}, \code{"Genus"}, 
#'            \code{"Species"}, \code{"Strain"}, or \code{"OTU"}. \cr
#'            Supported ranks will vary by biom. \cr
#'            Run \code{taxa_ranks(biom)} to see the available options.
#'          }
#'        }
#'        \bold{\emph{Required}}
#'           
#' @param rarefy   Should the dataset be rarefied first? When \bold{metric} is 
#'        an alpha diversity metric, this 'rarefy' parameter is passed on 
#'        directly to \code{adiv_table()}. Default: \code{FALSE}
#'        
#' @param weighted   When \bold{metric} is a beta diversity metric, should it 
#'        be run in in weighted mode? Default: \code{TRUE}
#'        
#' @param long  Pivot the returned data to long format?
#'        \itemize{
#'          \item{\code{FALSE} - }{ Each metric has its own column. }
#'          \item{\code{TRUE} - }{ 
#'            "Sample", "Metric" and "Diversity" are the columns returned. Rows 
#'            are added to attain all combinations of samples x metrics. }
#'        }
#'        Default: \code{TRUE}
#'        
#' @param md   Include metadata in the output data frame? Options are: 
#'        \itemize{
#'          \item{\code{FALSE} - }{ Don't include metadata. }
#'          \item{\code{TRUE} - }{ Include all metadata. }
#'          \item{\emph{character vector} - }{
#'            Include only the specified metadata columns. }
#'        }
#'        Default: \code{FALSE}
#'
#' @return A \code{data.frame} object. The first column will be named 
#'         \bold{Sample}, and when possible, the rownames will hold the sample 
#'         name as well. For beta diversity metrics the first two columns will 
#'         be named \bold{Sample1} and \bold{Sample2}. The remaining columns 
#'         will be largely dependent on \code{metric}, \code{long}, and 
#'         \code{md}. See examples below.
#' 
#' @export
#' @seealso [biom_stats()]
#' @examples
#'     library(rbiom)
#'     
#'     biom_distill(hmp50, "UniFrac", md=c("==Body Site", "!=Sex")) %>% head()
#'     
#'     biom_distill(hmp50, c("Shannon", "OTUs"), md=c("BMI", "Sex")) %>% head()
#'     
#'     biom_distill(hmp50, "Phylum", long=FALSE, md=FALSE)[1:4,1:4]
#'     
#'     biom_distill(hmp50, "Phylum", long=FALSE, md=c("Age", "Body Site"))[1:4,1:6]
#'
biom_distill <- function (biom, metric, weighted = TRUE, rarefy = FALSE, long = TRUE, md = TRUE) {
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- lapply(as.list(environment()), eval)
  cache_file <- get_cache_file("biom_distill", params)
  if (!is.null(cache_file) && Sys.setFileTime(cache_file, Sys.time()))
    return (readRDS(cache_file))
  
  
  if (!is(biom, 'BIOM'))
    stop ("Input for 'biom' must be a BIOM object.")
  
  
  metric <- validate_metrics(biom, metric, multi=TRUE)
  mode   <- head(attr(metric, 'mode', exact = TRUE), 1)
  
  if (isTRUE(rarefy) && !mode %in% c('adiv'))
    biom <- sample_rarefy(biom)
  
  
  if (mode == "adiv") {
    
    df <- adiv_table(
      biom   = biom, 
      adiv   = metric, 
      rarefy = rarefy, 
      long   = long, 
      md     = md )
    
    
  } else if (mode == "rank") {
    
    # Return specific rank(s) or all ranks
    ranks <- as.vector(metric)
    if (identical(ranks, "Rank"))
      ranks <- taxa_ranks(biom)
    
    df      <- taxa_rollup(biom, rank = ranks[[1]], long = long, md = md)
    df_rank <- rep(ranks[[1]], nrow(df))
    taxa_in <- attr(df, 'taxa_in', exact = TRUE)
    attrs   <- attributes(df)
    
    for (rank in ranks[-1]) {
      x <- taxa_rollup(biom, rank = rank, long = long, md = md)
    
      if (identical(taxa_in, "cols")) {
        # Don't duplicate metadata columns
        df <- cbind(df, x[,sort(unique(otu_taxonomy(biom, rank)[,1])),drop=F])
        
      } else {
        df      <- rbind(df, x)
        df_rank <- c(df_rank, rep(rank, nrow(x)))
      }
    }
    
    # Put a 'Rank' column in the second position
    if (!identical(taxa_in, "cols")) {
      df[['.rank']] <- factor(df_rank, levels=ranks)
      df <- df[,order(!colnames(df) %in% c(".sample", ".rank", "Sample", "Rank")),drop=FALSE]
    }
    
    for (i in setdiff(names(attrs), names(attributes(df))))
      attr(df, i) <- attrs[[i]]
    
    
  } else if (mode == "taxon") {
    
    rank  <- names(which.max(apply(otu_taxonomy(biom), 2L, function (x) sum(x == metric))))
    df    <- taxa_rollup(biom, rank = rank, long = long, md = md)
    attrs <- attributes(df)
    
    # Only return the single requested taxon
    if (long) {
      df <- df[df[[".taxa"]] == metric,,drop=FALSE]
    } else {
      keep <- c('Sample', '.sample', colnames(sample_metadata(biom)), metric)
      df   <- df[,colnames(df) %in% keep,drop=FALSE]
    }
    
    for (i in setdiff(names(attrs), names(attributes(df))))
      attr(df, i) <- attrs[[i]]
    
    
  } else if (mode == "bdiv") {
    
    if (is.character(md)) md <- paste0(attr(md, 'op', exact = TRUE), md)
    df <- bdiv_table(biom, bdiv = metric, weighted = weighted, md = md)
    
  } else {
    stop("Don't know how to biom_distill metric '", metric, "'.")
  }
  
  attr(df, 'mode') <- mode
  
  
  
  #________________________________________________________
  # Update the response variable for consistency.
  #________________________________________________________
  if (isTRUE(long))
    df <- rename_response(df, ".value")
  
  
  set_cache_value(cache_file, df)
  return (df)
}



