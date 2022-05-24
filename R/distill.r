#' Convert to a data.frame with metadata and diversity/abundance values.
#' 
#' @name distill
#' 
#' @param biom   A BIOM object, as returned from \link{read.biom}.
#' 
#' @param metric   The diversity/abundance values of interest. Options are:
#'        \describe{
#'            \item{Alpha Diversity Metrics}{
#'              \bold{OTUs}, \bold{Shannon}, \bold{Chao1}, \bold{Simpson}, 
#'              and/or \bold{InvSimpson}.
#'            }
#'            \item{Beta Diversity Metrics}{
#'              \bold{manhattan}, \bold{euclidean}, \bold{bray-curtis}, 
#'              \bold{jaccard}, or \bold{unifrac}. Use in combination with the 
#'              \code{weighted} parameter. Metadata column names can be 
#'              prefixed with \bold{==} or \bold{!=} to limit distance
#'              calculations to \emph{within} or \emph{between}, respectively, 
#'              those categories. See examples below.
#'            }
#'            \item{Taxa Abundances}{
#'              \bold{Kingdom}, \bold{Phylum}, \bold{Class}, \bold{Order}, 
#'              \bold{Family}, \bold{Genus}, \bold{Species}, \bold{Strain}, or 
#'              \bold{OTU}. Supported ranks will vary by biom. Run
#'              \code{taxa.ranks(biom)} to see the available options.
#'            }
#'           }
#'           
#' @param rarefy   Should the dataset be rarefied first? When \bold{metric} is 
#'        an alpha diversity metric, this 'rarefy' parameter is passed on 
#'        directly to \code{alpha.div()}. (Default: FALSE)
#'        
#' @param weighted   When \bold{metric} is a beta diversity metric, should it 
#'        be run in in weighted mode? (Default: TRUE)
#'        
#' @param long  Pivot the returned data to long format?
#'        \describe{
#'          \item{\bold{FALSE}}{ Each metric has its own column. }
#'          \item{\bold{TRUE}}{  "Sample", "Metric" and "Diversity" are the 
#'                               columns returned. Rows are added to attain all 
#'                               combinations of samples x metrics. (Default) }
#'        }
#'        
#' @param md   Include metadata in the output data frame? Options are: 
#'        \describe{
#'          \item{\code{FALSE}}{ Don't include metadata. (Default) }
#'          \item{\code{TRUE}}{ Include all metadata. }
#'          \item{\emph{character vector}}{ Include only the specified metadata 
#'                                          columns. }
#'        } 
#'        
#' @param safe   Should autogenerated columns be prefixed with a "." to avoid 
#'        conflicting with metadata column names? (Default: FALSE)
#'
#' @return A \code{data.frame} object. The first column will be named 
#'         \bold{Sample}, and when possible, the rownames will hold the sample 
#'         name as well. For beta diversity metrics the first two columns will 
#'         be named \bold{Sample1} and \bold{Sample2}. The remaining columns 
#'         will be largely dependent on \code{metric}, \code{long}, and 
#'         \code{md}. See examples below.
#'         
#' @export
#' @seealso \code{\link{stats.table}}
#' @examples
#'     library(rbiom)
#'     
#'     distill(hmp50, "UniFrac", md=c("==Body Site", "!=Sex")) %>% head()
#'     
#'     distill(hmp50, c("Shannon", "OTUs"), md=c("BMI", "Sex")) %>% head()
#'     
#'     distill(hmp50, "Phylum", long=FALSE, md=FALSE)[1:4,1:4]
#'     
#'     distill(hmp50, "Phylum", long=FALSE, md=c("Age", "Body Site"))[1:4,1:6]
#'
distill <- function (biom, metric, weighted = TRUE, rarefy = FALSE, long = TRUE, md = TRUE, safe = FALSE) {
  
  if (!is(biom, 'BIOM'))
    stop ("Input for 'biom' must be a BIOM object.")
  
  
  metric <- validate_metrics(biom, metric, multi=TRUE)
  mode   <- head(attr(metric, 'mode', exact = TRUE), 1)
  
  if (isTRUE(rarefy) && !mode %in% c('adiv'))
    biom <- rbiom::rarefy(biom)
  
  
  if (mode == "adiv") {
    df <- alpha.div(biom, metrics = metric, rarefy = rarefy, long = long, md = md, safe = safe)
    
  } else if (mode == "rank") {
    
    # Return specific rank(s) or all ranks
    ranks <- as.vector(metric)
    if (identical(ranks, "Rank"))
      ranks <- taxa.ranks(biom)
    
    df      <- taxa.rollup(biom, rank = ranks[[1]], long = long, md = md, safe = safe)
    df_rank <- rep(ranks[[1]], nrow(df))
    taxa_in <- attr(df, 'taxa_in', exact = TRUE)
    
    for (rank in ranks[-1]) {
      x <- taxa.rollup(biom, rank = rank, long = long, md = md, safe = safe)
    
      if (identical(taxa_in, "cols")) {
        # Don't duplicate metadata columns
        df <- cbind(df, x[,sort(unique(taxonomy(biom, rank))),drop=F])
        
      } else {
        df      <- rbind(df, x)
        df_rank <- c(df_rank, rep(rank, nrow(x)))
      }
    }
    
    # Put a 'Rank' column in the second position
    if (!identical(taxa_in, "cols") && (isTRUE(safe) || length(ranks) > 1)) {
      df[[ifelse(safe, '.rank', 'Rank')]] <- factor(df_rank, levels=ranks)
      df <- df[,order(!colnames(df) %in% c(".sample", ".rank", "Sample", "Rank")),drop=FALSE]
    }
    
  } else if (mode == "taxon") {
    
    rank <- names(which.max(apply(taxonomy(biom), 2L, function (x) sum(x == metric))))
    df   <- taxa.rollup(biom, rank = rank, long = long, md = md, safe = safe)
    
    # Only return the single requested taxon
    if (long) {
      df <- df[df[[ifelse(safe, ".taxa", "Taxa")]] == metric,,drop=FALSE]
    } else {
      keep <- c('Sample', '.sample', colnames(metadata(biom)), metric)
      df   <- df[,colnames(df) %in% keep,drop=FALSE]
    }
    
  } else if (mode == "bdiv") {
    
    if (is.character(md)) md <- paste0(attr(md, 'op', exact = TRUE), md)
    df <- beta.div(biom, method = metric, weighted = weighted, long = long, md = md, safe = safe)
    
  } else {
    stop("Don't know how to distill metric '", metric, "'.")
  }
  
  
  #--------------------------------------------------------------
  # Flag the response variable
  #--------------------------------------------------------------
  if (isTRUE(long)) {
    if        (mode == "adiv")  { attr(df, 'response') <- "Diversity"
    } else if (mode == "bdiv")  { attr(df, 'response') <- "Distance"
    } else                      { attr(df, 'response') <- "Abundance" }
    
    if (isTRUE(safe))
      attr(df, 'response') <- ".value"
  }
  
  
  return (df)
}



