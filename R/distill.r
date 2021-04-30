#' Convert to a data.frame with metadata and diversity/abundance values.
#' 
#' @name distill
#' @param biom      A BIOM object, as returned from \link{read.biom}.
#' @param metric    The diversity/abundance values of interest. Options are:
#'     \describe{
#'         \item{Alpha Diversity Metrics (one or more)}{
#'           \bold{OTUs}, \bold{Shannon}, \bold{Chao1}, \bold{Simpson}, and/or \bold{InvSimpson}.
#'           You may also set \code{metric = "Diversity"} to get all five metrics.
#'         }
#'         \item{Beta Diversity Metrics (one only)}{
#'           \bold{manhattan}, \bold{euclidean}, \bold{bray-curtis}, \bold{jaccard}, 
#'           or \bold{unifrac}. Use in combination with the \code{weighted} parameter.
#'           Metadata column names can be prefixed with \bold{==} or \bold{!=} to limit distance
#'           calculations to \emph{within} or \emph{between}, respectively, those categories. See
#'           examples below. Setting \code{metric = "Distance"} will use \bold{unifrac} if a 
#'           phylogenetic tree is present, or \bold{bray-curtis} otherwise.
#'         }
#'         \item{Taxa Abundances (one only)}{
#'           \bold{Kingdom}, \bold{Phylum}, \bold{Class}, \bold{Order}, \bold{Family}, \bold{Genus}, 
#'           \bold{Species}, \bold{Strain}, or \bold{OTU}. Supported ranks will vary by biom. Run
#'           \code{taxa.ranks(biom)} to see the available options. Specifying 
#'           \code{metric = "Abundance"} will default to the most precise rank possible.
#'         }
#'        }
#' @param weighted   When \bold{metric} is a beta diversity metric, should it be run in in weighted
#'                   mode? (Default: TRUE)
#' @param long  Pivot the returned data to long format?
#'     \describe{
#'       \item{\bold{FALSE}}{ Each metric has its own column. (Default) }
#'       \item{\bold{TRUE}}{ "Sample", "Metric" and "Diversity" are the columns 
#'                       returned. Rows are added to attain all combinations of 
#'                       samples x metrics. }
#'     }
#' @param md  Include metadata in the output data frame? Options are: 
#'     \describe{
#'       \item{\code{FALSE}}{ Don't include metadata. (Default) }
#'       \item{\code{TRUE}}{ Include all metadata. }
#'       \item{\emph{character vector}}{ Include only the specified metadata columns. }
#'     } 
#' @return A \code{data.frame} object. The first column will be named \bold{Sample}, 
#'         and when possible, the rownames will hold the sample name as well. For
#'         beta diversity metrics the first two columns will be named \bold{Sample1}
#'         and \bold{Sample2}. The remaining columns will be largely dependent on
#'         \code{metric}, \code{long}, and \code{md}. See examples below.
#' @export
#' @seealso \code{\link{stats.table}}
#' @examples
#'     library(rbiom) 
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     distill(biom, "UniFrac", md=c("==Body Site", "!=Sex")) %>% head()
#'     
#'     distill(biom, c("Shannon", "OTUs"), md=c("BMI", "Sex")) %>% head()
#'     
#'     distill(biom, "Phylum", long=FALSE, md=FALSE)[1:4,1:4]
#'     
#'     distill(biom, "Phylum", long=FALSE, md=c("Age", "Body Site"))[1:4,1:6]
#'
distill <- function (biom, metric, weighted = TRUE, long = TRUE, md = TRUE) {
  
  if (!is(biom, 'BIOM'))
    stop ("Input for 'biom' must be a BIOM object.")
  
  vals <- validate_metrics(biom, metric)
  mode <- attr(vals, 'mode')
  
  if (mode == "adiv") {
    df <- alpha.div(biom, metrics = vals, long = long, md = md)
    attr(df, 'response') <- "Diversity"
    
  } else if (mode == "taxa") {
    df <- taxa.rollup(biom, rank = vals, long = long, md = md)
    
    if (isTRUE(long)) {
      names(df)[2] <- vals
      if (!is.null(attr(df, 'facet')))
        attr(df, 'facet') <- vals
    }
    if (isFALSE(long) && isFALSE(md)) {
      df <- data.frame(
        check.names = FALSE, 
        Sample = colnames(df),
        t(df))
    }
    
    attr(df, 'response') <- "Abundance"
    
  } else if (mode == "bdiv") {
    df <- beta.div(biom, method = vals, weighted = weighted, long = long, md = md)
    
    if (isFALSE(long) && isFALSE(md))
      df <- df %>% as.matrix() %>% as.data.frame()
    
    attr(df, 'response') <- "Distance"
  }
  
  
  return (df)
}



