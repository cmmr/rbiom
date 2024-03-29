% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bdiv_stats.r
\name{bdiv_stats}
\alias{bdiv_stats}
\title{Test beta diversity vs categorical metadata.}
\usage{
bdiv_stats(
  biom,
  stat.by,
  bdiv = "Bray-Curtis",
  weighted = TRUE,
  split.by = NULL,
  within = NULL,
  between = NULL,
  tree = NULL,
  trans = "none",
  test = "means",
  p.adj = "fdr"
)
}
\arguments{
\item{biom}{An rbiom-class object, or data coercible with
\code{\link[=as_rbiom]{as_rbiom()}}.}

\item{stat.by}{The categorical or numeric metadata field over which
statistics should be calculated. Required.}

\item{bdiv}{Beta diversity distance algorithm(s) to use. Options are:
\code{"Bray-Curtis"}, \code{"Manhattan"}, \code{"Euclidean"},
\code{"Jaccard"}, and \code{"UniFrac"}. For \code{"UniFrac"}, a
phylogenetic tree must be present in \code{biom} or explicitly
provided via \verb{tree=}. Default: \code{"Bray-Curtis"} \cr\cr
Multiple/abbreviated values allowed.}

\item{weighted}{Take relative abundances into account. When
\code{weighted=FALSE}, only presence/absence is considered.
Default: \code{TRUE} \cr\cr
Multiple values allowed.}

\item{split.by}{Name(s) of metadata columns that the data should be split
by prior to any calculations. Default: \code{NULL}}

\item{within, between}{Metadata field(s) for intra- or inter- sample
comparisons. Default: \code{NULL}}

\item{tree}{A \code{phylo} object representing the phylogenetic
relationships of the taxa in \code{biom}. Only required when
computing UniFrac distances. Default: \code{biom$tree}}

\item{trans}{Transformation to apply. Options are:
\code{c("none", "rank", "log", "log1p", "sqrt")}. \code{"rank"} is useful for
correcting for non-normally distributions before applying regression
statistics. Default: \code{"none"}}

\item{test}{Permutational test for accessing significance. Options are:
\itemize{
\item{\code{"adonis2"} - }{ Permutational MANOVA; \code{\link[vegan:adonis]{vegan::adonis2()}}. }
\item{\code{"mrpp"} - }{ Multiple response permutation procedure; \code{\link[vegan:mrpp]{vegan::mrpp()}}. }
\item{\code{"none"} - }{ Don't run any statistics. }
}
Default: \code{"adonis2"} \cr\cr
Abbreviations are allowed.}

\item{p.adj}{Method to use for multiple comparisons adjustment of
p-values. Run \code{p.adjust.methods} for a list of available
options. Default: \code{"fdr"}}
}
\value{
A tibble data frame with summary statistics. \cr
The R code or generating these statistics is in \verb{$code}.
}
\description{
Test beta diversity vs categorical metadata.
}
\examples{
    library(rbiom)
    
    biom <- rarefy(hmp50)
    
    bdiv_stats(biom, stat.by = "Body Site", split.by = "Sex")
      
    bdiv_stats(biom, stat.by = "Sex", bdiv = c("bray", "unifrac"))
    
    # The R code used to compute the stats is in $code.
    tbl <- bdiv_stats(biom, stat.by = "Sex")
    tbl$code
}
\seealso{
Other beta_diversity: 
\code{\link{bdiv_boxplot}()},
\code{\link{bdiv_corrplot}()},
\code{\link{bdiv_heatmap}()},
\code{\link{bdiv_ord_plot}()},
\code{\link{bdiv_ord_table}()},
\code{\link{bdiv_table}()},
\code{\link{distmat_stats}()}

Other stats_tables: 
\code{\link{adiv_stats}()},
\code{\link{distmat_stats}()},
\code{\link{taxa_stats}()}
}
\concept{beta_diversity}
\concept{stats_tables}
