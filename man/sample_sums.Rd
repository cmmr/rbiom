% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adiv_table.r
\name{sample_sums}
\alias{sample_sums}
\title{Sum the observations in each sample.}
\usage{
sample_sums(biom)
}
\arguments{
\item{biom}{An \link[=rbiom_objects]{rbiom object}, such as from \code{\link[=as_rbiom]{as_rbiom()}}.
Any value accepted by \code{\link[=as_rbiom]{as_rbiom()}} can also be given here.}
}
\value{
A named numeric vector of the number of observations in each
sample. The names are the sample IDs.
}
\description{
Sum the observations in each sample.
}
\examples{
    library(rbiom)
    
    sample_sums(hmp50) \%>\% sort() \%>\% head()
    
    hist(sample_sums(hmp50))
}
\seealso{
Other samples: 
\code{\link{pull.rbiom}()}

Other rarefaction: 
\code{\link{rare_corrplot}()},
\code{\link{rare_multiplot}()},
\code{\link{rare_stacked}()},
\code{\link{rarefy}()},
\code{\link{rarefy_cols}()}
}
\concept{rarefaction}
\concept{samples}
