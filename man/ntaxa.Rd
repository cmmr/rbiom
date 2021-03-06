% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessors.r
\name{ntaxa}
\alias{ntaxa}
\title{Number of taxa in a BIOM.}
\usage{
ntaxa(biom)
}
\arguments{
\item{biom}{A \code{BIOM} object, as returned from \link{read.biom}.}
}
\value{
The number of taxa present.
}
\description{
Number of taxa in a BIOM.
}
\examples{
    library(rbiom)
    
    infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
    biom <- read.biom(infile)
    
    ntaxa(biom)

}
\seealso{
Other accessor functions: 
\code{\link{counts}()},
\code{\link{info}()},
\code{\link{metadata}()},
\code{\link{nsamples}()},
\code{\link{phylogeny}()},
\code{\link{sample.names}()},
\code{\link{sequences}()},
\code{\link{taxa.names}()},
\code{\link{taxa.ranks}()},
\code{\link{taxonomy}()}
}
\concept{accessor functions}
