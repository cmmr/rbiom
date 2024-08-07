% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_biom.r
\name{biom_file_format}
\alias{biom_file_format}
\title{Inspects a file to see if it's most likely hdf5, json, or tab-delimited.
Also reports gzip or bzip2 compression.}
\usage{
biom_file_format(file)
}
\arguments{
\item{file}{The path to an rbiom file.}
}
\value{
One of \code{c("tsv", "tsv.gz", "tsv.bz2", "json", "json.gz", 
        "json.bz2", "hdf5")}.
}
\description{
Inspects a file to see if it's most likely hdf5, json, or tab-delimited.
Also reports gzip or bzip2 compression.
}
\examples{
    fp <- system.file("extdata", "hmp50.bz2", package = "rbiom")
    biom_file_format(file = fp)

}
\keyword{internal}
