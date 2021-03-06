% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write.biom.r
\name{write.biom}
\alias{write.biom}
\title{Write counts, metadata, taxonomy, and phylogeny to a biom file.}
\usage{
write.biom(biom, file, format = "json")
}
\arguments{
\item{biom}{The BIOM object to save to the file.}

\item{file}{Path to the output file. If the file name ends in \code{.gz} 
or \code{.bz2}, the file contents will be compressed accordingly.}

\item{format}{Options are \bold{\dQuote{tab}}, 
\bold{\dQuote{json}}, and \bold{\dQuote{hdf5}}, 
corresponding to classic tabular format, biom format version 1.0 and 
biom version 2.1, respectively. Abbreviations are also accepted. See
\url{http://biom-format.org/documentation/} for details. NOTE: to write
HDF5 formatted BIOM files, the BioConductor R package \code{rhdf5} must
be installed.}
}
\value{
On success, returns \code{NULL} invisibly.
}
\description{
Write counts, metadata, taxonomy, and phylogeny to a biom file.
}
