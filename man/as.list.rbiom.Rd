% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_methods.r
\name{as.list.rbiom}
\alias{as.list.rbiom}
\title{Convert an rbiom object to a base R list.}
\usage{
\method{as.list}{rbiom}(biom)
}
\arguments{
\item{biom}{An \link[=rbiom_objects]{rbiom object}, such as from \code{\link[=as_rbiom]{as_rbiom()}}.}
}
\value{
A list with names
\code{c('counts', 'metadata', 'taxonomy', 'tree', 'sequences', 'id', 'comment', 'date', 'generated_by')}.
}
\description{
Convert an rbiom object to a base R list.
}
\concept{conversion}
