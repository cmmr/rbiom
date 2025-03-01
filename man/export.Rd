% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_biom.r
\name{export}
\alias{export}
\alias{write_mothur}
\alias{write_qiime2}
\title{Export data to QIIME 2 or mothur.}
\usage{
write_mothur(biom, dir = tempfile(), prefix = "biom_")

write_qiime2(biom, dir = tempfile(), prefix = "biom_")
}
\arguments{
\item{biom}{An \link[=rbiom_objects]{rbiom object}, such as from \code{\link[=as_rbiom]{as_rbiom()}}.
Any value accepted by \code{\link[=as_rbiom]{as_rbiom()}} can also be given here.}

\item{dir}{Where to save the files. If the directory doesn't exist, it will
be created. Default: \code{tempfile()}}

\item{prefix}{A string to prepend to each file name. Default: \code{'biom_'}}
}
\value{
The normalized directory path that was written to (invisibly).
}
\description{
Populates a directory with the following files, formatted according to
QIIME 2 or mothur's specifications.
\itemize{
\item \code{biom_counts.tsv}
\item \code{biom_metadata.tsv}
\item \code{biom_taxonomy.tsv}
\item \code{biom_tree.nwk}
\item \code{biom_seqs.fna}
}

\code{biom_counts.tsv} will always be created. The others are dependent on
whether the content is present in the \code{biom} argument.
}
\examples{
    library(rbiom)
    
    tdir <- tempfile()
    
    write_qiime2(hmp50, tdir, 'qiime2_')
    write_mothur(hmp50, tdir, 'mothur_')
    
    list.files(tdir)
    
    readLines(file.path(tdir, 'qiime2_metadata.tsv'), n = 4)
    
    readLines(file.path(tdir, 'mothur_taxonomy.tsv'), n = 3)
    
    unlink(tdir, recursive = TRUE)

}
