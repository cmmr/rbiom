#' Write sequences from a BIOM object to a file in fasta format.
#'
#' @param biom     A BIOM object containing sequences.
#' @param outfile  Path to the output fasta file. Files ending in \code{.gz} or
#'                   \code{.bz2} will be compressed accordingly.
#' @return On success, returns \code{NULL} invisibly.
#' 
#' @export


write.fasta <- function (biom, outfile) {
  
  if (!is(biom, 'BIOM'))
    return (simpleError('In write.fasta(), biom must be a BIOM-class object.'))
  
  if (!is.character(outfile) || !identical(nchar(outfile) >= 3, TRUE) || length(outfile) != 1)
    return (simpleError('In write.fasta(), outfile must be a single filename at least 3 characters in length.'))
  
  seqs <- rbiom::sequences(biom)
  if (length(seqs) == 0)
    return (simpleError("There are no sequences in this BIOM object."))
  
  
  suffix <- substr(outfile, nchar(outfile) - 2, nchar(outfile))
  if        (identical(suffix, ".gz")) { conn <- gzfile(outfile, "w")
  } else if (identical(suffix, "bz2")) { conn <- bzfile(outfile, "w")
  } else                               { conn <-   file(outfile, "w") }
  
  cat(sprintf(">%s\n%s", names(seqs), seqs), file=conn, sep="\n")
  
  close(conn)
  
  
  return (invisible(NULL))
}