#' Write sequences from a BIOM object to a file in fasta format.
#'
#' @param seqs     A named character vector where names are sequence names and
#'                   values are the sequences. Also accepts a \code{BIOM}
#'                   object which contains sequences.
#' @param outfile  Path to the output fasta file. Files ending in \code{.gz} or
#'                   \code{.bz2} will be compressed accordingly.
#' @return On success, returns \code{NULL} invisibly.
#' 
#' @export


write.fasta <- function (seqs, outfile) {
  
  if (is(seqs, 'BIOM'))
    seqs <- rbiom::sequences(seqs)
  
  if (!is(seqs, 'character') || !is(names(seqs), 'character'))
    stop (simpleError('In write.fasta(), seqs must be a named character vector or a BIOM-class object.'))
  
  if (!is.character(outfile) || !identical(nchar(outfile) >= 3, TRUE) || length(outfile) != 1)
    stop (simpleError('In write.fasta(), outfile must be a single filename at least 3 characters in length.'))
  
  
  if (length(seqs) == 0)
    stop (simpleError("There are no sequences in this object."))
  
  
  suffix <- substr(outfile, nchar(outfile) - 2, nchar(outfile))
  if        (identical(suffix, ".gz")) { conn <- gzfile(outfile, "w")
  } else if (identical(suffix, "bz2")) { conn <- bzfile(outfile, "w")
  } else                               { conn <-   file(outfile, "w") }
  
  cat(sprintf(">%s\n%s", names(seqs), seqs), file=conn, sep="\n")
  
  close(conn)
  
  
  return (invisible(NULL))
}