
#' Parse a fasta file into a named character vector.
#'
#' @param file   A file path, URL, or literal string with fasta-formatted
#'        sequences. Files can optionally be compressed with gzip or bzip2.
#'
#' @param ids   Character vector of IDs to retrieve. The default, `NULL`,
#'        will retrieve everything.
#'
#' @return A named character vector in which names are the fasta headers
#'         (the text after `>` up to the first whitespace) and values are the
#'         sequences. Sequences spanning multiple lines are joined, and
#'         whitespace within sequences is removed.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     fasta <- ">OTU1 first sequence\nATCG\nGGCC\n>OTU2\nTTAA\n"
#'     read_fasta(fasta)
#'
#'     outfile <- write_fasta(hmp50, tempfile(fileext = ".fna.gz"))
#'     seqs    <- read_fasta(outfile)
#'     unlink(outfile)
#'
#'     identical(seqs, hmp50$sequences)
#'

read_fasta <- function (file, ids = NULL) {

  if (!is_scalar_character(file) || is_na(file))
    cli_abort("`file` must be a file path, URL, or fasta string, not {.type {file}}.")

  if (!is.null(ids) && !is.character(ids))
    cli_abort("`ids` must be NULL or a character vector, not {.type {ids}}.")


  #________________________________________________________
  # Fetch URLs; write literal strings to a temp file.
  #________________________________________________________
  fpath <- as_filepath(file)
  on.exit(fpath$cleanup(), add = TRUE)

  if (!identical(fpath$format, 'text'))
    cli_abort("{.file {fpath$id}} is not a fasta file ({fpath$format} format detected).")


  #________________________________________________________
  # Headers start with '>'; everything else is sequence.
  #________________________________________________________
  lines <- readLines(fpath$path, warn = FALSE, encoding = "UTF-8")
  lines <- sub("\r$", "", strip_bom(lines))
  lines <- lines[nzchar(trimws(lines))]

  is_hdr <- startsWith(lines, ">")

  if (!any(is_hdr))
    cli_abort("No fasta headers (lines beginning with {.val >}) found in {.file {fpath$id}}.")

  if (!is_hdr[[1]])
    cli_abort("Sequence data precedes the first fasta header in {.file {fpath$id}}.")

  n_hdr <- sum(is_hdr)
  grp   <- cumsum(is_hdr)
  nms   <- sub("\\s.*$", "", trimws(substring(lines[is_hdr], 2)))

  seqs <- gsub("[[:space:]]", "", lines[!is_hdr])
  seqs <- split(seqs, factor(grp[!is_hdr], levels = seq_len(n_hdr)))
  seqs <- vapply(seqs, paste, "", collapse = "", USE.NAMES = FALSE)

  res <- setNames(seqs, nms)


  #________________________________________________________
  # IDs must be unique.
  #________________________________________________________
  if (length(x <- unique(nms[duplicated(nms)])) > 0)
    cli_abort(c(
      '!' = "Fasta file {.file {fpath$id}}",
      'x' = 'Has multiple entries for {.val {x}}.' ))


  #________________________________________________________
  # Ensure we return as many IDs as requested
  #________________________________________________________
  if (!is_null(ids)) {

    if (length(x <- setdiff(ids, nms)))
      cli_abort(c(
        '!' = "In fasta file {.file {fpath$id}}",
        'x' = 'Expected IDs are missing: {.val {x}}.' ))

    res <- res[ids]
  }


  return (res)
}
