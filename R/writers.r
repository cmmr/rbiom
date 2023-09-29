


#' Common code for writers.
#' @noRd
#'
#' @param outfile   The file to open a connection to.
#' 
#' @param callback   A function that takes a connection as its only
#'        argument and writes content to it.
#' 
#' @return The normalized filepath that was written to (invisibly).
#' 
write_wrapper <- function (outfile, callback) {
  
  stopifnot(is_scalar_character(outfile) && !is_na(outfile))
  stopifnot(isTRUE(nchar(outfile) > 0))
  
  
  outfile <- normalizePath(outfile, winslash = "/", mustWork = FALSE)
  
  tryCatch(
    error = function (e) stop("Can't write to '", outfile, "'.\n", e),
    expr  = local({
      
      if (!dir.exists(dirname(outfile)))
        dir.create(dirname(outfile), recursive = TRUE)
      
      con <- if (endsWith(outfile, ".gz"))  { gzfile(outfile, "w")
      } else if (endsWith(outfile, ".bz2")) { bzfile(outfile, "w")
      } else                                { file(outfile, "w") }
      
      callback(con)
      
      close(con)
    }))
  
  return (invisible(outfile))
  
  
}




#' Write BIOM metadata to a tab-separated value (tsv) file.
#'
#' @param biom   A BIOM object
#' 
#' @param file   Path to the output file. Filenames ending in 
#'        \code{.gz} or \code{.bz2} will be compressed accordingly.
#' 
#' @return The normalized filepath that was written to (invisibly).
#' 
#' @export

write_metadata <- function (biom, file) {
  
  stopifnot(is(biom, 'BIOM'))
  
  write_wrapper(file, function (con) {
    md     <- biom[['metadata']]
    header <- c("SampleID", colnames(md))
    cat(file = con, paste(collapse = "\t", header), "\n")
    write.table(
      x         = md, 
      file      = con,
      col.names = FALSE,
      append    = TRUE,
      quote     = FALSE, 
      sep       = "\t" )
  })
}




#' Write BIOM counts to a tab-separated value (tsv) file.
#'
#' @param biom   A BIOM object
#' 
#' @param file   Path to the output file. File names ending in 
#'        \code{.gz} or \code{.bz2} will be compressed accordingly.
#' 
#' @return The normalized filepath that was written to (invisibly).
#' 
#' @export

write_counts <- function (biom, file) {
  write_biom(as.matrix(biom[['counts']]), file)
}




#' Write BIOM taxonomy map to a tab-separated value (tsv) file.
#'
#' @param biom   A BIOM object
#' 
#' @param file   Path to the output file. Filenames ending in 
#'        \code{.gz} or \code{.bz2} will be compressed accordingly.
#' 
#' @return The normalized filepath that was written to (invisibly).
#' 
#' @export

write_taxonomy <- function (biom, file) {
  
  stopifnot(is(biom, 'BIOM'))
  
  write_wrapper(file, function (con) {
    write.table(
      x         = biom[['taxonomy']], 
      file      = con, 
      quote     = FALSE, 
      sep       = "\t", 
      col.names = FALSE )
  })
}




#' Write DNA sequences to a file in fasta format.
#'
#' @param seqs   A named character vector where names are sequence names and
#'        values are the sequences. Also accepts a \code{BIOM} object which 
#'        contains sequences.
#' 
#' @param file   Path to the output fasta file. Filenames ending in 
#'        \code{.gz} or \code{.bz2} will be compressed accordingly.
#' 
#' @return The normalized filepath that was written to (invisibly).
#' 
#' @export

write_fasta <- function (seqs, file = NULL) {
  
  
  if (is(seqs, 'BIOM'))
    seqs <- seqs[['sequences']]
  
  if (!is(seqs, 'character') || !is(names(seqs), 'character'))
    stop("In write_fasta(), seqs must be a named character vector or a BIOM-class object.")
  
  if (length(seqs) == 0)
    warning("In write_fasta(), seqs is empty.")
  
  
  write_wrapper(file, function (con) {
    cat(file=con, sep="\n", sprintf(">%s\n%s", names(seqs), seqs))
  })
}




#' Write a newick formatted phylogenetic tree.
#' 
#' @param tree  A \code{phylo} object, as returned from [read_tree()]. Also 
#'        accepts a \code{BIOM} object if it has a phylogentic tree.
#'         
#' @param file  Filename or connection to write the newick file to (optional).
#' 
#' @return If file is NULL, the newick string as a character vector. Otherwise,
#'         the normalized filepath that was written to (invisibly).
#'         
#' @export
#' @examples
#'     library(rbiom) 
#'     
#'     infile <- system.file("extdata", "newick.tre", package = "rbiom")
#'     tree   <- read_tree(infile)
#'     newick <- write_tree(tree)
#'
write_tree <- function (tree, file=NULL) {
  
  if (is(tree, "BIOM"))
    tree <- tree[['phylogeny']]
  
  if (!is(tree, "phylo"))
    stop("Provided tree is not a 'phylo' or 'BIOM' class object.")
  
  
  rootNode <- setdiff(tree$edge[,1], tree$edge[,2])
  parentAt <- aggregate(1:nrow(tree$edge), by=list(tree$edge[,1]), c, simplify=FALSE)
  parentAt <- setNames(lapply(parentAt[,2], unlist), parentAt[,1])
  
  fx <- function (root=NULL) ({
    
    nodes <- parentAt[[as.character(root)]]
    
    if (length(nodes) == 0) {
      
      nodeLabel <- tree$tip.label[root]
      
      if (any(grepl(" ", nodeLabel, fixed=TRUE))) {
        if (any(grepl("_", nodeLabel, fixed=TRUE))) {
          nodeLabel <- paste0("'", nodeLabel, "'")
        } else {
          nodeLabel <- gsub(" ", "_", nodeLabel)
        }
      }
      return (nodeLabel)
    }
    
    children <- tree$edge[nodes, 2]
    children <- sapply(children, fx)
    
    if (!is_null(tree$edge.length))
      children <- paste(sep=":", children, tree$edge.length[nodes])
    
    sprintf("(%s)", paste(collapse=",", children))
  })
  
  newick <- paste0(fx(rootNode), ";")
  
  
  if (is_null(file)) return (newick)
  
  
  write_wrapper(file, function (con) {
    writeLines(text=newick, con=con, sep="")
  })
  
}

