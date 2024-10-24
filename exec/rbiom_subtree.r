#!/usr/bin/env Rscript


#============================================================================================
#=-------------------------------------------------------------------------------------------
# Daniel Smith
# October 24th, 2024
# Baylor College of Medicine
# ==========================
#
#   Subsets a newick file to match IDs from a fasta file.
#
#=-------------------------------------------------------------------------------------------
#============================================================================================


if (nchar(system.file(package = "optparse")) == 0)
  stop("CRAN R package 'optparse' must be installed to use this script.")


opt_parser <- optparse::OptionParser(
  
  usage = "\n   %prog -i input.tre -f seqs.fna -o output.tre",
  
  option_list=list(
    
    optparse::make_option(c("-i", "--infile"), type="character", default=NULL,
                          help="Existing newick tree file to subset. Can be a URL."),
    
    optparse::make_option(c("-o", "--outfile"), type="character", default=NULL,
                          help="Where to save the new newick tree file."),
    
    optparse::make_option(c("-f", "--fasta"), type="character", default=NULL,
                          help="A fasta file with IDs matching those in 'infile'.")
  ),
  
  epilogue = '
Subsets a newick file to match IDs from a fasta file.

Example:
  
  rbiom_subtree.r -i silva.tre -f silva_V4.fna -o silva_V4.tre
  

Author: Daniel Smith, 2024 Baylor College of Medicine
'
)
opt <- optparse::parse_args(opt_parser)

for (i in c('infile', 'fasta', 'outfile'))
  if (is.null(opt[[i]])) {
    optparse::print_help(opt_parser)
    stop(paste("--", i, "is required.\n\n", call.=FALSE))
  }


library(rbiom)

tree <- read_tree(src = opt$infile)
seqs <- read_fasta(file = opt$fasta)
tree <- tree_subset(tree = tree, tips = names(seqs))
write_tree(biom = tree, file = opt$outfile)
