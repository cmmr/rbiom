#!/bin/env Rscript


#============================================================================================
#=-------------------------------------------------------------------------------------------
# Daniel Smith
# May 1st, 2017
# Baylor College of Medicine
# ==========================
#
#   A standalone script for running UniFrac from the command line via R and rbiom.
#
#=-------------------------------------------------------------------------------------------
#============================================================================================




#========================================================
# Parameters provided by the user.
#========================================================
cat("Parsing command line options.\n")

opt_parser = optparse::OptionParser(

  usage = "\n\n  %prog -i in.biom -t in.tre -m weighted -o out.txt",

  description = "
  Creates a UniFrac distance matrix from a provided biom and tree file.",

  option_list=list(
    optparse::make_option(c("-i", "--infile"), type="character", default=NULL,
                          help="Input file in BIOM format."),

    optparse::make_option(c("-t", "--tree"), type="character", default=NULL,
                          help="Phylogenetic tree in newick format. Will attempt to load from BIOM."),

    optparse::make_option(c("-m", "--method"), type="character", default="weighted",
                          help="Options are 'weighted' and 'unweighted' [default: %default]"),

    optparse::make_option(c("-o", "--outfile"), type="character", default="distmat.txt",
                          help="Where to write distance matrix. [default: %default]")
  ),

  epilogue = "Author: Daniel Smith, 2017 Baylor College of Medicine\n\n"
)
opt = optparse::parse_args(opt_parser)

if (is.null(opt$infile)) {
  optparse::print_help(opt_parser)
  stop("An input BIOM file is required.\n\n", call.=FALSE)
}

if (!file.exists(opt$infile))
  stop(sprintf("The specified file (%s) does not exist.\n\n", opt$file), call.=FALSE)

if (!opt$method %in% c('weighted', 'unweighted'))
  stop("Method must be 'weighted' or 'unweighted'.\n\n", call.=FALSE)


library(rbiom)

biom     <- read.biom(opts$infile)
weighted <- ifelse(identical(opts$method, "weighted"), TRUE, FALSE)
dm       <- unifrac(biom, weighted, opts$tree)
write.table(dm, opts$outfile, sep="\t", quote=FALSE)


