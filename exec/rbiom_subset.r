#!/usr/bin/env Rscript


#============================================================================================
#=-------------------------------------------------------------------------------------------
# Daniel Smith
# October 26th, 2022
# Baylor College of Medicine
# ==========================
#
#   Subsets a biom file, either using an include/exclude sample ID list, or by metadata.
#
#=-------------------------------------------------------------------------------------------
#============================================================================================


if (nchar(system.file(package = "optparse")) == 0)
  stop("CRAN R package 'optparse' must be installed to use this script.")


opt_parser <- optparse::OptionParser(
  
  usage = "\n   %prog -i input.biom -o output.biom [options]",
  
  option_list=list(
    
    optparse::make_option(c("-i", "--infile"), type="character", default=NULL,
      help="Existing biom file to subset. Can be a URL."),
    
    optparse::make_option(c("-o", "--outfile"), type="character", default=NULL,
      help="Where to save the new biom file."),
    
    optparse::make_option(c("-f", "--field"), type="character", default="SampleID",
      help="Metadata field to use when --keep and/or --drop is specified. [default: %default]"),
    
    optparse::make_option(c("-k", "--keep"), type="character", default=NULL,
      help="Values to keep. Either comma-separated list, or a filename containing one value per line."),
    
    optparse::make_option(c("-d", "--drop"), type="character", default=NULL,
      help="Values to drop. Either comma-separated list, or a filename containing one value per line."),
    
    optparse::make_option(c("-m", "--metadata"), type="character", default=NULL,
      help="A metadata file to import into the new biom file. First column must be SampleID. 
            Samples missing from this file will be excluded from the output biom file.
            Metadata file format is assumed to be tab-delimited except when file's extension is '.csv'.")
  ),
  
  epilogue = '
Subsets a biom file by sample IDs or metadata.

Examples:
  
  SampleID Filtering:
  ------------------------------------------
  rbiom_subset.r -i hmp50.biom -o 02.biom -k HMP01,HMP02
  rbiom_subset.r -i hmp50.biom -o 48.biom -d HMP01,HMP02
  
  
  When the SampleIDs are in a File:
  ------------------------------------------
  cat << EOF > sample_ids.txt
  HMP01
  HMP02
  HMP03
  EOF
  
  rbiom_subset.r -i hmp50.biom -o subset.biom -k sample_ids.txt
  
  
  Metadata Filtering:
  ------------------------------------------
  rbiom_subset.r -i hmp50.biom -o males.biom   -f Sex -k Male
  rbiom_subset.r -i hmp50.biom -o females.biom -f Sex -d Male
  
  
  When the Metadata is in a File:
  ------------------------------------------
  cat << EOF > metadata.txt
  SampleID  AB
  HMP01     A
  HMP02     B
  HMP03     A
  HMP04     B
  EOF
  
  rbiom_subset.r -i hmp50.biom -o 1234.biom -m metadata.txt
  rbiom_subset.r -i hmp50.biom -o A.biom    -m metadata.txt -f AB -k A
  rbiom_subset.r -i hmp50.biom -o B.biom    -m metadata.txt -f AB -d A
  

Author: Daniel Smith, 2022 Baylor College of Medicine
'
)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$infile) || is.null(opt$outfile)) {
  optparse::print_help(opt_parser)
  stop("--infile and --outfile are required.\n\n", call.=FALSE)
}

if (is.null(opt$keep) && is.null(opt$drop) && is.null(opt$metadata)) {
  optparse::print_help(opt_parser)
  stop("At least one of --keep, --drop, or --metadata is required.\n\n", call.=FALSE)
}


#________________________________________________________
# Import the biom file.
#________________________________________________________
library(rbiom)

cat("\n")
cat("Reading biom file '", opt$infile, "'\n", sep = '')
biom <- read_biom(src = opt$infile)

cat("Initial sample count: ", nsamples(biom), "\n", sep = '')


#________________________________________________________
# Make a comma-separated list for error messages.
#________________________________________________________
oxford <- function (vals) {
  vals <- paste0("'", vals, "'")
  if (length(vals) == 1) return (vals)
  if (length(vals) == 2) return (paste(vals, collapse = " and "))
  if (length(vals) > 6)
    vals <- c(vals[1:4], sprintf("%i more", length(vals) - 4))
  vals[[length(vals)]] <- sprintf("and %s", vals[[length(vals)]])
  return (paste(vals, collapse = ", "))
}


#________________________________________________________
# Optionally add new metadata to the biom object.
#________________________________________________________
if (!is.null(f <- opt$metadata)) {
  biom <- local({
    
    new_md <- if (endsWith(tolower(f), ".csv")) {
      utils::read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
    } else {
      utils::read.delim(f, check.names=FALSE, stringsAsFactors=FALSE)
    }
    
    
    # Check for duplicate IDs in metadata file.
    #________________________________________________________
    if (length(i <- which(duplicated(new_md[[1]]))) > 0)
      stop("Duplicated sample IDs in first column of metadata file: ", oxford(new_md[[1]][i]))
    rownames(new_md) <- new_md[[1]]
    new_md  <- new_md[,seq_len(ncol(new_md) - 1) + 1,drop=FALSE]
    new_ids <- rownames(new_md)
    
    
    # Intersect the old and new sample IDs
    #________________________________________________________
    old_ids <- sample_names(biom)
    
    if (length(x <- setdiff(new_ids, old_ids)) > 0)
      warning("Sample IDs from metadata file are not in input biom file: ", oxford(x))
    
    keep_ids <- intersect(new_ids, old_ids)
    if (length(keep_ids) == 0)
      stop("No sample IDs in common between metadata file and input biom file.")
    
    drop_ids <- setdiff(old_ids, keep_ids)
    if (length(x <- setdiff(new_ids, old_ids)) > 0)
      cat("Metadata file excluded ", length(drop_ids), " samples from input biom file.\n", sep = '')
    
    biom <- select(biom, keep_ids)
    
    
    # Intersect the old and new metadata.
    #________________________________________________________
    md <- metadata(biom)
    for (i in colnames(new_md))
      md[[i]] <- new_md[rownames(md), i]
    metadata(biom) <- md
    
    return (biom)
  })
}


#________________________________________________________
# Sanity check metadata field to limit by.
#________________________________________________________
if (!identical(opt$field, "SampleID")) {
  if (!opt$field %in% colnames(metadata(biom)))
    stop("Metadata field ", opt$field, " does not exist.")
}


#________________________________________________________
# Keep/drop particular SampleIDs/values
#________________________________________________________
if (!is.null(opt$keep) || !is.null(opt$drop)) {
  
  biom <- local({
    
    # Read SampleIDs/values from csv string or file.
    #________________________________________________________
    for (v in c("keep", "drop")) {
      vals <- opt[[v]]
      if (is.null(vals))            { vals <- NULL
      } else if (file.exists(vals)) { vals <- readLines(con = vals)
      } else                        { vals <- strsplit(vals, ",")[[1]] }
      assign(v, vals[nzchar(trimws(vals))])
      remove("vals")
    }
    
    
    # Convert metadata values to SampleIDs
    #________________________________________________________
    if (!identical(opt$field, "SampleID")) {
      ids  <- sample_names(biom)
      vals <- as.character(metadata(biom, opt$field))
      
      if (length(x <- setdiff(c(keep, drop), vals)) > 0)
        warning(opt$field, " values are never used: ", oxford(x))
      
      if (!is.null(keep)) keep <- ids[vals %in% keep]
      if (!is.null(drop)) drop <- ids[vals %in% drop]
      
      remove("ids", "vals")
    }
    
    
    # Drop SampleIDs from the biom.
    #________________________________________________________
    ids <- sample_names(biom)
    if (length(x <- setdiff(c(keep, drop), ids)) > 0)
      warning("Sample IDs are not in input files: ", oxford(x))
    
    if (!is.null(keep)) ids <- intersect(ids, keep)
    if (!is.null(drop)) ids <- setdiff(ids, drop)
    
    if (length(ids) == 0)
      stop("No samples retained.")
    
    biom <- select(biom = biom, samples = ids)
    return (biom)
  })
}


if (length(nsamples(biom)) == 0)
  stop("No samples retained.")

cat("Final sample count: ", nsamples(biom),   "\n",  sep = '')
cat("Saving new biom file to '", opt$outfile, "'\n", sep = '')

write_biom(biom = biom, file = opt$outfile)

cat("Finished.\n\n", sep = '')


