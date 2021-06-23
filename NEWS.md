# rbiom 1.0.2.9025

* Can now disable rarefaction in write.xlsx().
* Support for 'Decontam' and 'MicrobiomeDB' generated BIOM files.
* Rarefy now handles samples and taxa with zero observations.
* write.xlsx exports 'dist' and 'matrix' attributes in addition to 'data.frame's
* Optionally change e.g. "Bacteria; Gracilibacteria; c" into "Bacteria; Gracilibacteria; Phylum Gracilibacteria (c)"
* Fixed crash when rarefying biom files and last taxa is not dropped first.
* Easily pull a single column of metadata with metadata(biom, field).
* Added apcoa function for covariate adjusted principal coordinates analysis.
* Switched from rjson to jsonlite for better handling of non-UTF characters.
* Support for writing out biom files compressed with gzip or bzip2.
* Memoised the functions that do heavy-lifting.
* Set memoise to be off by default.
* Can now download files with odd characters in their name.
* Prevent NUL from being appended to read.biom's text inputs.
* alpha.div can now return metadata and/or subset of adiv metrics in wide or long format.
* Automagically remove NAs and update factors when using subset().
* taxa.rollup can now return with metadata in wide or long format.
* beta.div can now return with metadata in long format.
* Added stats.table function to compute p-vals for alpha, beta, and taxa metrics.
* Added sample.sums function to count the number of observations per sample.
* New functions: depth(), is.rarefied(), has.phylogeny(), and has.sequences().
* Added plot function - provides stat brackets, patterned fills, overlaid geoms and more.


# rbiom 1.0.2

* Improved compatibility with Debian and Solaris.
* The 'rhdf5' package is now an optional dependency.
* The select() and subset() functions now subset sequences too.


# rbiom 1.0.0

* Initial Release
