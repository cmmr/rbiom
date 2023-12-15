#' Speed Ups.
#' 
#' When working with very large datasets, you can make use of these tips and 
#' tricks to speed up operations on rbiom objects.
#' 
#' @name speed_ups
#' @keywords internal
#' 
#' 
#' @section Skip Cloning:
#' 
#' Functions that modify rbiom objects, like [subset()] and [rarefy()], will
#' automatically clone the object before modifying it. This is to make these 
#' functions behave as most R users would expect - but at a performance trade 
#' off.
#' 
#' Rather than:
#' ```r
#' biom <- subset(biom, ...)
#' biom <- rarefy(biom)
#' ```
#' 
#' Modify `biom` in place like this:
#' ```r
#' subset(biom, clone = FALSE, ...)
#' rarefy(biom, clone = FALSE)
#' 
#' # Or:
#' biom$metadata %<>% subset(...)
#' biom$counts %<>% rarefy_cols()
#' ```
#' \cr
#' 
#' 
#' 
#' @section Drop Components:
#' 
#' ## Sequences
#' 
#' Reference sequences for OTUs will be imported along with the rest of your 
#' dataset and stored in `$sequences`. However, rbiom doesn't currently use 
#' these sequences for anything (except writing them back out with 
#' [write_biom()] or [write_fasta()]).
#' 
#' You can delete them from your rbiom object with:
#' 
#' ```r
#' biom$sequences <- NULL
#' ```
#' 
#' 
#' ## Tree
#' 
#' The phylogenetic reference tree for OTUs is only used for calculating 
#' UniFrac distances. If you aren't using UniFrac, the tree can be dropped 
#' from the rbiom object with:
#' 
#' ```r
#' biom$tree <- NULL
#' ```
#' 
#' Alternatively, you can store the tree separately from the rbiom object and 
#' provide it to just the functions that use it. For example:
#' 
#' ```r
#' tree <- biom$tree
#' biom$tree <- NULL
#' dm <- bdiv_distmat(biom, 'unifrac', tree = tree)
#' ```
#' \cr
#' 
#' 
#' 
#' @section Increase Caching:
#' 
#' Caching is enabled by default - up to 20 MB per R session.
#' 
#' For large datasets, increasing the cache size can help. The size is 
#' specified in bytes by an R option or environment variable.
#' 
#' ```r
#' options(rbiom.cache_size=200 * 1024 ^ 2) # 200 MB
#' Sys.setenv(RBIOM_CACHE_SIZE=1024 ^ 3)    # 1 GB
#' ```
#' \cr
#' 
#' You can also specify a cache directory where results can be preserved from 
#' one R session to the next.
#' 
#' ```r
#' options(rbiom.cache_dir=tools::R_user_dir("rbiom", "cache"))
#' Sys.setenv(RBIOM_CACHE_DIR="~/rbiom_cache")
#' ```
#' \cr
#' 
#' Other quick notes about caching:
#' 
#' * Setting the cache directory to `"FALSE"` will disable caching.
#' * R options will override environment variables.
#' * The key hash algorithm can be set with `options(rbiom.cache_hash=rlang::hash)`.
#' 
NULL
