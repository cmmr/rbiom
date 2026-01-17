# rbiom 3.0.0

## Announcements

The code for `rbiom`'s very fast implementations of diversity metrics,
rarefaction, and tree importing has been spun off into its own dependency-free
package, [`ecodive`](https://github.com/cmmr/ecodive), to encourage adoption by
other R packages. This also means `rbiom` will be easier to install going
forward as it will no longer need to be compiled.


## Breaking Changes

* `weighted` and `normalized` arguments are no longer used by any functions.
* Diversity metric names will now be all lowercase in returned matrices, tables, etc.

The following functions have been deprecated since `rbiom` v2.0 and are now
removed completely. They are generally functions which have been renamed for
improved organization.

`alpha.div`, `as.percent`, `beta.div`, `comments`, `comments<-`, `counts`,
`counts<-`, `depth`, `depths_barplot`, `has.phylogeny`, `has.sequences`, `id`,
`id<-`, `info`, `is.rarefied`, `metadata`, `metadata<-`, `nsamples`, `ntaxa`,
`phylogeny`, `phylogeny<-`, `rarefy`, `read.biom`, `read_biom`, `read.fasta`,
`read.tree`, `repair`, `sample.names`, `sample.names<-`, `sample_subset`,
`sample.sums`, `sequences`, `sequences<-`, `subtree`, `taxa.means`,
`taxa.names`, `taxa.names<-`, `taxa.ranks`, `taxa.ranks<-`, `taxa.rollup`,
`taxa.sums`, `taxonomy`, `taxonomy<-`, `tips`, `top.taxa`, `top_taxa`,
`unifrac`, `write.biom`, `write.fasta`, `write.tree`, `write.xlsx`


## Additions

* `adiv_vector()` returns a simple named vector of alpha diversity values.
* Replaced `rhdf5` dependency with `h5lite`.
* Updated list of alpha diversity metrics: `c("ace", "berger", "brillouin", "chao1", "faith", "fisher", "simpson", "inv_simpson", "margalef", "mcintosh", "menhinick", "observed", "shannon", "squares")`
* Updated list of beta diversity metrics: `c("aitchison", "bhattacharyya", "bray", "canberra", "chebyshev", "chord", "clark", "sorensen", "divergence", "euclidean", "generalized_unifrac", "gower", "hamming", "hellinger", "horn", "jaccard", "jensen", "jsd", "lorentzian", "manhattan", "matusita", "minkowski", "morisita", "motyka", "normalized_unifrac", "ochiai", "psym_chisq", "soergel", "squared_chisq", "squared_chord", "squared_euclidean", "topsoe", "unweighted_unifrac", "variance_adjusted_unifrac", "wave_hedges", "weighted_unifrac")`



# rbiom 2.2.1

* Compatibility fixes for ggplot2 4.0.0 #30


# rbiom 2.2.0

* Satisfy new R CMD check: no setting attributes on built-in functions. #28


# rbiom 2.1.2

* Fix for `write_biom(format='json')` where `{...}` should be `[...]`.
* Restore `rline` functionality in `rare_corrplot()`.


# rbiom 2.1.1

* `unifrac()` no longer normalizes weighted values, for back-compatibility.
* Option `underscores` allows keeping underscores in tree IDs.


# rbiom 2.0.13

* Major release with significant new features.
* Plotting added.
* Statistics added.
* Caching added.
* Clearer naming scheme for functions. Previous names still work but are deprecated.
* Generation of human-readable ggplot2 code for modifying plots outside of rbiom.
* Provenance tracking for BIOM objects and derivatives.


# rbiom 1.0.3

* Corrects for breaking changes in 'rhdf5' package.


# rbiom 1.0.2

* Improved compatibility with Debian and Solaris.
* The 'rhdf5' package is now an optional dependency.
* The select() and subset() functions now subset sequences too.


# rbiom 1.0.0

* Initial Release
