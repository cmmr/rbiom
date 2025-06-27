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
