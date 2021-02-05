# rbiom 1.0.2.9006

* Can now disable rarefaction in write.xlsx().
* Support for 'Decontam' generated BIOM files.
* Rarefy now handles samples and taxa with zero observations.
* write.xlsx exports 'dist' and 'matrix' attributes in addition to 'data.frame's
* Optionally change e.g. "Bacteria; Gracilibacteria; c" into "Bacteria; Gracilibacteria; Phylum Gracilibacteria (c)"
* Fixed crash when rarefying biom files and last taxa is not dropped first.


# rbiom 1.0.2

* Improved compatibility with Debian and Solaris.
* The 'rhdf5' package is now an optional dependency.
* The select() and subset() functions now subset sequences too.


# rbiom 1.0.0

* Initial Release
