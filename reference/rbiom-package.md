# rbiom: Integrated Analysis and Visualization of Microbiome Data

A toolkit for working with Biological Observation Matrix (BIOM) files.
Features include reading/writing all BIOM formats, rarefaction, alpha
diversity, beta diversity (including UniFrac), summarizing counts by
taxonomic level, sample subsetting, statistics, and visualization.
Standalone functions for reading, writing, and subsetting phylogenetic
trees are also provided.

## Reproducibility

Every rbiom object records the steps that produced it, and every table,
statistic, and plot carries the R code needed to regenerate it. See
`biom$code`, `table$code`, `stats$code`, and `plot$code`. The [rbiom
objects](https://cmmr.github.io/rbiom/reference/rbiom_objects.md) page
has details.

## Multithreading

Diversity calculations are performed by the 'ecodive' package, which
runs on all available cores by default. Use the `cpus` argument of
individual functions to limit the number of threads, e.g.
`adiv_table(biom, cpus = 1)`.

The linear algebra library (BLAS) that R is linked to may run threads of
its own: OpenBLAS, MKL, and BLIS start up to one per core for operations
such as the eigendecomposition behind a PCoA, whatever `cpus` says. So
when `cpus` is less than the number of threads the library is set to
use, rbiom runs it on one thread until the call returns, and then puts
the previous setting back, also after an error or an interrupt. One
thread rather than `cpus`: in rbiom's computations extra BLAS threads
were at most about 20% faster on an idle machine, and slower on a busy
one or for model fits. With the default `cpus`, all available cores, the
library is left as it is. The number of BLAS threads can change results
only by rounding in the last digits (around 1e-15).

This needs a version of 'ecodive' that provides `blas_threads()`. R's
reference BLAS (the default on Windows and macOS) is single-threaded, so
there is nothing to limit; Apple's Accelerate can't be limited this way.
To limit the library for a whole R process instead, set an environment
variable such as `OPENBLAS_NUM_THREADS=1` before R starts.

## See also

Useful links:

- <https://cmmr.github.io/rbiom/>

- <https://github.com/cmmr/rbiom>

- Report bugs at <https://github.com/cmmr/rbiom/issues>

## Author

**Maintainer**: Daniel P. Smith <dansmith01@gmail.com>
([ORCID](https://orcid.org/0000-0002-2479-2044))

Authors:

- Daniel P. Smith <dansmith01@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-2479-2044))

Other contributors:

- Alkek Center for Metagenomics and Microbiome Research \[copyright
  holder, funder\]
