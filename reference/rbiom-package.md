# rbiom: Read/Write, Transform, and Summarize BIOM Data

A toolkit for working with Biological Observation Matrix (BIOM) files.
Features include reading/writing all BIOM formats, rarefaction, alpha
diversity, beta diversity (including UniFrac), summarizing counts by
taxonomic level, and sample subsetting. Standalone functions for
reading, writing, and subsetting phylogenetic trees are also provided.
All CPU intensive operations are encoded in C with multi-thread support.

## Multithreading

Many rbiom functions support multithreading:

The default behavior of these function is to run on as many cores as are
available in the local compute environment. If you wish to limit the
number of simultaneous threads, set `RcppParallel`'s `numThreads`
option. For instance:

        RcppParallel::setThreadOptions(numThreads = 4)

## See also

Useful links:

- <https://cmmr.github.io/rbiom/>

- <https://github.com/cmmr/rbiom>

- Report bugs at <https://github.com/cmmr/rbiom/issues>

## Author

**Maintainer**: Daniel P. Smith <dansmith01@gmail.com>
([ORCID](https://orcid.org/0000-0002-2479-2044))

Other contributors:

- Alkek Center for Metagenomics and Microbiome Research \[copyright
  holder, funder\]
