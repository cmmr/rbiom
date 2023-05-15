---
title: 'rbiom: An R package to for microbiome analysis'
authors:
  - name: Daniel P Smith
    orcid: 0000-0002-2479-2044
    corresponding: true
    affiliation: 1
  - name: Joseph F Petrosino
    orcid: 0000-0002-4046-6898
    affiliation: 1
affiliations:
 - name: Alkek Center for Metagenomics and Microbiome Research, Department of Molecular Virology and Microbiology, Baylor College of Medicine, Houston, Texas, USA
   index: 1
date: 21 April 2023
bibliography: paper.bib
---

# Summary

Microbes live all around us, on us, and even inside our bodies. Their influence
on health and disease is profound, and only beginning to be fully understood.
Studying microbial populations is becoming easier with modern DNA sequencing
technology; examining trends across thousands of samples today is common. The
bottleneck is longer collecting data, but rather analyzing and interpreting the
results.

`rbiom` is an R package for working with abundance datasets, such as OTU or ASV 
counts from 16S amplicon sequencing. It enables importing/exporting all BIOM
formats, subsetting, rarefying, manipulation of metadata/taxonomy/phylogeny,
computation of alpha and beta diversity metrics, and summarizing counts per
taxonomic rank. Computationally intensive tasks (including UniFrac [@unifrac]) 
have been implemented with multithreaded C++ to greatly reduce calculation 
time.

Visualization is a key component of `rbiom`. Rarefaction curves, taxa 
abundances, alpha diversity, and beta diversity can all be plotted in a variety 
of graphical formats, including correlation, heatmap, ordination, stacked bar, 
and box plots. In `rbiom`, box plots can be any combination of box, bar, 
violin, dot, strip, and/or range layers. Each plot includes provenance and 
modification history as attributes, as well as the `ggplot2` [@ggplot2] call 
used to render it to encourage downstream user customization.

Correlations between sample metadata and microbiome structure can be identified 
by mapping one or more metadata variables of interest to a plot's axes, facets, 
and/or aesthetics. These mappings can optionally define color/shape/pattern 
assignments, category ordering, or subsetting parameters. When metadata is 
associated with a axis or aesthetic, `rbiom` will automatically run the 
appropriate statistical test, correct for multiple comparisons, and display 
significant differences on the plot, captioning it with a brief methodology. 

Currently, `rbiom` can perform four types of significance testing. On 
correlation plots with a numeric metadata variable on the x-axis (e.g., Age, 
BMI), linear regression will be computed with R's `lm` linear model function.
For plots with two categories (e.g. Male vs Female), a Mann-Whitney test 
[@Mann1947] is run with R's `wilcox.test`. When three or more categories are 
compared, the Kruskal-Wallis rank sum test [@Kruskal1952] is used instead via 
R's `kruskal.test` function. P-values for ordinations are derived using the 
`adonis2` function from the `vegan` R package [@vegan], which randomly 
re-categorizes samples 1,000 times to estimate the significance of the observed 
clustering. P-values are corrected for multiple comparisons using the method 
described by @Benjamini1995 via R's `p.adjust` function to control for the 
false discovery rate. 

`QIIME2` [@qiime2], `mothur` [@mothur], and `Phyloseq` [@phyloseq] offer 
overlapping functionality with rbiom, but with important distinctions. The 
first two are designed for command-line interaction, making them difficult to 
integrate into R projects. `Phyloseq` has been a staple of R bioinformatics for 
a decade, but is frustratingly slow for studies with thousands of samples.

This package is designed for users of all experience levels. Novice R users
will appreciate that a couple commands will produce publication-ready figures.
Advanced R users can use `rbiom` to complement their existing pipelines with
faster and more flexible functions. `rbiom` is cross-platform compatible and
available from CRAN and conda-forge. The latest development version is on
GitHub.


# References
