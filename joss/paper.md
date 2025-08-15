---
title: 'rbiom: An R Package for Streamlined Microbiome Analysis and Automated Statistical Visualization'
authors:
- name: Daniel P Smith
  orcid: "0000-0002-2479-2044"
  corresponding: true
  affiliation: 1, 2
- name: Sara J Javornik Cregeen
  orcid: "0009-0000-2698-6478"
  affiliation: 1, 2
- name: Joseph F Petrosino
  orcid: "0000-0002-4046-6898"
  affiliation: 1, 2
date: "15 August 2025"
output:
  word_document: default
  pdf_document: default
tags:
- bioinformatics
- bacteria
- ecology
- metadata
- statistics
- visualization
- regression
- correlation
- heatmap
- ordination
- community
- reproducible
- unifrac
affiliations:
- name: The Alkek Center for Metagenomics and Microbiome Research, Department of Molecular
    Virology and Microbiology, Baylor College of Medicine, Houston, TX 77030, USA
  index: 1
- name: Department of Molecular Virology and Microbiology, Baylor College of Medicine,
    Houston, TX, USA
  index: 2
  ror: 02pttbw34
bibliography: paper.bib
---



# Summary

`rbiom` is a software toolkit for the R programming language that simplifies the analysis of complex ecological data. Its main purpose is to help researchers uncover hidden trends and patterns by connecting what they observe in a biological sample (like the types and numbers of bacteria) with information they know about that sample (like where it was collected or the health of the individual). `rbiom` makes this process easy by focusing on creating clear, publication-quality figures that include statistical test results right on the chart. This allows scientists to quickly visualize and understand the important relationships within their data.


# Statement of Need

A core challenge in microbiome analysis is identifying meaningful associations between a community's composition and its associated metadata. While established packages like `phyloseq` [@phyloseq] exist for handling the large and complex datasets often stored in Biological Observation Matrix (BIOM) files, there is a persistent need for a tool that streamlines the entire workflow from data ingestion to the creation of publication-quality figures with integrated statistical results. `rbiom` addresses this need by providing a unified interface that connects the analysis of taxonomic abundance and diversity with powerful statistical methods and visualizations. Its design prioritizes speed and reproducibility, making it particularly suitable for both routine data exploration and rigorous statistical reporting.



# Key Features

The core design principles of rbiom focus on providing a seamless and reproducible workflow for microbiome data analysis. Key features include:

* **Integrated statistical visualization:** Statistical tests are integrated directly into plotting functions, automatically annotating figures with results like p-values.

* **Reproducibility:** The package provides the underlying R code used to generate figures and statistics, allowing for easy customization and extension.

* **Workflow simplification:** A unified interface connects data import, manipulation, analysis, and visualization in a single toolkit.



# Related Works

A recent review by @Wen2023 highlighted the vast number of R packages available for microbiome analysis. While many of these tools offer some level of statistical capability, few make a concerted effort to display statistical results directly on generated figures. For example, the foundational `phyloseq` package lacks built-in significance testing. The closest peer to `rbiom` is `microeco` [@microeco], which offers a comprehensive array of statistical tests and integrates their output into figures through the `ggpubr` package [@ggpubr]. `rbiom` extends this functionality by automatically positioning significance brackets, displaying the statistical method on the plot, and providing the underlying R code used to generate the figures and statistics. This approach enhances reproducibility and allows users to easily customize or extend the analysis.


# Functionality

The `rbiom` package offers the following key functionalities:

* **Data Management and Manipulation:** `rbiom` seamlessly imports various data formats, including BIOM files, `QIIME2` [@QIIME2] and `mothur` [@mothur] outputs, and objects from other popular R packages using the `phyloseq` or `SummarizedExperiment` classes. The package includes a robust set of tools for rarefaction, filtering, and summarization, enabling users to prepare their data for downstream analysis.

* **Diversity and Composition Analysis:** The package focuses on three key community features: **alpha diversity**, **beta diversity**, and **taxa abundance**. For each, it can compute statistics against categorical or continuous metadata using appropriate methods such as **Kruskal-Wallis**, **Mann-Whitney**, **PERMANOVA**, and **estimated marginal means** [@Kruskal1952; @Mann1947; @Anderson2001; @emmeans].

* **Statistical Visualization:** A core feature of `rbiom` is its ability to directly overlay statistical test results onto `ggplot2`-based figures [@ggplot2]. Functions like `adiv_boxplot()` and `stats_corrplot()` generate visualizations with automated annotations for p-values, trend lines, confidence intervals, and methodology. Customization is supported through parameters like `p.label` to display only significant results, ensuring clean, publication-ready graphics.

* **High-Level Plotting:** The package provides a diverse set of highly customizable plot types, including **boxplots**, **heatmaps**, **stacked bar charts**, and **ordination plots**. Users have extensive control over aesthetics, with options for specifying metadata variables for the x-axis, statistical groups, and plot faceting.



# Example Figures

``` r
library(rbiom)
biom <- rarefy(hmp50) # hmp50 dataset is included with rbiom
bdiv_ord_plot(biom, bdiv = "UniFrac", stat.by = "Body Site", facet.by = "Sex")
```

![A beta diversity ordination plot. Samples cluster significantly by body site (p = 0.001) and are characterized by different bacterial genera.](figures/bdiv_ord_plot.png)

``` r
adiv_boxplot(biom, x = "Sex", adiv = c("simp", "shan"), stat.by = "Body Site")
```

![An alpha diversity box plot. Observed OTUs and shannon diversity indices vary significantly by body site for both males (p = 2e-04) and females (p = 0.003).](figures/adiv_boxplot.png)


``` r
subset(biom, `Body Site` == 'Buccal mucosa') %>% 
  taxa_corrplot("Age", taxa = 2, layers = 'ptc', fit = 'lm', test = 'emtrends') +
  ggplot2::theme_classic()
```

![A taxa correlation plot using an alternative theme from ggplot2. The two most abundant buccal mucosa-associated genera show weak correlations with age.](figures/taxa_corrplot.png)



# Conclusion

By unifying data management, statistical analysis, and automated visualization, rbiom provides a powerful and accessible tool for microbiome research. The package's focus on reproducible, publication-ready graphics with built-in statistical results makes it a valuable addition to the ecosystem of R packages for ecological and biological data science.



# Acknowledgements

This study was supported by NIH/NIAD (Grant number U19 AI44297), and Baylor
College of Medicine and Alkek Foundation Seed.

The authors would like to thank Gemini for its assistance in drafting and
refining this manuscript.


# References
