test_that("rbiom-deprecated", {
  
  skip_on_cran()
  
  rlang::local_options(lifecycle_verbosity = "warning")
  
  biom <- hmp5$clone()
  
  expect_warning(file <- write.biom(biom,            tempfile()))
  expect_warning(xlsx <- write.xlsx(biom,            tempfile()))
  expect_warning(tree <- write.tree(biom$tree,       tempfile()))
  expect_warning(seqs <- write.fasta(biom$sequences, tempfile()))
  
  expect_warning(read.biom(file))
  expect_warning(read_biom(file))
  expect_warning(read.tree(tree))
  expect_warning(read.fasta(seqs))
  
  unlink(c(file, xlsx, tree, seqs))
  
  expect_warning(alpha.div(biom, rarefy = 'multi'))
  expect_warning(beta.div(biom, long = FALSE))
  expect_warning(beta.div(biom, long = TRUE))
  expect_warning(info(biom))
  expect_warning(nsamples(biom))
  expect_warning(ntaxa(biom))
  expect_warning(select(biom, nRandom = 3))
  expect_warning(subtree(biom$tree, 1:3))
  expect_warning(taxa.rollup(biom, long = FALSE))
  expect_warning(taxa.rollup(biom, long = TRUE))
  expect_warning(tips(biom$tree))
  expect_warning(unifrac(biom))
  expect_warning(as.percent(biom))
  expect_warning(depth(biom))
  expect_warning(depths_barplot(biom))
  expect_warning(has.phylogeny(biom))
  expect_warning(has.sequences(biom))
  expect_warning(is.rarefied(biom))
  expect_warning(repair(biom))
  expect_warning(sample_subset(biom))
  expect_warning(sample.sums(biom, long = FALSE))
  expect_warning(sample.sums(biom, long = TRUE))
  expect_warning(taxa_max(biom))
  expect_warning(taxa.means(biom))
  expect_warning(taxa.sums(biom))
  expect_warning(top.taxa(biom))
  expect_warning(top_taxa(biom))
  
  expect_warning(comments(biom))
  expect_warning(counts(biom))
  expect_warning(id(biom))
  expect_warning(metadata(biom))
  expect_warning(phylogeny(biom))
  expect_warning(sample.names(biom))
  expect_warning(sequences(biom))
  expect_warning(taxa.names(biom))
  expect_warning(taxa.ranks(biom))
  expect_warning(taxonomy(biom))
  
  expect_warning(comments(biom)     <- biom$comment)
  expect_warning(counts(biom)       <- biom$counts)
  expect_warning(id(biom)           <- biom$id)
  expect_warning(metadata(biom)     <- biom$metadata)
  expect_warning(phylogeny(biom)    <- biom$phylogeny)
  expect_warning(sample.names(biom) <- biom$samples)
  expect_warning(sequences(biom)    <- biom$sequences)
  expect_warning(taxa.names(biom)   <- biom$otus)
  expect_warning(taxa.ranks(biom)   <- biom$ranks)
  expect_warning(taxonomy(biom)     <- biom$taxonomy)
  
})