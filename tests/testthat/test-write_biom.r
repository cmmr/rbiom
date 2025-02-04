test_that(desc = "write_biom", code = {
  
  tfile <- tempfile(fileext = '.biom')
  
  expect_silent(write_biom(hmp5, tfile, format = "tab"))
  expect_silent(read_biom(tfile))
  expect_error(write_biom(hmp5, tfile))
  unlink(tfile)
  
  expect_silent(write_biom(hmp5, tfile, format = "json"))
  expect_silent(read_biom(tfile))
  expect_silent(readChar(tfile, nchars = file.size(tfile)))
  unlink(tfile)
  
  skip_on_cran()
  
  
  expect_silent(write_counts(  hmp5, tfile)); unlink(tfile)
  expect_silent(write_metadata(hmp5, tfile)); unlink(tfile)
  expect_silent(write_taxonomy(hmp5, tfile)); unlink(tfile)
  expect_silent(write_fasta(   hmp5, tfile)); unlink(tfile)
  expect_error(write_fasta(    min5, tfile))
  expect_error(write_tree(     min5, tfile))
  
  tdir <- tempfile()
  expect_silent(write_qiime2(hmp5, tdir, 'qiime2_'))
  expect_silent(write_mothur(hmp5, tdir, 'mothur_'))
  expect_setequal(
    object   = list.files(tdir), 
    expected = c('mothur_counts.tsv', 'mothur_metadata.tsv', 'mothur_taxonomy.tsv', 
                 'mothur_tree.nwk', 'mothur_seqs.fna', 'qiime2_counts.tsv', 'qiime2_metadata.tsv', 
                 'qiime2_taxonomy.tsv', 'qiime2_tree.nwk', 'qiime2_seqs.fna'))
  unlink(tdir, recursive = TRUE)
  
  gzfile <- tempfile(fileext = '.gz')
  bzfile <- tempfile(fileext = '.bz2')
  
  expect_silent(write_biom(min5, gzfile, format = "tab"))
  expect_silent(write_biom(min5, bzfile, format = "tab"))
  unlink(c(gzfile, bzfile))
  
  expect_silent(write_biom(min5, gzfile, format = "json"))
  expect_silent(write_biom(min5, bzfile, format = "json"))
  unlink(c(gzfile, bzfile))
  
  expect_silent(write_metadata(min5, gzfile))
  expect_silent(write_metadata(min5, bzfile))
  unlink(c(gzfile, bzfile))
  
  
  skip_if_not_installed('rhdf5')
  expect_silent(write_biom(hmp5, tfile, format = "hdf5"))
  expect_silent(read_biom(tfile))
  unlink(tfile)
  
})
