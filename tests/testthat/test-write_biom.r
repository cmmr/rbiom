test_that(desc = "write_biom", code = {

  tfile <- tempfile(fileext = '.biom')

  expect_silent(write_biom(hmp5, tfile, format = "tab"))
  expect_silent(read_biom(tfile))
  expect_error(write_biom(hmp5, tfile), "overwrite")
  expect_silent(write_biom(hmp5, tfile, overwrite = TRUE))
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

  # All writers refuse to overwrite unless asked.
  expect_silent(write_counts(hmp5, tfile))
  expect_error(write_counts(  hmp5, tfile), "overwrite")
  expect_error(write_metadata(hmp5, tfile), "overwrite")
  expect_error(write_fasta(   hmp5, tfile), "overwrite")
  expect_error(write_tree(    hmp5, tfile), "overwrite")
  expect_silent(write_tree(   hmp5, tfile, overwrite = TRUE))
  expect_identical(read_tree(tfile)$tip.label, hmp5$tree$tip.label)
  unlink(tfile)

  # A failed write leaves nothing behind.
  expect_error(write_wrapper(tfile, function (con) stop("boom")), "boom")
  expect_false(file.exists(tfile))
  expect_length(list.files(dirname(tfile), pattern = basename(tfile), all.files = TRUE), 0)

  tdir <- tempfile()
  expect_silent(write_qiime2(hmp5, tdir, 'qiime2_'))
  expect_silent(write_mothur(hmp5, tdir, 'mothur_'))
  expect_setequal(
    object   = list.files(tdir),
    expected = c('mothur_counts.tsv', 'mothur_metadata.tsv', 'mothur_taxonomy.tsv',
                 'mothur_tree.nwk', 'mothur_seqs.fna', 'qiime2_counts.tsv', 'qiime2_metadata.tsv',
                 'qiime2_taxonomy.tsv', 'qiime2_tree.nwk', 'qiime2_seqs.fna'))
  expect_error(write_qiime2(hmp5, tdir, 'qiime2_'), "overwrite")

  # mothur taxonomy uses the OTU ids from the counts file.
  tax <- read.delim(file.path(tdir, 'mothur_taxonomy.tsv'), check.names = FALSE, quote = "", comment.char = "")
  cts <- read.delim(file.path(tdir, 'mothur_counts.tsv'),   check.names = FALSE, quote = "", comment.char = "")
  expect_identical(tax$OTU, cts$Representative_Sequence)
  expect_identical(tax$OTU, hmp5$otus)
  unlink(tdir, recursive = TRUE)

  # QIIME 2 column types for logical and Date fields.
  write_qiime2(hostile5, tdir)
  lines <- readLines(file.path(tdir, 'biom_metadata.tsv'), encoding = "UTF-8")
  expect_identical(strsplit(lines[[2]], "\t")[[1]], c("#q2:types", "numeric", "numeric", "categorical", "categorical", "categorical", "categorical", "numeric"))
  expect_true(grepl("2020-01-02", lines[[3]], fixed = TRUE))
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


  skip_if_not_installed('h5lite')
  expect_silent(write_biom(hmp5, tfile, format = "hdf5"))
  expect_silent(read_biom(tfile))
  expect_error(write_biom(hmp5, tfile, format = "hdf5"), "overwrite")
  expect_silent(write_biom(hmp5, tfile, format = "hdf5", overwrite = TRUE))
  unlink(tfile)

  # A failed HDF5 write leaves no partial file.
  bad <- hmp5$clone()
  bad$metadata$`a/b` <- 1:5
  expect_error(write_biom(bad, tfile, format = "hdf5"), "slash")
  expect_false(file.exists(tfile))
  expect_length(list.files(dirname(tfile), pattern = basename(tfile), all.files = TRUE), 0)
})


test_that("write_fasta", {

  seqs <- c(a = "ATCG", b = "GGCC")
  txt  <- write_fasta(seqs)
  expect_identical(txt, ">a\nATCG\n>b\nGGCC\n")
  expect_identical(read_fasta(txt), seqs)

  # NA sequences (e.g. from biom_merge) are skipped.
  expect_identical(write_fasta(c(a = "ATCG", b = NA)), ">a\nATCG\n")
  expect_error(write_fasta(c(a = NA_character_)))

  expect_type(write_fasta(hmp5), "character")
  expect_identical(read_fasta(write_fasta(hmp5)), hmp5$sequences)
})


test_that("write_tree quotes newick metacharacters", {

  tree <- read_tree("((A:1,B:2)n1:0.5,(C:1,D:3):0.25)root;")
  tree$tip.label  <- c("A b", "c,d", "e:f", "g;h")
  tree$node.label <- c("the root", "(n1)", "x_y")

  nwk <- write_tree(tree)
  expect_identical(nwk, "(('A b':1,'c,d':2)'(n1)':0.5,('e:f':1,'g;h':3)'x_y':0.25)'the root';")

  back <- read_tree(nwk)
  expect_identical(back$tip.label,  tree$tip.label)
  expect_identical(back$node.label, tree$node.label)
  expect_identical(back$edge,       tree$edge)
  expect_equal(back$edge.length,    tree$edge.length)

  # Embedded quotes are doubled; brackets are quoted.
  tree$tip.label[1] <- "it's [1]"
  expect_true(grepl("'it''s [1]'", write_tree(tree), fixed = TRUE))

  # Unicode, underscores, spaces, apostrophes, brackets, and every other
  # newick metacharacter survive a round trip, in tip and node labels.
  tree <- hmp5$tree
  hostile <- c("ünï code", "x_y", "a'b", "o__[Sapro]", "'quoted'", "c:d", "(e)", "f;g", "a,b", "[]")
  tree$tip.label[seq_along(hostile)] <- hostile
  tree$node.label <- rep_len(c("it's", "n[1]", "node 1"), tree$Nnode)
  back <- read_tree(write_tree(tree))
  expect_identical(back$tip.label,  tree$tip.label)
  expect_identical(back$node.label, tree$node.label)
  expect_identical(back$edge,       tree$edge)
  expect_equal(back$edge.length,    tree$edge.length)

  # No branch lengths; no node labels.
  tree <- read_tree("((A,B),(C,D));")
  expect_identical(write_tree(tree), "((A,B),(C,D));")

  # A 2000-tip ladder (caterpillar) tree writes quickly and exactly.
  n   <- 2000
  nwk <- paste0(paste(rep("(", n - 1), collapse = ""), "t1:1", paste0(",t", 2:n, ":1)", collapse = ""), ";")
  lad <- read_tree(nwk)
  t0  <- Sys.time()
  out <- write_tree(lad)
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 1)
  back <- read_tree(out)
  expect_identical(back$tip.label, lad$tip.label)
  expect_identical(back$edge,      lad$edge)
  expect_equal(back$edge.length,   lad$edge.length)

  skip_on_cran()

  # Compressed output
  f <- tempfile(fileext = ".nwk.gz")
  write_tree(hmp5, f)
  expect_identical(read_tree(f)$tip.label, hmp5$tree$tip.label)
  unlink(f)
})
