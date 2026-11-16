test_that("read_biom", {

  skip_on_cran()

  f <- expect_silent(write_biom(hmp5, tempfile(), 'json'))
  g <- expect_silent(write_biom(min5, tempfile(), 'json'))
  expect_silent(read_biom(src = f))
  expect_silent(read_biom(src = g))
  unlink(c(f, g))

  f <- expect_silent(write_biom(hmp5, tempfile(fileext = '.gz'), 'tab'))
  g <- expect_silent(write_biom(min5, tempfile(fileext = '.gz'), 'tab'))
  expect_silent(read_biom(src = f))
  expect_silent(read_biom(src = g))
  unlink(c(f, g))

  expect_error(read_biom_hdf5(tempfile()))
  expect_error(parse_json_metadata(list(columns = list(list(metadata = c('a' = 1, 'a' = 1))))))

  json <- list(matrix_type = 'dense', rows = list(list(id='t1')), columns = list(list(id='s1')), data = list(1))
  expect_silent(parse_json_counts(json))
  expect_error(parse_json_counts(list(matrix_type = 'invalid')))
  expect_error(parse_json_counts(list(matrix_type = 'dense', data = list())))


  json <- list(generated_by = 'MicrobiomeDB', rows = list(list(id='Archaea;Euryarchaeota', metadata = list())))
  expect_identical(unname(parse_json_taxonomy(json)[1,]), c('Archaea', 'Euryarchaeota'))
  json <- list(rows = list(list(id='Unc01pdq', metadata = list(list("Bacteria", "__Fusobacteriota")))))
  expect_identical(unname(parse_json_taxonomy(json)[1,]), c('Bacteria', '__Fusobacteriota'))
  expect_null(parse_json_taxonomy(list(rows = list(list(id='t1')))))
  json <- list(rows = list(list(metadata=list(taxonomy='t1')), list(metadata=list(taxonomy='t2;t3'))))
  expect_silent(parse_json_taxonomy(json))

  # A JSON null inside the taxonomy array is NA, not dropped.
  json <- list(rows = list(list(id='x', metadata=list(taxonomy=list('a', NULL, 'c')))))
  expect_identical(unname(parse_json_taxonomy(json)[1,]), c('a', NA, 'c'))

  expect_null(parse_json_sequences(list(rows = list(list(metadata=list())))))


  skip_if_not_installed('h5lite')

  # https://github.com/biocore/biom-format/blob/master/examples/rich_sparse_otu_table_hdf5.biom
  x <- expect_silent(read_biom(test_path('inputs/biom.hdf5')))
  expect_identical(x$samples, paste0('Sample', 1:6))
  expect_identical(x$otus,    paste0('GG_OTU_', 1:5))
  expect_identical(x$ranks,   c('.otu', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
  expect_identical(x$fields,  c('.sample', 'BODY_SITE', 'BarcodeSequence', 'Description', 'LinkerPrimerSequence'))
  expect_match(x$date,         "^2014-07-29T16:16:36")
  expect_identical(x$generated_by, "example")

  f <- expect_silent(write_biom(hmp5, tempfile(), 'hdf5'))
  g <- expect_silent(write_biom(min5, tempfile(), 'hdf5'))
  expect_silent(read_biom(src = f))
  expect_silent(read_biom(src = g))
  unlink(c(f, g))
})


test_that("file metadata is kept", {

  f <- system.file("extdata", "hmp50.bz2", package = "rbiom")
  b <- read_biom(f)

  expect_identical(b$generated_by, "rbiom 2.0.4")
  expect_identical(b$date,         "2023-09-22T00:39:26Z")
  expect_identical(b$history[[1]], call('read_biom', f))

  # Built-in datasets replay from their names. (helper.r rebuilds `hmp50`
  # from the file, so look at the package's own copy.)
  expect_identical(rbiom::hmp50$history[[1]], as.name('hmp50'))

  # Properties given to read_biom() override the file.
  b <- read_biom(f, id = "custom", metadata = NULL)
  expect_identical(b$id, "custom")
  expect_identical(b$fields, ".sample")
})


test_that("round trips preserve hostile names and column types", {

  skip_on_cran()

  # JSON
  f <- tempfile(fileext = ".biom")
  expect_silent(write_biom(hostile5, f, 'json'))
  x <- expect_silent(read_biom(f))
  expect_equal_rbiom(x, hostile5)
  expect_same_types(x, hostile5)
  expect_identical(x$metadata$when, hostile5$metadata$when)
  expect_identical(x$metadata$flag, hostile5$metadata$flag)
  unlink(f)

  # JSON string
  json <- write_biom(hostile5, NULL, 'json')
  expect_true(grepl('"1e-07"', json, fixed = TRUE) || grepl('1e-07', json, fixed = TRUE))
  x <- expect_silent(read_biom(json))
  expect_equal_rbiom(x, hostile5)

  # Tabular: counts, OTU names, sample names and taxonomy only.
  f <- tempfile(fileext = ".tsv")
  expect_silent(write_biom(hostile5, f, 'tab'))
  x <- expect_silent(read_biom(f))
  expect_identical(x$otus,    hostile5$otus)
  expect_identical(x$samples, hostile5$samples)
  expect_equal(as.matrix(x$counts), as.matrix(hostile5$counts))
  expect_identical(lapply(x$taxonomy, as.character), lapply(hostile5$taxonomy, as.character))
  expect_false(any(grepl("NA", readLines(f, warn = FALSE)[-1], fixed = TRUE)))
  unlink(f)

  # Relative abundances are not rounded in JSON.
  rel <- biom_relativize(hostile5)
  f   <- tempfile(fileext = ".biom")
  write_biom(rel, f, 'json')
  x   <- read_biom(f)
  expect_equal(as.matrix(x$counts), as.matrix(rel$counts), tolerance = 1e-12)
  expect_true(x$relative)
  unlink(f)

  # HDF5
  skip_if_not_installed('h5lite')
  f <- tempfile(fileext = ".biom")
  expect_silent(write_biom(hostile5, f, 'hdf5'))
  x <- expect_silent(read_biom(f))
  expect_equal_rbiom(x, hostile5)
  expect_same_types(x, hostile5)
  expect_identical(x$date, hostile5$date)
  expect_identical(x$metadata$when, hostile5$metadata$when)
  expect_identical(x$metadata$flag, hostile5$metadata$flag)
  unlink(f)
})


test_that("tabular reader handles awkward files", {

  skip_on_cran()

  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f), add = TRUE)

  # No leading '#', a '#' inside an ID, apostrophes, backslashes, a trailing
  # tab, a comment line before the header, and a blank line.
  writeLines(c(
    "# Constructed by hand",
    "OTU_ID\tS1\tS2\tS3\ttaxonomy\t",
    "#1\t1\t2\t3\tBacteria; Firmicutes\t",
    "",
    "it's\t4\t5\t6\tBacteria; \t",
    "a\\b\t7\t8\t9\t; \t" ), f)

  x <- expect_silent(read_biom(f))
  expect_identical(x$otus,    c("#1", "it's", "a\\b"))
  expect_identical(x$samples, c("S1", "S2", "S3"))
  expect_equal(unname(as.matrix(x$counts)), matrix(1:9, 3, byrow = TRUE))
  expect_identical(as.character(x$taxonomy$Level.2), c("Firmicutes", NA, NA))
  expect_identical(as.character(x$taxonomy$Level.1), c("Bacteria", "Bacteria", NA))

  # Blank cells become 0 with a warning; non-numeric cells are an error
  # naming the cell; text columns are ignored with a warning.
  writeLines(c("#OTU ID\tS1\tS2\tDescription", "o1\t1\t\tfoo", "o2\t\t3\tbar"), f)
  expect_warning(expect_warning(x <- read_biom(f), "blank"), "Description")
  expect_equal(unname(as.matrix(x$counts)), matrix(c(1, 0, 0, 3), 2))

  writeLines(c("#OTU ID\tS1\tS2", "o1\t1\t2", "o2\tabc\t3"), f)
  expect_error(read_biom(f), "abc")

  # Comma separated with quoted fields.
  writeLines(c('"#OTU ID","S,1","S2"', '"o,1",1,2', 'o2,3,4'), f)
  x <- expect_silent(read_biom(f))
  expect_identical(x$samples, c("S,1", "S2"))
  expect_identical(x$otus,    c("o,1", "o2"))

  # JSON with a BOM, leading whitespace and a blank line is still JSON
  # (the content is sniffed when the extension isn't .txt/.tsv/.csv).
  json <- write_biom(min5, NULL, 'json')
  g <- tempfile(fileext = ".biom")
  on.exit(unlink(g), add = TRUE)
  writeLines(c("﻿", "  ", json), g)
  expect_identical(file_type(g), 'json')
  expect_equal(as.matrix(read_biom(g)$counts), as.matrix(min5$counts))
  expect_identical(file_type(test_path('inputs/other.zip')), 'zip')
  expect_identical(unclass(file_type(paste0("﻿\n  ", json))), 'json')
})


test_that("QIIME 2 artifacts (.qza) are unpacked", {

  # Fixtures hold a 4 OTU x 3 sample table, its taxonomy, and its tree.
  tree <- expect_silent(read_tree(test_path('inputs/tree.qza')))
  expect_identical(tree$tip.label, paste0("otu", 1:4))

  # Visualizations, unsupported types, and other zip files.
  expect_error(read_biom(test_path('inputs/viz.qzv')),   "visualization")
  expect_error(read_biom(test_path('inputs/alpha.qza')), "SampleData[AlphaDiversity]", fixed = TRUE)
  expect_error(read_biom(test_path('inputs/other.zip')), "readme.txt")

  # The wrong kind of artifact for the job.
  expect_error(read_biom(test_path('inputs/tree.qza')),    "Phylogeny[Rooted]", fixed = TRUE)
  expect_error(import_table(test_path('inputs/tree.qza')), "Phylogeny[Rooted]", fixed = TRUE)
  expect_error(read_biom(test_path('inputs/seqs.qza')),    "not a feature table")
  expect_error(read_fasta(test_path('inputs/tree.qza')),   "not sequences")

  # Sequences.
  seqs <- expect_silent(read_fasta(test_path('inputs/seqs.qza')))
  expect_identical(names(seqs), paste0("otu", 1:4))
  expect_identical(seqs[['otu2']], "ACGTTGCAAC")
  expect_identical(read_fasta(test_path('inputs/seqs.qza'), ids = "otu3"), c(otu3 = "GGGTTTAAAC"))

  # Temporary files are removed afterwards.
  before <- list.files(tempdir(), pattern = "^qza")
  read_tree(test_path('inputs/tree.qza'))
  try(read_biom(test_path('inputs/alpha.qza')), silent = TRUE)
  expect_identical(list.files(tempdir(), pattern = "^qza"), before)

  skip_if_not_installed('h5lite')

  expect_error(read_tree(test_path('inputs/table.qza')), "FeatureTable[Frequency]", fixed = TRUE)

  biom <- expect_silent(read_biom(test_path('inputs/table.qza')))
  expect_identical(biom$samples, paste0("S", 1:3))
  expect_identical(biom$otus,    paste0("otu", 1:4))
  expect_identical(unname(as.matrix(biom$counts)[2,]), c(2, 8, 0))

  biom <- expect_silent(as_rbiom(
    biom     = test_path('inputs/table.qza'),
    taxonomy = test_path('inputs/taxonomy.qza'),
    tree     = test_path('inputs/tree.qza'),
    sequences = test_path('inputs/seqs.qza') ))
  expect_identical(biom$sequences, seqs)

  # The sequences setter takes a fasta file, a fasta string, or a vector.
  f <- tempfile(fileext = ".fna")
  on.exit(unlink(f), add = TRUE)
  write_fasta(rev(seqs), f)
  biom$sequences <- NULL
  biom$sequences <- f
  expect_identical(biom$sequences, seqs)
  biom$sequences <- paste0(">", names(seqs), "\n", toupper(seqs), "\n", collapse = "")
  expect_identical(biom$sequences, seqs)
  biom$sequences <- rev(seqs)
  expect_identical(biom$sequences, seqs)
  expect_error(biom$sequences <- tempfile(fileext = ".fna"), "File not found")
  expect_error(biom$sequences <- test_path('inputs/tree.qza'), "not sequences")
  expect_error(biom$sequences <- 1:4, "fasta file")

  # A file path is replayed as a file path.
  x <- as_rbiom(test_path('inputs/table.qza'), sequences = f)
  expect_true(any(grepl("read_biom|as_rbiom", x$code)))
  expect_identical(x$sequences, seqs)
  expect_identical(names(biom$taxonomy), c(".otu", "Level.1", "Level.2", "Level.3"))
  expect_identical(as.character(biom$taxonomy$Level.3[[4]]), "c__Methanobacteria")
  expect_identical(biom$tree$tip.label, biom$otus)
})


test_that("counts tables: decimal commas, mothur shared, binary, and Excel", {

  f <- tempfile(fileext = ".tsv")
  on.exit(unlink(f), add = TRUE)

  # Decimal commas are reported, not silently dropped or converted.
  writeLines(c("#OTU ID\tA\tB", "o1\t1,5\t2,0", "o2\t3\t4,25"), f)
  expect_error(read_biom(f), "decimal commas")

  # mothur shared files have samples in rows.
  shared <- c("label\tGroup\tnumOtus\tOtu1\tOtu2\tOtu3", "0.03\tA\t3\t5\t0\t2", "0.03\tB\t3\t1\t4\t0")
  writeLines(shared, f)
  x <- expect_silent(read_biom(f))
  expect_identical(x$samples, c("A", "B"))
  expect_identical(x$otus,    c("Otu1", "Otu2", "Otu3"))
  expect_identical(unname(as.matrix(x$counts)[,'B']), c(1, 4, 0))

  writeLines(c(shared, "0.05\tA\t3\t5\t0\t2"), f)
  expect_error(read_biom(f), "2 labels")

  # Binary and UTF-16 files.
  writeBin(as.raw(rep(c(0, 200, 7, 255, 31, 129), 50)), f)
  expect_error(read_biom(f), "binary file")
  writeBin(c(as.raw(c(0xFF, 0xFE)), iconv("#OTU ID\tA\no1\t1\n", to = "UTF-16LE", toRaw = TRUE)[[1]]), f)
  expect_error(read_biom(f), "UTF-16")

  # Latin-1 text is read rather than failing on invalid UTF-8.
  writeBin(c(charToRaw("#OTU ID\tcaf"), as.raw(0xE9), charToRaw("\no1\t1\n")), f)
  x <- expect_silent(read_biom(f))
  expect_identical(x$samples, "café")

  # Excel workbooks, with or without the file extension.
  skip_if_not_installed('openxlsx')
  xl <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xl), add = TRUE)
  openxlsx::write.xlsx(data.frame(otu = c("a", "b"), S1 = 1:2), xl)
  expect_error(read_biom(xl), "Excel file")
  noext <- tempfile()
  on.exit(unlink(noext), add = TRUE)
  file.copy(xl, noext)
  expect_error(read_biom(noext), "Excel file")
})


test_that("HDF5 round trip keeps factor level order", {

  skip_if_not_installed('h5lite')

  biom <- hmp50$clone()
  md   <- biom$metadata
  md$Sex    <- factor(md$Sex, levels = c("Male", "Female", "Unused"))
  md$Dose   <- factor(rep(c("10", "2"), 25), levels = c("10", "2"))
  md$Sex[3] <- NA
  biom$metadata <- md

  expect_identical(levels(biom$metadata$Sex), c("Male", "Female"))

  f <- tempfile(fileext = ".h5")
  on.exit(unlink(f), add = TRUE)
  write_biom(biom, f, 'hdf5')

  x <- expect_silent(read_biom(f))
  expect_identical(levels(x$metadata$Sex),  c("Male", "Female"))
  expect_identical(levels(x$metadata$Dose), c("10", "2"))
  expect_identical(x$metadata$Sex,          biom$metadata$Sex)
  expect_identical(levels(x$metadata$`Body Site`), levels(biom$metadata$`Body Site`))
  expect_identical(names(x$metadata), names(biom$metadata))

  # The values are stored as plain strings, for other BIOM readers.
  expect_type(h5lite::h5_read(f, "sample/metadata/Sex"), "character")

  # The tree is a 1-element array; biom-format (Python) can't read a scalar.
  expect_identical(as.integer(h5lite::h5_dim(f, "observation/group-metadata/phylogeny")), 1L)
  expect_identical(x$tree, biom$tree)

  # Without the attribute (files from older versions): as before.
  h5lite::h5_delete(f, "sample/metadata/Sex",  attr = "rbiom_levels")
  h5lite::h5_delete(f, "sample/metadata/Dose", attr = "rbiom_levels")
  x <- expect_silent(read_biom(f))
  expect_identical(levels(x$metadata$Sex), c("Female", "Male"))
  expect_type(x$metadata$Dose, "double")

  # Levels that don't cover the values are ignored.
  write_biom(biom, f, 'hdf5', overwrite = TRUE)
  h5lite::h5_delete(f, "sample/metadata/Sex", attr = "rbiom_levels")
  h5lite::h5_write(c("Male", "Other"), f, "sample/metadata/Sex", attr = "rbiom_levels")
  x <- expect_silent(read_biom(f))
  expect_identical(levels(x$metadata$Sex), c("Female", "Male"))
})
