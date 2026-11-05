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
  expect_identical(unclass(file_type(paste0("﻿\n  ", json))), 'json')
})
