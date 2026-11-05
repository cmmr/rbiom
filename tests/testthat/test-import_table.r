test_that("import_table", {

  file_csv  <- tempfile(fileext = '.csv')
  file_tsv  <- tempfile(fileext = '.tsv')
  on.exit(unlink(c(file_csv, file_tsv)))

  write_delim_file(hmp5$metadata, file_csv, sep = ",")
  expect_equal_tibble(import_metadata(file_csv, sids = hmp5$samples), hmp5$metadata)

  write_delim_file(hmp5$taxonomy, file_tsv, sep = "\t")
  expect_equal_tibble(import_taxonomy(file_tsv, otus = hmp5$otus), hmp5$taxonomy)

  skip_on_cran()


  file_xlsx <- tempfile(fileext = '.xlsx')
  file_json <- tempfile(fileext = '.json')
  on.exit(unlink(c(file_xlsx, file_json)), add = TRUE)

  if (requireNamespace('openxlsx', quietly = TRUE) && requireNamespace('readxl', quietly = TRUE)) {
    openxlsx::write.xlsx(hmp5$metadata, file_xlsx, overwrite = TRUE)
    expect_equal_tibble(import_metadata(file_xlsx, sids = hmp5$samples), hmp5$metadata)
  }

  jsonlite::write_json(hmp5$metadata, file_json)
  expect_equal_tibble(import_metadata(file_json, sids = hmp5$samples), hmp5$metadata)

  # as literal strings of data instead of file names
  file_csv  <- readChar(file_csv,  file.size(file_csv))
  file_tsv  <- readChar(file_tsv,  file.size(file_tsv))
  file_json <- readChar(file_json, file.size(file_json))
  expect_equal_tibble(import_metadata(file_csv,  sids = hmp5$samples), hmp5$metadata)
  expect_equal_tibble(import_taxonomy(file_tsv,  otus = hmp5$otus),    hmp5$taxonomy)
  expect_equal_tibble(import_metadata(file_json, sids = hmp5$samples), hmp5$metadata)

  # Duplicated column names.
  expect_error(import_table(sub("BMI", "Age", file_csv, fixed = TRUE)))

  # Return tibble with just sample or OTU names.
  expect_equal_tibble(import_metadata(NULL, sids = hmp5$samples), hmp5$metadata[,1])
  expect_equal_tibble(import_taxonomy(NULL, otus = hmp5$otus),    hmp5$taxonomy[,1])

  # Sample / OTU names under a non-standard column name / position.
  x <- relocate(rename(hmp5$metadata, sample = .sample), sample, .after = BMI)
  y <- relocate(rename(hmp5$taxonomy, otu    = .otu),    otu,    .after = Kingdom)
  expect_equal_tibble(import_metadata(x, sids = hmp5$samples), hmp5$metadata)
  expect_equal_tibble(import_taxonomy(y, otus = hmp5$otus),    hmp5$taxonomy)

  # Sample / OTU names in rownames()
  x <- as.data.frame(hmp5$metadata[,-1]); rownames(x) <- hmp5$samples
  y <- as.data.frame(hmp5$taxonomy)[,-1]; rownames(y) <- hmp5$otus
  expect_equal_tibble(import_metadata(x, sids = hmp5$samples), hmp5$metadata)
  expect_equal_tibble(import_taxonomy(y, otus = hmp5$otus),    hmp5$taxonomy)

  # Missing or overabundant row ids.
  x <- as.data.frame(hmp5$metadata); rownames(x) <- rev(x[[1]])
  y <- as.data.frame(hmp5$taxonomy); rownames(y) <- rev(y[[1]])
  expect_error(import_metadata(x, sids = hmp5$samples))
  expect_error(import_taxonomy(y, otus = hmp5$otus))
  expect_error(import_metadata(hmp5$metadata[,-1], sids = hmp5$samples))
  expect_error(import_taxonomy(hmp5$taxonomy[,-1], otus = hmp5$otus))

  # Bad column name.
  x <- rename(hmp5$metadata, .BMI     = BMI)
  y <- rename(hmp5$taxonomy, .Kingdom = Kingdom)
  expect_error(import_metadata(x, sids = hmp5$samples))
  expect_error(import_taxonomy(y, otus = hmp5$otus))

  # No overlap in row ids.
  x <- mutate(hmp5$metadata, .sample = paste0('x', .sample))
  y <- mutate(hmp5$taxonomy,    .otu = paste0('x', .otu))
  expect_error(import_metadata(x, sids = hmp5$samples))
  expect_error(import_taxonomy(y, otus = hmp5$otus))

  # Incomplete overlap in row ids.
  expect_warning(import_metadata(hmp5$metadata, sids = hmp5$samples[-1]))
  expect_warning(import_taxonomy(hmp5$taxonomy, otus = hmp5$otus[-1]))
  expect_warning(import_metadata(hmp5$metadata, sids = c('x', hmp5$samples)))
  expect_warning(import_taxonomy(hmp5$taxonomy, otus = c('x', hmp5$otus)))

  # Duplicate row ids with different metadata.
  x <- hmp5$metadata; x$.sample[[1]] <- x$.sample[[2]]
  y <- hmp5$taxonomy; y$.otu[[1]]    <- y$.otu[[2]]
  expect_warning(import_metadata(x, sids = hmp5$samples[-1]))
  expect_warning(import_taxonomy(y, otus = hmp5$otus[-1]))

  # Duplicate row ids with identical metadata, quietly merged.
  x <- rbind(hmp5$metadata, hmp5$metadata)
  y <- rbind(hmp5$taxonomy, hmp5$taxonomy[1:5,])
  expect_equal_tibble(import_metadata(x,     sids = hmp5$samples), hmp5$metadata)
  expect_equal_tibble(import_taxonomy(y,     otus = hmp5$otus),    hmp5$taxonomy)
  expect_equal_tibble(import_metadata(x[,1], sids = hmp5$samples), hmp5$metadata[,1])
  expect_equal_tibble(import_taxonomy(y[,1], otus = hmp5$otus),    hmp5$taxonomy[,1])

})


test_that("text files: separators, typing, and sample column detection", {

  skip_on_cran()

  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f), add = TRUE)
  sids <- c("1.10", "2.20", "007")

  # Numeric-looking sample IDs stay character; other columns are typed.
  # A UTF-8 BOM, 'Body Site' header, a '#' in a value, and an apostrophe.
  writeLines(c(
    "﻿sample_id\tAge\tBody Site\tflag\twhen\tcode",
    "1.10\t22.0\tStool #1\tTRUE\t2020-01-02\t007",
    "2.20\t31\tO'Neil\tfalse\t2021-12-31\t008",
    "007\tNA\t\tNA\t\t009" ), f, useBytes = TRUE)

  x <- expect_silent(import_metadata(f, sids = sids))
  expect_identical(x$.sample, sids)
  expect_identical(names(x), c(".sample", "Age", "Body Site", "flag", "when", "code"))
  expect_identical(x$Age,  c(22, 31, NA))
  expect_identical(x$flag, c(TRUE, FALSE, NA))
  expect_identical(x$when, as.Date(c("2020-01-02", "2021-12-31", NA)))
  expect_identical(as.character(x$`Body Site`), c("Stool #1", "O'Neil", NA))
  expect_identical(as.character(x$code), c("007", "008", "009"))

  # Comma and semicolon separators; '#SampleID' header; QIIME 2 types line.
  writeLines(c("#SampleID,Age", "1.10,1", "2.20,2", "007,3"), f)
  expect_identical(import_metadata(f, sids = sids)$Age, c(1, 2, 3))

  writeLines(c("sample-id;Age", "#q2:types;numeric", "1.10;1", "2.20;2", "007;3"), f)
  expect_identical(import_metadata(f, sids = sids)$Age, c(1, 2, 3))

  # Quoted fields with embedded separators and doubled quotes.
  writeLines(c('id,Site', '"1.10","a, b"', '"2.20","say ""hi"""', '007,c'), f)
  x <- import_metadata(f, sids = sids)
  expect_identical(as.character(x$Site), c("a, b", 'say "hi"', "c"))

  # The exact-name column wins over an earlier overlapping column.
  df <- data.frame(other = c("1.10", "x", "y"), Sample = sids, v = 1:3, check.names = FALSE)
  expect_identical(import_metadata(df, sids = sids)$.sample, sids)

  # Ambiguous: two unnamed candidates.
  df <- data.frame(a = sids, b = rev(sids), check.names = FALSE)
  expect_error(import_metadata(df, sids = sids), "Candidate")

  # Row names in a header with one fewer field.
  writeLines(c("Age\tSex", "1.10\t5\tF", "2.20\t6\tM", "007\t7\tF"), f)
  x <- expect_silent(import_metadata(f, sids = sids))
  expect_identical(x$.sample, sids)
  expect_identical(x$Age, c(5, 6, 7))
})


test_that("QIIME 2 taxonomy with a single lineage column", {

  otus <- c("asv1", "asv2", "asv3")

  tx <- data.frame(
    'Feature ID' = otus,
    'Taxon'      = c("k__Bacteria; p__Firmicutes; c__Bacilli", "k__Bacteria; p__Bacteroidota", "Unassigned"),
    'Confidence' = c(0.99, 0.8, 0.1),
    check.names  = FALSE )

  x <- expect_silent(import_taxonomy(tx, otus = otus))
  expect_identical(names(x), c(".otu", "Level.1", "Level.2", "Level.3"))
  expect_identical(as.character(x$Level.1), c("k__Bacteria", "k__Bacteria", "Unassigned"))
  expect_identical(as.character(x$Level.2), c("p__Firmicutes", "p__Bacteroidota", NA))
  expect_identical(as.character(x$Level.3), c("c__Bacilli", NA, NA))

  f <- tempfile(fileext = ".tsv")
  on.exit(unlink(f), add = TRUE)
  write_delim_file(tx, f)
  expect_equal_tibble(import_taxonomy(f, otus = otus), x)

  # A lone non-lineage column keeps its name.
  tx <- data.frame(.otu = otus, Genus = c("a", "b", NA))
  expect_identical(names(import_taxonomy(tx, otus = otus)), c(".otu", "Genus"))

  # split_lineages keeps trailing empty ranks and NA rows.
  m <- split_lineages(c("a;b;", "a", NA, "", " a | b "))
  expect_identical(dim(m), c(5L, 3L))
  expect_identical(m[1,], c("a", "b", NA))
  expect_identical(m[3,], c(NA_character_, NA, NA))
  expect_identical(m[5,], c("a", "b", NA))
})
