
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Biom File Read/Write
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Newick File Reading")


#----------------------------------------------------------
# Read Newick File
#----------------------------------------------------------
curr <- rbiom::read.tree(file="inputs/newick.tre")
prev <- readRDS("outputs/newick.rds")

test_that("Newick from file", {
  expect_equal(curr, prev)
})

remove("curr", "prev")



