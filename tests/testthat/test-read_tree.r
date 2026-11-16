test_that("read_tree", {

  tree <- "(t9,((t5,t2),(((t10,(t7,t4)),(t6,(t3,t1))),t8)));"
  tree <- expect_s3_class(read_tree(tree), 'phylo')
  expect_silent(tree_subset(hmp5$tree, 1:50))

  skip_on_cran()

  expect_error(tree_subset(tree, c('a', 'b')))
  expect_error(tree_subset(tree, c(TRUE, TRUE)))
  expect_silent(tree_subset(tree, !logical(10)))
  expect_error(tree_subset(tree, 1:20))
  expect_error(tree_subset(tree, -1:5))
  expect_error(tree_subset(tree, 2:5 / 2))
  expect_error(tree_subset(tree, tree))
  expect_error(tree_subset(tree, 1))

  # Files, compressed files, and a newick string with a leading newline.
  f <- tempfile(fileext = ".tre.gz")
  write_tree(hmp5, f)
  expect_identical(read_tree(f)$tip.label, hmp5$tree$tip.label)
  unlink(f)
  expect_identical(read_tree("\n (A:1,B:2); \n")$tip.label, c("A", "B"))
  expect_error(read_tree(1))
  expect_error(read_tree(tempfile(fileext = ".tre")))
})


test_that("read_tree follows the newick quoting rules", {

  # '' is an apostrophe; brackets inside quotes are literal; brackets
  # outside quotes are comments, which may themselves contain apostrophes.
  expect_identical(read_tree("('it''s':1,b:2);")$tip.label,               c("it's", "b"))
  expect_identical(read_tree("('x[1]':1,b[a comment]:2);")$tip.label,     c("x[1]", "b"))
  expect_identical(read_tree("(A[it's],B);")$tip.label,                   c("A", "B"))
  expect_identical(read_tree("((A,B)'node 1',C)'it''s';")$node.label,     c("it's", "node 1"))

  # Quoted labels keep underscores; unquoted ones follow `underscores`.
  expect_identical(read_tree("('a_b':1,c_d:2);")$tip.label,                       c("a_b", "c d"))
  expect_identical(read_tree("('a_b':1,c_d:2);", underscores = TRUE)$tip.label,   c("a_b", "c_d"))

  # Malformed input is an error, not a garbled tree.
  expect_error(read_tree("('abc:1,b:2);"), "quote")

  # The hostile OTU names in `hostile5` survive a write/read cycle.
  tree <- hostile5$tree
  expect_true(any(grepl("'", tree$tip.label, fixed = TRUE)))
  back <- read_tree(write_tree(tree))
  expect_identical(back$tip.label, tree$tip.label)
  expect_identical(back$edge,      tree$edge)
})


test_that("tree_subset", {

  tree <- read_tree("(((A:1,B:2)ab:0.5,C:3)abc:0.25,((D:1,E:1)de:2,F:4)def:0.75)root;")

  # Unsorted numeric tips give the same result as sorted ones.
  a <- tree_subset(tree, c(5, 1, 3))
  b <- tree_subset(tree, c(1, 3, 5))
  expect_identical(a, b)
  expect_identical(a$tip.label, c("A", "C", "E"))

  # Node labels are kept for surviving nodes; single-child nodes are
  # collapsed and their branch lengths summed.
  expect_identical(a$node.label, c("root", "abc"))
  expect_equal(a$edge.length[a$edge[,2] == 3], 0.75 + 2 + 1)   # root -> E
  expect_equal(a$edge.length[a$edge[,2] == 1], 0.5 + 1)        # abc -> A

  # Character and logical specifications agree.
  expect_identical(tree_subset(tree, c("E", "A", "C")), a)
  expect_identical(tree_subset(tree, c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)), a)

  # The new root is the MRCA of the kept tips.
  x <- tree_subset(tree, c("A", "B"))
  expect_identical(x$node.label, "ab")
  expect_equal(x$edge.length, c(1, 2))

  # Round trip through newick is exact.
  back <- read_tree(write_tree(a))
  expect_identical(back$edge, a$edge)
  expect_identical(back$tip.label, a$tip.label)
  expect_identical(back$node.label, a$node.label)
  expect_equal(back$edge.length, a$edge.length)

  # Keeping everything returns the tree unchanged.
  expect_identical(tree_subset(tree, 1:6), tree)

  # Trees without branch lengths or node labels.
  x <- tree_subset(read_tree("((A,B),(C,D));"), c("A", "D"))
  expect_null(x$edge.length)
  expect_null(x$node.label)
  expect_identical(x$tip.label, c("A", "D"))

  # Large subsets stay consistent with the full tree's tip order.
  x <- tree_subset(hmp50$tree, c(400, 2, 77, 150))
  expect_identical(x$tip.label, hmp50$tree$tip.label[c(2, 77, 150, 400)])
  expect_equal(sum(x$edge.length) > 0, TRUE)

  # Caching on and off give the same answer.
  on  <- rlang::with_options(tree_subset(hmp50$tree, 1:10), rbiom.cache_dir = tempfile())
  off <- rlang::with_options(tree_subset(hmp50$tree, 1:10), rbiom.cache_dir = "FALSE")
  expect_identical(on, off)
})


test_that("read_tree explains files that are not Newick", {

  f <- tempfile(fileext = ".tre")
  on.exit(unlink(f), add = TRUE)

  nexus <- c("#NEXUS", "begin trees;", "  tree t = (A:1,B:2);", "end;")
  writeLines(nexus, f)
  expect_error(read_tree(f), "NEXUS")
  expect_error(read_tree(f), "ape::read.nexus", fixed = TRUE)
  expect_error(read_tree(paste(nexus, collapse = "\n")), "NEXUS")

  writeLines("just some text", f)
  expect_error(read_tree(f), "not a Newick tree")

  writeLines("(A:1,B:x);", f)
  expect_error(read_tree(f), "Can't parse")
  expect_error(read_tree("(A:1,(B:2,C:3);"), "Can't parse")

  writeBin(as.raw(rep(c(0, 200, 7, 255, 31, 129), 50)), f)
  expect_error(read_tree(f), "binary file")

  expect_error(read_tree(test_path('inputs/biom.hdf5')), "not a Newick tree")
  expect_error(read_tree(test_path('inputs/viz.qzv')),   "visualization")

  # A file whose path has parentheses is still a file.
  d <- tempfile("tree (1)")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  writeLines("(A:1,B:2);", file.path(d, "t.nwk"))
  expect_identical(read_tree(file.path(d, "t.nwk"))$tip.label, c("A", "B"))
})


test_that("tree/OTU mismatches are counted and explained", {

  mtx <- matrix(1:8, nrow = 4, dimnames = list(paste0("otu_", 1:4), c("A", "B")))

  # Unquoted underscores become spaces by default.
  tree <- "((otu_1:1,otu_2:1):1,(OTU_3:1,otu_4:1):1,extra:1);"
  expect_error(as_rbiom(mtx, tree = tree), "OTUs missing from tree: 4 of 4")
  expect_error(as_rbiom(mtx, tree = tree), "not in the counts table: 5 of 5")
  expect_error(as_rbiom(mtx, tree = tree), "differ only by")
  expect_error(as_rbiom(mtx, tree = tree), "underscores")

  # Case still differs with underscores = TRUE.
  expect_error(as_rbiom(mtx, tree = tree, underscores = TRUE), "OTUs missing from tree: 1 of 4")
  expect_error(as_rbiom(mtx, tree = tree, underscores = TRUE), '"otu_3" vs "OTU_3"')

  # Unrelated names: no hint.
  err <- tryCatch(as_rbiom(mtx, tree = "((a:1,b:1):1,(c:1,d:1):1);"), error = conditionMessage)
  expect_match(err, "OTUs missing from tree: 4 of 4")
  expect_false(grepl("differ only by", err))

  # Extra tips alone are fine.
  tree <- "((otu_1:1,otu_2:1):1,(otu_3:1,otu_4:1):1,extra:1);"
  x <- expect_silent(as_rbiom(mtx, tree = tree, underscores = TRUE))
  expect_identical(sort(x$tree$tip.label), x$otus)
})
