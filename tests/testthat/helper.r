
rare50 <- rarefy(hmp50)
# min50  <- as_rbiom(list(counts = hmp50$counts))

hmp5 <- hmp50[1:5]

# Make sure sample names can handle any characters, including unicode.
hmp5$samples <- c("~4L<  vl%", "\U1F60A]<?^@D'", "x\\T@;`\U26F9/", '{1"[~/q;$8>', "k_s \U2B50[|L:")

rare5 <- rarefy(hmp5)
min5  <- as_rbiom(list(counts = hmp5$counts))



expect_equal_rbiom <- function (a, b) {
  
  expect_identical(class(a), c('rbiom', 'R6'))
  expect_identical(class(b), c('rbiom', 'R6'))
  
  expect_identical(a$id,             b$id)
  expect_identical(a$comment,        b$comment)
  expect_identical(a$sequences,      b$sequences)
  expect_identical(a$tree$tip.label, b$tree$tip.label)
  
  expect_identical(lapply(a$metadata, as.character), lapply(b$metadata, as.character))
  expect_identical(lapply(a$taxonomy, as.character), lapply(b$taxonomy, as.character))
  
  expect_equal(a$depth,              b$depth)
  expect_equal(as.matrix(a$counts),  as.matrix(b$counts))
  expect_equal(a$tree$edge,          b$tree$edge)
  expect_equal(a$tree$edge.length,   b$tree$edge.length)
  
}
