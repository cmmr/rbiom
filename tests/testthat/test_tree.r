
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Working with Newick trees
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Phylogenetic Trees")

#----------------------------------------------------------
# Read/Write Newick File
#----------------------------------------------------------

original <- "((((((('Pa5 Bac_29':0.05208,Atl_Porci:0.07758):0.38379,(AciSp313:0.51812,(MxlBact8:0.01897,MxlBacte:0.05666):0.11191):0.12761):0.02218,PseS1107:0.68041):0.00872,((BurCep29:0.6125,HrbSpe65:0.13211):0.17282,((SteMa290:0.13815,SteSpe52:0.06709):0.05649,XnhBac23:0.14576):0.06986):0.04759):0.00947,LegYabu2:0.23091):0.11389,((RhbSph13:0.48056,BraSp452:0.22189):0.01105,((SppAlask:0.13478,SphSp203:0.1652):0.04423,SphSp482:0.26706):0.07172):0.08156):0.03387,(HelPy142:0.21121,(CmpHomi2:0.09831,CmpConc9:0.0047):0.21009):0.13918);"

f1 <- tempfile()
f2 <- tempfile()
on.exit(file.remove(c(f1,f2)))

writeChar(original, f1)

t1 <- rbiom::read.tree(file=f1)
t2 <- rbiom::read.tree(text=original)
newick <- rbiom::write.tree(tree=t1)
rbiom::write.tree(tree=t1, file=f2)

test_that("Newick File IO", {
  expect_equal(t1, t2)
  expect_equal(original, newick)
  expect_equal(original, readChar(f2, file.size(f2)) )
})


# t3 <- rbiom::subtree(tree=t1, tips=1:10)
# test_that("Subset Tree", {
#   expect_equal(t1, t2)
#   expect_equal(original, newick)
#   expect_equal(original, readChar(f2, file.size(f2)) )
# })




remove("t1", "t2", "original", "newick")






