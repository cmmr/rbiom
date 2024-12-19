# 
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Working with Newick trees
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# context("Phylogenetic Trees")
# 
# #----------------------------------------------------------
# # Read/Write Newick File
# #----------------------------------------------------------
# 
# original <- "((((((('Pa5 Bac_29':0.05208,Atl_Porci:0.07758):0.38379,(AciSp313:0.51812,(MxlBact8:0.01897,MxlBacte:0.05666):0.11191):0.12761):0.02218,PseS1107:0.68041):0.00872,((BurCep29:0.6125,HrbSpe65:0.13211):0.17282,((SteMa290:0.13815,SteSpe52:0.06709):0.05649,XnhBac23:0.14576):0.06986):0.04759):0.00947,LegYabu2:0.23091):0.11389,((RhbSph13:0.48056,BraSp452:0.22189):0.01105,((SppAlask:0.13478,SphSp203:0.1652):0.04423,SphSp482:0.26706):0.07172):0.08156):0.03387,(HelPy142:0.21121,(CmpHomi2:0.09831,CmpConc9:0.0047):0.21009):0.13918);"
# 
# f1 <- tempfile()
# f2 <- tempfile()
# on.exit(file.remove(c(f1,f2)))
# 
# writeChar(original, f1)
# 
# t1 <- rbiom::read_tree(src=f1)
# t2 <- rbiom::read_tree(src=original)
# newick <- rbiom::write_tree(tree=t1)
# rbiom::write_tree(tree=t1, file=f2)
# 
# test_that("Newick File IO", {
#   expect_equal(t1, t2)
#   expect_equal(original, newick)
#   expect_equal(original, readChar(f2, file.size(f2)) )
# })
# 
# 
# t1  <- "((a,(b,c,d)),(e,f,g,(h,i)));"
# phy <- rbiom::read_tree(src=t1)
# t2  <- rbiom::write_tree(tree=phy)
# 
# test_that("Newick Short Names", {
#   expect_equal(t1, t2)
#   expect_equal(phy$tip.label, letters[1:9])
# })
# 
# 
# 
# phy <- rbiom::read_tree(src="((a,(b,c,d)),(e,f,g,(h,i)));")
# t1  <- rbiom::write_tree(subtree(phy, c('b', 'c', 'f', 'g', 'i')))
# t2  <- rbiom::write_tree(subtree(phy, c('a', 'b', 'c', 'f', 'g', 'h', 'i')))
# 
# phy <- rbiom::read_tree(src=original)
# t3  <- rbiom::write_tree(rbiom::subtree(phy, c('Atl Porci', 'SteMa290', 'HelPy142')))
# t4  <- rbiom::write_tree(rbiom::subtree(phy, c('AciSp313',  'MxlBact8', 'SphSp203')))
# 
# test_that("Newick Subtree", {
#   expect_equal(t1, "((b,c),(f,g,i));")
#   expect_equal(t2, "((a,(b,c)),(f,g,(h,i)));")
#   expect_equal(t3, "((Atl_Porci:0.49227,SteMa290:0.31209):0.15723,HelPy142:0.35039);")
#   expect_equal(t4, "((AciSp313:0.51812,MxlBact8:0.13088):0.28187,SphSp203:0.36271);")
# })
# 
# 
# 
# 
# remove("t1", "t2", "t3", "t4", "original", "newick", "phy")
