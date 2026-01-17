# Create a subtree by specifying tips to keep.

Create a subtree by specifying tips to keep.

## Usage

``` r
tree_subset(tree, tips, underscores = FALSE)
```

## Arguments

- tree:

  A phylo object, as returned from
  [`read_tree()`](https://cmmr.github.io/rbiom/reference/read_tree.md).

- tips:

  A character, numeric, or logical vector of tips to keep.

- underscores:

  When parsing the tree, should underscores be kept as is? By default
  they will be converted to spaces (unless the entire ID is quoted).
  Default `FALSE`

## Value

A `phylo` object for the subtree.

## See also

Other phylogeny:
[`read_tree()`](https://cmmr.github.io/rbiom/reference/read_tree.md)

## Examples

``` r
    library(rbiom)
    
    infile <- system.file("extdata", "newick.tre", package = "rbiom")
    tree <- read_tree(infile)
    tree
#> 
#> Phylogenetic tree with 20 tips and 19 internal nodes.
#> 
#> Tip labels:
#>   Pa5Bac29, AtlPorci, AciSp313, MxlBact8, MxlBacte, PseS1107, ...
#> 
#> Rooted; includes branch length(s).
    
    subtree <- tree_subset(tree, tips = head(tree$tip.label))
    subtree
#> 
#> Phylogenetic tree with 6 tips and 5 internal nodes.
#> 
#> Tip labels:
#>   Pa5Bac29, AtlPorci, AciSp313, MxlBact8, MxlBacte, PseS1107
#> 
#> Rooted; includes branch length(s).
```
