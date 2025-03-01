% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_tree.r
\name{tree_subset}
\alias{tree_subset}
\title{Create a subtree by specifying tips to keep.}
\usage{
tree_subset(tree, tips, underscores = FALSE)
}
\arguments{
\item{tree}{A phylo object, as returned from \code{\link[=read_tree]{read_tree()}}.}

\item{tips}{A character, numeric, or logical vector of tips to keep.}

\item{underscores}{When parsing the tree, should underscores be kept as
is? By default they will be converted to spaces (unless the entire ID
is quoted). Default \code{FALSE}}
}
\value{
A \code{phylo} object for the subtree.
}
\description{
Create a subtree by specifying tips to keep.
}
\examples{
    library(rbiom)
    
    infile <- system.file("extdata", "newick.tre", package = "rbiom")
    tree <- read_tree(infile)
    tree
    
    subtree <- tree_subset(tree, tips = head(tree$tip.label))
    subtree

}
\seealso{
Other phylogeny: 
\code{\link{read_tree}()}
}
\concept{phylogeny}
