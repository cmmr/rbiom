% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_heatmap.r
\name{plot_heatmap}
\alias{plot_heatmap}
\title{Create a heatmap with tracks and dendrograms from any matrix.}
\usage{
plot_heatmap(
  mtx,
  grid = list(label = "Grid Value", colors = "imola"),
  tracks = NULL,
  label = TRUE,
  label_size = NULL,
  rescale = "none",
  trees = TRUE,
  clust = "complete",
  dist = "euclidean",
  asp = 1,
  tree_height = 10,
  track_height = 10,
  legend = "right",
  title = NULL,
  xlab.angle = "auto",
  ...
)
}
\arguments{
\item{mtx}{A numeric \code{matrix} with named rows and columns.}

\item{grid}{Color palette name, or a list with entries for \code{label},
\code{colors}, \code{range}, \code{bins}, \code{na.color}, and/or
\code{guide}. See the Track Definitions section for details.
Default: \code{list(label = "Grid Value", colors = "imola")}}

\item{tracks}{List of track definitions. See details below.
Default: \code{NULL}.}

\item{label}{Label the matrix rows and columns. You can supply a list
or logical vector of length two to control row labels and column
labels separately, for example
\code{label = c(rows = TRUE, cols = FALSE)}, or simply
\code{label = c(TRUE, FALSE)}. Other valid options are \code{"rows"},
\code{"cols"}, \code{"both"}, \code{"bottom"}, \code{"right"},
and \code{"none"}.
Default: \code{TRUE}}

\item{label_size}{The font size to use for the row and column labels. You
can supply a numeric vector of length two to control row label sizes
and column label sizes separately, for example
\code{c(rows = 20, cols = 8)}, or simply \code{c(20, 8)}.
Default: \code{NULL}, which computes:
\code{pmax(8, pmin(20, 100 / dim(mtx)))}}

\item{rescale}{Rescale rows or columns to all have a common min/max.
Options: \code{"none"}, \code{"rows"}, or \code{"cols"}.
Default: \code{"none"}}

\item{trees}{Draw a dendrogram for rows (left) and columns (top). You can
supply a list or logical vector of length two to control the row tree
and column tree separately, for example
\code{trees = c(rows = TRUE, cols = FALSE)},
or simply \code{trees = c(TRUE, FALSE)}.
Other valid options are \code{"rows"}, \code{"cols"}, \code{"both"},
\code{"left"}, \code{"top"}, and \code{"none"}.
Default: \code{TRUE}}

\item{clust}{Clustering algorithm for reordering the rows and columns by
similarity. You can supply a list or character vector of length two to
control the row and column clustering separately, for example
\code{clust = c(rows = "complete", cols = NA)}, or simply
\code{clust = c("complete", NA)}. Options are:
\describe{
\item{\code{FALSE} or \code{NA} - }{ Disable reordering. }
\item{An \code{hclust} class object}{ E.g. from \code{\link[stats:hclust]{stats::hclust()}}. }
\item{A method name - }{ \code{"ward.D"},
\code{"ward.D2"}, \code{"single"}, \code{"complete"},
\code{"average"}, \code{"mcquitty"}, \code{"median"}, or
\code{"centroid"}. }
}
Default: \code{"complete"}}

\item{dist}{Distance algorithm to use when reordering the rows and columns
by similarity. You can supply a list or character vector of length
two to control the row and column clustering separately, for example
\code{dist = c(rows = "euclidean", cols = "maximum")}, or simply
\code{dist = c("euclidean", "maximum")}. Options are:
\describe{
\item{A \code{dist} class object}{ E.g. from \code{\link[stats:dist]{stats::dist()}} or \code{\link[=bdiv_distmat]{bdiv_distmat()}}. }
\item{A method name - }{ \code{"euclidean"},
\code{"maximum"}, \code{"manhattan"}, \code{"canberra"},
\code{"binary"}, or \code{"minkowski"}. }
}
Default: \code{"euclidean"}}

\item{asp}{Aspect ratio (height/width) for entire grid.
Default: \code{1} (square)}

\item{tree_height, track_height}{The height of the dendrogram or annotation
tracks as a percentage of the overall grid size. Use a numeric vector
of length two to assign \code{c(top, left)} independently.
Default: \code{10} (10\% of the grid's height)}

\item{legend}{Where to place the legend. Options are: \code{"right"} or
\code{"bottom"}. Default: \code{"right"}}

\item{title}{Plot title. Default: \code{NULL}.}

\item{xlab.angle}{Angle of the labels at the bottom of the plot.
Options are \code{"auto"}, \code{'0'}, \code{'30'}, and \code{'90'}.
Default: \code{"auto"}.}

\item{...}{Additional arguments to pass on to ggplot2::theme().}
}
\value{
A \code{ggplot2} plot. The computed data points and ggplot
command are available as \verb{$data} and \verb{$code},
respectively.
}
\description{
Create a heatmap with tracks and dendrograms from any matrix.
}
\section{Track Definitions}{


One or more colored tracks can be placed on the left and/or top of the
heatmap grid to visualize associated metadata values.

\preformatted{## Categorical ----------------------------
cat_vals <- sample(c("Male", "Female"), 10, replace = TRUE)
tracks   <- list('Sex' = cat_vals)
tracks   <- list('Sex' = list(values = cat_vals, colors = "bright"))
tracks   <- list('Sex' = list(
  values = cat_vals, 
  colors = c('Male' = "blue", 'Female' = "red")) )

## Numeric --------------------------------
num_vals <- sample(25:40, 10, replace = TRUE)
tracks   <- list('Age' = num_vals)
tracks   <- list('Age' = list(values = num_vals, colors = "greens"))
tracks   <- list('Age' = list(values = num_vals, range = c(0,50)))
tracks   <- list('Age' = list(
  label  = "Age (Years)",
  values = num_vals, 
  colors = c("azure", "darkblue", "darkorchid") ))

## Multiple Tracks ------------------------
tracks <- list('Sex' = cat_vals, 'Age' = num_vals)
tracks <- list(
  list(label = "Sex", values = cat_vals, colors = "bright"),
  list(label = "Age", values = num_vals, colors = "greens") )
  
  
mtx           <- matrix(sample(1:50), ncol = 10)
dimnames(mtx) <- list(letters[1:5], LETTERS[1:10])
plot_heatmap(mtx = mtx, tracks = tracks)
}

The following entries in the track definitions are understood:

\describe{
\item{\code{values} - }{ The metadata values. When unnamed, order must match matrix. }
\item{\code{range} - }{ The c(min,max) to use for scale values. }
\item{\code{label} - }{ Label for this track. Defaults to the name of this list element. }
\item{\code{side} - }{ Options are \code{"top"} (default) or \code{"left"}. }
\item{\code{colors} - }{ A pre-defined palette name or custom set of colors to map to. }
\item{\code{na.color} - }{ The color to use for \code{NA} values. }
\item{\code{bins} - }{ Bin a gradient into this many bins/steps. }
\item{\code{guide} - }{ A list of arguments for guide_colorbar() or guide_legend(). }
}

All built-in color palettes are colorblind-friendly. See
\href{https://cmmr.github.io/rbiom/articles/aes.html#discrete-palettes}{Mapping Metadata to Aesthetics}
for images of the palettes.

Categorical palette names: \code{"okabe"}, \code{"carto"}, \code{"r4"},
\code{"polychrome"}, \code{"tol"}, \code{"bright"}, \code{"light"},
\code{"muted"}, \code{"vibrant"}, \code{"tableau"}, \code{"classic"},
\code{"alphabet"}, \code{"tableau20"}, \code{"kelly"}, and \code{"fishy"}.

Numeric palette names: \code{"reds"}, \code{"oranges"}, \code{"greens"},
\code{"purples"}, \code{"grays"}, \code{"acton"}, \code{"bamako"},
\code{"batlow"}, \code{"bilbao"}, \code{"buda"}, \code{"davos"},
\code{"devon"}, \code{"grayC"}, \code{"hawaii"}, \code{"imola"},
\code{"lajolla"}, \code{"lapaz"}, \code{"nuuk"}, \code{"oslo"},
\code{"tokyo"}, \code{"turku"}, \code{"bam"}, \code{"berlin"},
\code{"broc"}, \code{"cork"}, \code{"lisbon"}, \code{"roma"},
\code{"tofino"}, \code{"vanimo"}, and \code{"vik"}.
}

\examples{
    library(rbiom) 
    
    set.seed(123)
    mtx <- matrix(runif(5*8), nrow = 5, dimnames = list(LETTERS[1:5], letters[1:8]))
    
    plot_heatmap(mtx)
    plot_heatmap(mtx, grid="oranges")
    plot_heatmap(mtx, grid=list(colors = "oranges", label = "Some \%", bins = 5))
    
    tracks <- list(
      'Number' = sample(1:ncol(mtx)),
      'Person' = list(
        values = factor(sample(c("Alice", "Bob"), ncol(mtx), TRUE)),
        colors = c('Alice' = "purple", 'Bob' = "darkcyan") ),
      'State' = list(
        side   = "left",
        values = sample(c("TX", "OR", "WA"), nrow(mtx), TRUE),
        colors = "bright" )
    )
    
    plot_heatmap(mtx, tracks=tracks)
    
}
\seealso{
Other visualization: 
\code{\link{adiv_boxplot}()},
\code{\link{adiv_corrplot}()},
\code{\link{bdiv_boxplot}()},
\code{\link{bdiv_corrplot}()},
\code{\link{bdiv_heatmap}()},
\code{\link{bdiv_ord_plot}()},
\code{\link{rare_corrplot}()},
\code{\link{rare_multiplot}()},
\code{\link{rare_stacked}()},
\code{\link{stats_boxplot}()},
\code{\link{stats_corrplot}()},
\code{\link{taxa_boxplot}()},
\code{\link{taxa_corrplot}()},
\code{\link{taxa_heatmap}()},
\code{\link{taxa_stacked}()}
}
\concept{visualization}
