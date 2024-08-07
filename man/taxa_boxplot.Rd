% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boxplot.r
\name{taxa_boxplot}
\alias{taxa_boxplot}
\title{Visualize BIOM data with boxplots.}
\usage{
taxa_boxplot(
  biom,
  x = NULL,
  rank = -1,
  layers = "x",
  taxa = 6,
  unc = "singly",
  other = FALSE,
  p.top = Inf,
  y.trans = "sqrt",
  stat.by = x,
  facet.by = NULL,
  colors = TRUE,
  shapes = TRUE,
  patterns = FALSE,
  flip = FALSE,
  stripe = NULL,
  ci = "ci",
  level = 0.95,
  p.adj = "fdr",
  outliers = NULL,
  xlab.angle = "auto",
  p.label = 0.05,
  trans = "none",
  caption = TRUE,
  ...
)
}
\arguments{
\item{biom}{An \link[=rbiom_objects]{rbiom object}, such as from \code{\link[=as_rbiom]{as_rbiom()}}.
Any value accepted by \code{\link[=as_rbiom]{as_rbiom()}} can also be given here.}

\item{x}{A categorical metadata column name to use for the x-axis. Or
\code{NULL}, which puts taxa along the x-axis. Default: \code{NULL}}

\item{rank}{What rank(s) of taxa to display. E.g. \code{"Phylum"},
\code{"Genus"}, \code{".otu"}, etc. An integer vector can also be
given, where \code{1} is the highest rank, \code{2} is the second
highest, \code{-1} is the lowest rank, \code{-2} is the second
lowest, and \code{0} is the OTU "rank". Run \code{biom$ranks} to
see all options for a given rbiom object. Default: \code{-1}.}

\item{layers}{One or more of
\code{c("bar", "box" ("x"), "violin", "dot", "strip", "crossbar", "errorbar", "linerange", "pointrange")}.
Single letter abbreviations are also accepted. For instance,
\code{c("box", "dot")} is equivalent to \code{c("x", "d")} and \code{"xd"}.
Default: \code{"x"}}

\item{taxa}{Which taxa to display. An integer value will show the top n
most abundant taxa. A value 0 <= n < 1 will show any taxa with that
mean abundance or greater (e.g. \code{0.1} implies >= 10\%). A
character vector of taxa names will show only those named taxa.
Default: \code{6}.}

\item{unc}{How to handle unclassified, uncultured, and similarly ambiguous
taxa names. Options are:
\describe{
\item{\code{"singly"} - }{ Replaces them with the OTU name. }
\item{\code{"grouped"} - }{ Replaces them with a higher rank's name. }
\item{\code{"drop"} - }{ Excludes them from the result. }
\item{\code{"asis"} - }{ To not check/modify any taxa names. }
}
Default: \code{"singly"} \cr\cr
Abbreviations are allowed.}

\item{other}{Sum all non-itemized taxa into an "Other" taxa. When
\code{FALSE}, only returns taxa matched by the \code{taxa}
argument. Specifying \code{TRUE} adds "Other" to the returned set.
A string can also be given to imply \code{TRUE}, but with that
value as the name to use instead of "Other".
Default: \code{FALSE}}

\item{p.top}{Only display taxa with the most significant differences in
abundance. If \code{p.top} is >= 1, then the \code{p.top} most
significant taxa are displayed. If \code{p.top} is less than one, all
taxa with an adjusted p-value <= \code{p.top} are displayed.
Recommended to be used in combination with the \code{taxa} parameter
to set a lower bound on the mean abundance of considered taxa.
Default: \code{Inf}}

\item{y.trans}{The transformation to apply to the y-axis. Visualizing
differences of both high- and low-abundance taxa is best done with
a non-linear axis. Options are:
\describe{
\item{\code{"sqrt"} - }{ square-root transformation }
\item{\code{"log1p"} - }{ log(y + 1) transformation }
\item{\code{NULL} - }{ no transformation }
}
These methods allow visualization of both high- and low-abundance
taxa simultaneously, without complaint about 'zero' count
observations. Default: \code{"sqrt"}}

\item{stat.by}{Dataset field with the statistical groups. Must be
categorical. Default: \code{NULL}}

\item{facet.by}{Dataset field(s) to use for faceting. Must be categorical.
Default: \code{NULL}}

\item{colors}{How to color the groups. Options are:
\describe{
\item{\code{TRUE} - }{ Automatically select colorblind-friendly colors. }
\item{\code{FALSE} or \code{NULL} - }{ Don't use colors. }
\item{a palette name - }{ Auto-select colors from this set. E.g. \code{"okabe"} }
\item{character vector - }{ Custom colors to use. E.g. \code{c("red", "#00FF00")} }
\item{named character vector - }{ Explicit mapping. E.g. \code{c(Male = "blue", Female = "red")} }
}
See "Aesthetics" section below for additional information.
Default: \code{TRUE}}

\item{shapes}{Shapes for each group.
Options are similar to \code{colors}'s: \code{TRUE}, \code{FALSE}, \code{NULL}, shape
names (typically integers 0 - 17), or a named vector mapping
groups to specific shape names.
See "Aesthetics" section below for additional information.
Default: \code{TRUE}}

\item{patterns}{Patterns for each group.
Options are similar to \code{colors}'s: \code{TRUE}, \code{FALSE}, \code{NULL}, pattern
names (\code{"brick"}, \code{"chevron"}, \code{"fish"}, \code{"grid"}, etc), or a named
vector mapping groups to specific pattern names.
See "Aesthetics" section below for additional information.
Default: \code{FALSE}}

\item{flip}{Transpose the axes, so that taxa are present as rows instead
of columns. Default: \code{FALSE}}

\item{stripe}{Shade every other x position. Default: \emph{same as flip}}

\item{ci}{How to calculate min/max of the \bold{crossbar},
\bold{errorbar}, \bold{linerange}, and \bold{pointrange} layers.
Options are: \code{"ci"} (confidence interval), \code{"range"},
\code{"sd"} (standard deviation), \code{"se"} (standard error), and
\code{"mad"} (median absolute deviation).
The center mark of \bold{crossbar} and \bold{pointrange} represents
the mean, except for \code{"mad"} in which case it represents the median.
Default: \code{"ci"}}

\item{level}{The confidence level for calculating a confidence interval.
Default: \code{0.95}}

\item{p.adj}{Method to use for multiple comparisons adjustment of
p-values. Run \code{p.adjust.methods} for a list of available
options. Default: \code{"fdr"}}

\item{outliers}{Show boxplot outliers? \code{TRUE} to always show.
\code{FALSE} to always hide. \code{NULL} to only hide them when
overlaying a dot or strip chart.  Default: \code{NULL}}

\item{xlab.angle}{Angle of the labels at the bottom of the plot.
Options are \code{"auto"}, \code{'0'}, \code{'30'}, and \code{'90'}.
Default: \code{"auto"}.}

\item{p.label}{Minimum adjusted p-value to display on the plot with a
bracket.
\describe{
\item{\code{p.label = 0.05} - }{ Show p-values that are <= 0.05. }
\item{\code{p.label = 0} - }{ Don't show any p-values on the plot. }
\item{\code{p.label = 1} - }{ Show all p-values on the plot. }
}
If a numeric vector with more than one value is
provided, they will be used as breaks for asterisk notation.
Default: \code{0.05}}

\item{trans}{Transformation to apply. Options are:
\code{c("none", "rank", "log", "log1p", "sqrt")}. \code{"rank"} is useful for
correcting for non-normally distributions before applying regression
statistics. Default: \code{"none"}}

\item{caption}{Add methodology caption beneath the plot.
Default: \code{TRUE}}

\item{...}{Additional parameters to pass along to ggplot2 functions.
Prefix a parameter name with a layer name to pass it to only that
layer. For instance, \code{d.size = 2} ensures only the points on the
\bold{dot} layer have their size set to \code{2}.}
}
\value{
A \code{ggplot2} plot. \cr The computed data points, ggplot2 command,
stats table, and stats table commands are available as \verb{$data},
\verb{$code}, \verb{$stats}, and \verb{$stats$code}, respectively.
}
\description{
Visualize BIOM data with boxplots.
}
\section{Aesthetics}{


All built-in color palettes are colorblind-friendly. The available
categorical palette names are: \code{"okabe"}, \code{"carto"}, \code{"r4"},
\code{"polychrome"}, \code{"tol"}, \code{"bright"}, \code{"light"},
\code{"muted"}, \code{"vibrant"}, \code{"tableau"}, \code{"classic"},
\code{"alphabet"}, \code{"tableau20"}, \code{"kelly"}, and \code{"fishy"}.

Patterns are added using the fillpattern R package. Options are \code{"brick"},
\code{"chevron"}, \code{"fish"}, \code{"grid"}, \code{"herringbone"}, \code{"hexagon"}, \code{"octagon"},
\code{"rain"}, \code{"saw"}, \code{"shingle"}, \code{"rshingle"}, \code{"stripe"}, and \code{"wave"},
optionally abbreviated and/or suffixed with modifiers. For example,
\code{"hex10_sm"} for the hexagon pattern rotated 10 degrees and shrunk by 2x.
See \code{\link[fillpattern:fill_pattern]{fillpattern::fill_pattern()}} for complete documentation of options.

Shapes can be given as per base R - numbers 0 through 17 for various shapes,
or the decimal value of an ascii character, e.g. a-z = 65:90; A-Z = 97:122 to use
letters instead of shapes on the plot. Character strings may used as well.
}

\examples{
    library(rbiom)
    
    biom <- rarefy(hmp50)
    
    taxa_boxplot(biom, stat.by = "Body Site", stripe = TRUE)
    taxa_boxplot(biom, layers = "bed", rank = c("Phylum", "Genus"), flip = TRUE)
    taxa_boxplot(
      biom    = subset(biom, `Body Site` \%in\% c('Saliva', 'Stool')), 
      taxa    = 3, 
      layers  = "ps", 
      stat.by = "Body Site",
      colors  = c('Saliva' = "blue", 'Stool' = "red") )
    
}
\seealso{
Other taxa_abundance: 
\code{\link{sample_sums}()},
\code{\link{taxa_clusters}()},
\code{\link{taxa_corrplot}()},
\code{\link{taxa_heatmap}()},
\code{\link{taxa_stacked}()},
\code{\link{taxa_stats}()},
\code{\link{taxa_sums}()},
\code{\link{taxa_table}()}

Other visualization: 
\code{\link{adiv_boxplot}()},
\code{\link{adiv_corrplot}()},
\code{\link{bdiv_boxplot}()},
\code{\link{bdiv_corrplot}()},
\code{\link{bdiv_heatmap}()},
\code{\link{bdiv_ord_plot}()},
\code{\link{plot_heatmap}()},
\code{\link{rare_corrplot}()},
\code{\link{rare_multiplot}()},
\code{\link{rare_stacked}()},
\code{\link{stats_boxplot}()},
\code{\link{stats_corrplot}()},
\code{\link{taxa_corrplot}()},
\code{\link{taxa_heatmap}()},
\code{\link{taxa_stacked}()}
}
\concept{taxa_abundance}
\concept{visualization}
