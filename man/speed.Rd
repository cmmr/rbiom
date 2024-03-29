% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/speed_ups.r
\name{speed}
\alias{speed}
\title{Speed Ups.}
\description{
When working with very large datasets, you can make use of these tips and
tricks to speed up operations on rbiom objects.
}
\section{Skip Cloning}{


Functions that modify rbiom objects, like \code{\link[=subset]{subset()}} and \code{\link[=rarefy]{rarefy()}}, will
automatically clone the object before modifying it. This is to make these
functions behave as most R users would expect - but at a performance trade
off.

Rather than:

\if{html}{\out{<div class="sourceCode r">}}\preformatted{biom <- subset(biom, ...)
biom <- rarefy(biom)
}\if{html}{\out{</div>}}

Modify \code{biom} in place like this:

\if{html}{\out{<div class="sourceCode r">}}\preformatted{subset(biom, clone = FALSE, ...)
rarefy(biom, clone = FALSE)

# Or:
biom$metadata \%<>\% subset(...)
biom$counts \%<>\% rarefy_cols()
}\if{html}{\out{</div>}}

\cr
}

\section{Drop Components}{

\subsection{Sequences}{

Reference sequences for OTUs will be imported along with the rest of your
dataset and stored in \verb{$sequences}. However, rbiom doesn't currently use
these sequences for anything (except writing them back out with
\code{\link[=write_biom]{write_biom()}} or \code{\link[=write_fasta]{write_fasta()}}).

You can delete them from your rbiom object with:

\if{html}{\out{<div class="sourceCode r">}}\preformatted{biom$sequences <- NULL
}\if{html}{\out{</div>}}
}

\subsection{Tree}{

The phylogenetic reference tree for OTUs is only used for calculating
UniFrac distances. If you aren't using UniFrac, the tree can be dropped
from the rbiom object with:

\if{html}{\out{<div class="sourceCode r">}}\preformatted{biom$tree <- NULL
}\if{html}{\out{</div>}}

Alternatively, you can store the tree separately from the rbiom object and
provide it to just the functions that use it. For example:

\if{html}{\out{<div class="sourceCode r">}}\preformatted{tree <- biom$tree
biom$tree <- NULL
dm <- bdiv_distmat(biom, 'unifrac', tree = tree)
}\if{html}{\out{</div>}}

\cr
}
}

\section{Increase Caching}{


Caching is enabled by default - up to 20 MB per R session.

For large datasets, increasing the cache size can help. The size is
specified in bytes by an R option or environment variable.

\if{html}{\out{<div class="sourceCode r">}}\preformatted{options(rbiom.cache_size=200 * 1024 ^ 2) # 200 MB
Sys.setenv(RBIOM_CACHE_SIZE=1024 ^ 3)    # 1 GB
}\if{html}{\out{</div>}}

\cr

You can also specify a cache directory where results can be preserved from
one R session to the next.

\if{html}{\out{<div class="sourceCode r">}}\preformatted{options(rbiom.cache_dir=tools::R_user_dir("rbiom", "cache"))
Sys.setenv(RBIOM_CACHE_DIR="~/rbiom_cache")
}\if{html}{\out{</div>}}

\cr

Other quick notes about caching:
\itemize{
\item Setting the cache directory to \code{"FALSE"} will disable caching.
\item R options will override environment variables.
\item The key hash algorithm can be set with \code{options(rbiom.cache_hash=rlang::hash)}.
\cr
}
}

\section{Summary Layers}{


The figure-generating functions allow you to display every data point.
However, when you have thousands of data points, rendering every single one
can be slow. Instead, set the \code{layers} parameter to use other options.

\if{html}{\out{<div class="sourceCode r">}}\preformatted{adiv_boxplot(biom, layers = "bl")  # bar, linerange
adiv_corrplot(biom, layers = "tc") # trend, confidence
bdiv_ord_plot(biom, layers = "e")  # ellipse
}\if{html}{\out{</div>}}

\cr
}

\keyword{internal}
