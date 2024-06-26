% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documentation.r
\name{documentation_cmp}
\alias{documentation_cmp}
\title{documentation_cmp}
\description{
documentation_cmp
}
\section{Metadata Comparisons}{


Prefix metadata fields with \code{==} or \code{!=} to limit comparisons to within or
between groups, respectively. For example, \code{stat.by = '==Sex'} will
run calculations only for intra-group comparisons, returning "Male" and
"Female", but NOT "Female vs Male". Similarly, setting
\code{stat.by = '!=Body Site'} will only show the inter-group comparisons, such
as "Saliva vs Stool", "Anterior nares vs Buccal mucosa", and so on.

The same effect can be achieved by using the \code{within} and \code{between}
parameters. \code{stat.by = '==Sex'} is equivalent to
\verb{stat.by = 'Sex', within = 'Sex'}.
}

\keyword{internal}
