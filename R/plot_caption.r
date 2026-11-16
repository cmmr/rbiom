
#________________________________________________________
# Sentences for plot captions (and, later, methods text) that state
# what was analysed and how. All are pure functions of their arguments
# and return a single string, or NULL when there is nothing to say.
#________________________________________________________


#________________________________________________________
# One sentence describing the abundance table that was analysed:
#
#   "Rarefied to 1,183 reads per sample."
#   "Relative abundances."
#   "Unrarefied counts (1,183 to 21,150 reads per sample)."
#
# `biom`    An rbiom object (or anything as_rbiom() accepts).
# `depth`   For functions that rarefy internally: the depth(s) the table
#           was rarefied to before the analysis. One value, or several
#           (e.g. the ten depths of a rarefaction curve). The sentence
#           then describes that, rather than the object passed in.
#________________________________________________________
input_table_caption <- function (biom, depth = NULL) {

  fmt <- function (x) {
    ifelse(
      test = abs(x - round(x)) < 1e-8,
      yes  = formatC(round(x), format = "d", big.mark = ","),
      no   = formatC(signif(x, 4), format = "fg", big.mark = ",") )
  }

  if (!is.null(depth)) {

    stopifnot(is.numeric(depth), length(depth) > 0, !anyNA(depth))
    depth <- sort(unique(depth))

    if (length(depth) == 1)
      return (sprintf("Rarefied to %s reads per sample.", fmt(depth)))

    return (sprintf(
      "Rarefied to %i depths, from %s to %s reads per sample.",
      length(depth), fmt(min(depth)), fmt(max(depth)) ))
  }


  biom <- as_rbiom(biom)
  sums <- as.vector(sample_sums(biom))
  sums <- sums[is.finite(sums)]

  if (length(sums) == 0) return ("No samples.")

  lo <- min(sums)
  hi <- max(sums)

  counts   <- biom$counts
  values   <- if (isS4(counts)) counts@x else as.vector(counts)   # sparse matrix: the non-zero values
  integers <- all(values %% 1 == 0)
  even     <- isTRUE(all.equal(lo, hi))

  if (isTRUE(biom$relative) || (even && !integers && isTRUE(all.equal(hi, 1))))
    return ("Relative abundances.")

  if (even && integers && hi > 1)
    return (sprintf("Rarefied to %s reads per sample.", fmt(hi)))

  if (even)
    return (sprintf("Normalized abundances (each sample sums to %s).", fmt(hi)))

  if (integers)
    return (sprintf("Unrarefied counts (%s to %s reads per sample).", fmt(lo), fmt(hi)))

  sprintf("Unnormalized abundances (sample sums from %s to %s).", fmt(lo), fmt(hi))
}



#________________________________________________________
# The parameters of beta diversity metrics that their names don't
# state: generalized UniFrac's `alpha`. NULL for the other metrics.
#________________________________________________________
bdiv_caption <- function (bdiv, alpha = 0.5) {

  ids <- vapply(as.character(bdiv), FUN.VALUE = '', function (id) {
    tryCatch(ecodive::match_metric(id, div = 'beta')$id, error = function (e) id)
  })

  if (!any(grepl("generalized|^gunifrac$", tolower(ids)))) return (NULL)

  sprintf("Generalized UniFrac with alpha = %s.", format(alpha))
}



#________________________________________________________
# Which sample pairs a beta diversity plot compares, e.g.
# "Comparing samples within Sex and between Body Site."
#________________________________________________________
within_between_caption <- function (within = NULL, between = NULL) {

  within  <- unique(within)
  between <- unique(between)

  parts <- c(
    if (length(within)  > 0) paste("within",  paste(within,  collapse = ", ")),
    if (length(between) > 0) paste("between", paste(between, collapse = ", ")) )

  if (length(parts) == 0) return (NULL)

  sprintf("Comparing samples %s.", paste(parts, collapse = " and "))
}



#________________________________________________________
# How a heat map's or bar chart's items were ordered, e.g.
# "'complete' clustering of euclidean distances". `dist` and `clust`
# are method names, or the dist / hclust objects used instead.
#________________________________________________________
clust_caption <- function (clust, dist, dist_label = "the given distances") {

  if (inherits(clust, 'hclust')) return ("the given clustering")

  sprintf(
    "'%s' clustering of %s", clust,
    if (inherits(dist, 'dist')) dist_label else paste(dist, "distances") )
}



#________________________________________________________
# Assemble a caption: each argument is one line (NULL and empty ones
# are dropped); the sentences within an argument share a line.
#________________________________________________________
caption_lines <- function (...) {

  lines <- lapply(list(...), function (x) {
    x <- as.character(unlist(x))
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) > 0) paste(x, collapse = " ")
  })

  lines <- unlist(lines)
  if (length(lines) == 0) return (NULL)

  paste(lines, collapse = "\n")
}
