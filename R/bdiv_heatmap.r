#' Display beta diversities in an all vs all grid.
#'
#' @inherit documentation_heatmap
#' @inherit documentation_default
#'
#' @family beta_diversity
#' @family visualization
#'
#' @param grid   Color palette name, or a list with entries for \code{label},
#'        \code{colors}, \code{range}, \code{bins}, \code{na.color}, and/or
#'        \code{guide}. See the Track Definitions section for details.
#'        Default: `"devon"`
#'
#' @param tracks   A character vector of metadata fields to display as tracks
#'        at the top of the plot. Or, a list as expected by the `tracks`
#'        argument of [plot_heatmap()]. Default: `NULL`
#'
#' @param clust   Clustering algorithm for reordering the samples by
#'        similarity, applied to the `bdiv` distances. Options are:
#'        `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"`,
#'        `"mcquitty"`, `"median"`, `"centroid"`, or `FALSE` to disable
#'        reordering. Default: `"complete"`
#'
#' @param ...   Additional arguments to pass on to ggplot2::theme().
#'        Arguments prefixed with `labs.` go to ggplot2::labs() instead,
#'        for example \code{labs.subtitle = "Plot subtitle"}.
#'
#'
#' @section Annotation Tracks:
#'
#' Metadata can be displayed as colored tracks above the heatmap. Common use
#' cases are provided below, with more thorough documentation available at
#' https://cmmr.github.io/rbiom .
#'
#' \preformatted{## Categorical ----------------------------
#' tracks = "Body Site"
#' tracks = list('Body Site' = "bright")
#' tracks = list('Body Site' = c('Stool' = "blue", 'Saliva' = "green"))
#'
#' ## Numeric --------------------------------
#' tracks = "Age"
#' tracks = list('Age' = "reds")
#'
#' ## Multiple Tracks ------------------------
#' tracks = c("Body Site", "Age")
#' tracks = list('Body Site' = "bright", 'Age' = "reds")
#' tracks = list(
#'   'Body Site' = c('Stool' = "blue", 'Saliva' = "green"),
#'   'Age'       = list('colors' = "reds") )
#' }
#'
#' The following entries in the track definitions are understood:
#'
#' \describe{
#'   \item{\code{colors} - }{ A pre-defined palette name or custom set of colors to map to. }
#'   \item{\code{range} - }{ The c(min,max) to use for scale values. }
#'   \item{\code{label} - }{ Label for this track. Defaults to the name of this list element. }
#'   \item{\code{side} - }{ Options are \code{"top"} (default) or \code{"left"}.
#'         `taxa_heatmap()` places all tracks on top. }
#'   \item{\code{na.color} - }{ The color to use for \code{NA} values. }
#'   \item{\code{bins} - }{ Bin a gradient into this many bins/steps. }
#'   \item{\code{guide} - }{ A list of arguments for guide_colorbar() or guide_legend(). }
#' }
#'
#' The default categorical colors are colorblind-friendly for up to 12
#' groups: Okabe and Ito's palette (\code{"okabe"}) for up to 8, and Paul
#' Tol's (\code{"tol"}) for 9-12.
#'
#' Categorical palette names: \code{"okabe"}, \code{"carto"}, \code{"r4"},
#' \code{"polychrome"}, \code{"tol"}, \code{"bright"}, \code{"light"},
#' \code{"muted"}, \code{"vibrant"}, \code{"tableau"}, \code{"classic"},
#' \code{"alphabet"}, \code{"tableau20"}, \code{"kelly"}, and \code{"fishy"}.
#'
#' Numeric palette names: \code{"reds"}, \code{"oranges"}, \code{"greens"},
#' \code{"purples"}, \code{"grays"}, \code{"acton"}, \code{"bamako"},
#' \code{"batlow"}, \code{"bilbao"}, \code{"buda"}, \code{"davos"},
#' \code{"devon"}, \code{"hawaii"}, \code{"imola"},
#' \code{"lajolla"}, \code{"lapaz"}, \code{"nuuk"}, \code{"oslo"},
#' \code{"tokyo"}, \code{"turku"}, \code{"bam"}, \code{"berlin"},
#' \code{"broc"}, \code{"cork"}, \code{"lisbon"}, \code{"roma"},
#' \code{"tofino"}, \code{"vanimo"}, and \code{"vik"}.
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     # Subset to 10 samples and rarefy them.
#'     hmp10 <- rarefy(hmp50[1:10])
#'
#'     bdiv_heatmap(hmp10, tracks=c("Body Site", "Age"))
#'
#'     bdiv_heatmap(hmp10, bdiv=c("u_unifrac", "n_unifrac"), tracks="sex")

bdiv_heatmap <- function (
    biom, bdiv = "bray", tree = NULL, tracks = NULL,
    grid = "devon", label = TRUE, label_size = NULL, rescale = "none",
    clust = "complete", trees = TRUE, asp = 1, tree_height = 10,
    track_height = 10, legend = "right", title = TRUE, xlab.angle = "auto",
    alpha = 0.5, cpus = n_cpus(), ...) {

  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)

  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('bdiv_heatmap', cmd_params(params))
  schema_check('bdiv_heatmap', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('bdiv_heatmap', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {

    #________________________________________________________
    # Validate and restructure user's arguments.
    #________________________________________________________
    with(params, {

      validate_bdiv(multiple = TRUE)

      if (biom$n_samples < 2)
        cli_abort("At least two samples are needed for a bdiv heatmap.")

      validate_title()

      if (!is_list(grid)) grid <- list(label = "Dissimilarity", colors = grid)

      .fn <- 'bdiv_heatmap()'
    })


    #________________________________________________________
    # Handle multiple metrics with one panel per metric.
    #________________________________________________________
    if (length(params$bdiv) > 1) {

      plots <- lapply(seq_along(params$bdiv), function (i) {
        pp <- rlang::env_clone(params)
        pp$bdiv <- params$bdiv[[i]]
        panel_params(pp, i)
        bdiv_heatmap_build(pp)
      })

      p <- plot_wrap(plots)

    } else {
      p <- bdiv_heatmap_build(params)
    }

    set_cache_value(cache_file, p)
  }

  return (plot_prov(p, cmd, biom))
}



#________________________________________________________
# One beta diversity heatmap (one metric).
#________________________________________________________
bdiv_heatmap_build <- function (params) {

  with(params, {

    if (isTRUE(title)) title <- ecodive::match_metric(bdiv, div = 'beta')$name


    #________________________________________________________
    # Samples are clustered on the plotted distances.
    #________________________________________________________
    dist <- bdiv_distmat(
      biom  = biom,
      bdiv  = bdiv,
      tree  = tree,
      alpha = alpha,
      cpus  = cpus )

    validate_finite_dist(dist, bdiv, "a heat map")

    mtx <- as.matrix(dist)

    .mtx_code <- c(
      sprintf("dm  <- %s", attr(dist, 'cmd', exact = TRUE)),
      "mtx <- as.matrix(dm)" )

    .dist_code <- "dm"
  })


  #________________________________________________________
  # Convert `tracks` into a named list of lists.
  #________________________________________________________
  biom_tracks(params)


  return (plot_heatmap_build(params))
}
