#' Display taxa abundances as a heatmap.
#'
#' @inherit documentation_heatmap
#' @inherit documentation_default
#' @inherit bdiv_heatmap sections
#'
#' @family taxa_abundance
#' @family visualization
#'
#' @param grid   Color palette name, or a list as expected [plot_heatmap()].
#'        Default: `"bilbao"`
#'
#' @param tracks   A character vector of metadata fields to display as tracks
#'        at the top of the plot. Or, a list as expected by the `tracks`
#'        argument of [plot_heatmap()]. Tracks describe samples, so they are
#'        always placed on top (`side = "left"` is ignored). Default: `NULL`
#'
#' @param ...   Additional arguments to pass on to ggplot2::theme().
#'        Arguments prefixed with `labs.` go to ggplot2::labs() instead,
#'        for example \code{labs.subtitle = "Plot subtitle"}.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     # Subset to 10 samples and rarefy them.
#'     hmp10 <- rarefy(hmp50[1:10])
#'
#'     taxa_heatmap(hmp10, rank = "Phylum", tracks = "Body Site")
#'
#'     taxa_heatmap(hmp10, rank = "Genus", tracks = c("sex", "bo"))
#'
#'     taxa_heatmap(hmp10, rank = "Phylum", tracks = list(
#'       'Sex'       = list(colors = c(m = "#0000FF", f = "violetred")),
#'       'Body Site' = list(colors = "muted", label = "Source") ))
#'
taxa_heatmap <- function (
    biom, rank = -1, taxa = 6, tracks = NULL, grid = "bilbao",
    other = FALSE, unc = "singly", lineage = FALSE,
    label = TRUE, label_size = NULL, rescale = "none", trees = TRUE,
    clust = "complete", dist = "euclidean",
    asp = 1, tree_height = 10, track_height = 10,
    legend = "right", title = TRUE, xlab.angle = "auto", ...) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('taxa_heatmap', params)
  schema_check('taxa_heatmap', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('taxa_heatmap', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {

    #________________________________________________________
    # Validate and restructure user's arguments.
    #________________________________________________________
    with(params, {

      validate_rank(max = Inf)
      validate_taxa()
      validate_var_choices('rescale', choices = CHOICES$rescale)
      validate_title(false_ok = TRUE)

      if (biom$n_samples < 1)
        cli_abort("At least one sample is needed for a heatmap.")

      if (!is_list(grid)) grid <- list(colors = grid)
      if (!is_scalar_character(grid$label))
        grid$label <- if (rescale != "none")      { "Rescaled Abundance"
        } else if (isTRUE(biom$relative))         { "Relative Abundance"
        } else if (is.null(biom$depth))           { "Unrarefied Reads"
        } else                                    { "Rarefied Reads" }

      if (length(clust) < 1) clust <- "complete"

      .fn <- 'taxa_heatmap()'
    })


    #________________________________________________________
    # Handle multiple ranks with one panel per rank.
    #________________________________________________________
    if (length(params$rank) > 1) {

      plots <- lapply(seq_along(params$rank), function (i) {
        pp <- rlang::env_clone(params)
        pp$rank <- params$rank[[i]]
        panel_params(pp, i)
        taxa_heatmap_build(pp)
      })

      p <- plot_wrap(plots, ncol = 1)

    } else {
      p <- taxa_heatmap_build(params)
    }

    set_cache_value(cache_file, p)
  }

  return (plot_prov(p, cmd, biom))
}



#________________________________________________________
# One taxa heatmap (one rank).
#________________________________________________________
taxa_heatmap_build <- function (params) {

  with(params, {

    if (isTRUE(title)) title <- paste(rank, "Abundance")


    #________________________________________________________
    # Matrix with samples in columns and taxa in rows.
    #________________________________________________________
    mtx <- taxa_matrix(
      biom    = biom,
      rank    = rank,
      taxa    = taxa,
      lineage = lineage,
      sparse  = FALSE,
      unc     = unc,
      other   = other )

    .mtx_code <- sprintf("mtx <- %s", attr(mtx, 'cmd', exact = TRUE))
  })


  #________________________________________________________
  # Convert `tracks` into a named list of lists. Tracks are per
  # sample, i.e. per column, so they all go on top.
  #________________________________________________________
  biom_tracks(params, side = "top")


  return (plot_heatmap_build(params))
}
