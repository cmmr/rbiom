#' Display taxa abundances as a stacked bar graph.
#'
#' @inherit documentation_default
#' @inherit documentation_plot_return return
#'
#' @family taxa_abundance
#' @family visualization
#'
#' @param colors,patterns   A character vector of colors or patterns to use in
#'        the graph. A named character vector can be used to map taxon names to
#'        specific colors or patterns. Set to `TRUE` to auto-select colors
#'        or patterns, or to `FALSE` to disable per-taxa colors or patterns.
#'        Patterns require the fillpattern package.
#'        Default: \code{colors=TRUE, patterns=FALSE}.
#'
#' @param label.by,order.by   What metadata column to use for labeling and/or
#'        sorting the samples across the x-axis. Set \code{label.by='.sample'}
#'        to display sample names. When \code{order.by=NULL}, samples are
#'        arranged based on \code{dist} and \code{clust}, below.
#'        Default: \code{label.by=NULL, order.by=NULL}.
#'
#' @param dist,clust   Distance ([stats::dist()]) and clustering
#'        ([stats::hclust()]) methods to use for automatically arranging
#'        samples along the x-axis to put samples with similar composition near
#'        one another. Default: \code{dist="euclidean", clust="complete"}.
#'
#' @param xlab.angle   Angle of the sample labels at the bottom of the plot.
#'        Options are `"auto"`, `"0"`, `"30"`, and `"90"` (the numbers 0, 30,
#'        and 90 are accepted too). Default: `"90"`
#'
#' @param ...   Parameters for underlying functions. Prefixing parameter names
#'        with a layer name ensures that a particular parameter is passed to,
#'        and only to, that layer. An argument that no layer accepts is an
#'        error.
#'
#' @details
#' If `biom` is rarefied (or holds relative abundances), then relative
#' abundance will be shown on the y-axis. Otherwise, raw abundance will be
#' displayed. Samples with `NA` values in `label.by`, `order.by`, or
#' `facet.by` are dropped.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(hmp50)
#'
#'     taxa_stacked(biom, rank="Phylum")
#'
#'     taxa_stacked(biom, rank = "genus", facet.by = "body site")
#'
taxa_stacked <- function (
    biom, rank = -1, taxa = 6, colors = TRUE, patterns = FALSE,
    label.by = NULL, order.by = NULL, facet.by = NULL,
    dist = "euclidean", clust = "complete", other = TRUE,
    unc = "singly", lineage = FALSE, xlab.angle = "90", ...) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('taxa_stacked', params)
  schema_check('taxa_stacked', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('taxa_stacked', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {

    #________________________________________________________
    # Sanity Checks
    #________________________________________________________
    with(params, {

      validate_rank(max = Inf)
      validate_dist()
      validate_clust()

      validate_biom_field('label.by', null_ok = TRUE)
      validate_biom_field('order.by', null_ok = TRUE, max = Inf)
      validate_biom_field('facet.by', null_ok = TRUE, max = Inf, col_type = "cat")

      xlab.angle <- as.character(xlab.angle)
      validate_var_choices('xlab.angle', choices = CHOICES$xlab.angle)


      #________________________________________________________
      # Drop samples with NA values in the metadata fields used.
      #________________________________________________________
      .fields <- unique(c(label.by, order.by, facet.by))
      if (length(.fields) > 0 && !all(stats::complete.cases(biom$metadata[.fields])))
        biom <- na.omit(biom, fields = .fields)
      remove(".fields")

      if (biom$n_samples < 1)
        cli_abort("At least one sample is needed for a stacked bar plot.")

      .fn <- 'taxa_stacked()'
    })


    #________________________________________________________
    # Handle multiple ranks with one panel per rank.
    #________________________________________________________
    if (length(params$rank) > 1) {

      plots <- lapply(seq_along(params$rank), function (i) {
        pp <- rlang::env_clone(params)
        pp$rank <- params$rank[[i]]
        pp$.dots$labs.title %<>% if.null(params$rank[[i]])
        panel_params(pp, i)
        taxa_stacked_build(pp)
      })

      p <- plot_wrap(plots, ncol = 1)

    } else {
      p <- taxa_stacked_build(params)
    }

    set_cache_value(cache_file, p)
  }

  return (plot_prov(p, cmd, params$biom))
}



#________________________________________________________
# One stacked bar plot (one rank).
#________________________________________________________
taxa_stacked_build <- function (params) {

  params$.dv <- params$.data_var %||% "data"


  #________________________________________________________
  # Compute taxa abundance values.
  #________________________________________________________
  with(params, {

    # The matrix will be needed for similarity clustering.
    #________________________________________________________
    mtx <- taxa_matrix(
      biom    = biom,
      rank    = rank,
      taxa    = taxa,
      lineage = lineage,
      sparse  = FALSE,
      unc     = unc,
      other   = other )


    # Pivot to long-form and add metadata.
    #________________________________________________________
    .ggdata <- tibble(
        '.rank'      = rank,
        '.sample'    = colnames(mtx)[col(mtx)],
        '.taxa'      = factor(rownames(mtx)[row(mtx)], levels = rownames(mtx)),
        '.abundance' = as.numeric(mtx) ) %>%
      dplyr::left_join(biom$metadata, by = '.sample')

    .data_code <- c(
      sprintf("mtx <- %s", attr(mtx, 'cmd', exact = TRUE)),
      "",
      sprintf("%s <- tibble::tibble(", .dv),
      sprintf("    .rank      = %s,", double_quote(rank)),
      "    .sample    = colnames(mtx)[col(mtx)],",
      "    .taxa      = factor(rownames(mtx)[row(mtx)], levels = rownames(mtx)),",
      "    .abundance = as.numeric(mtx) ) %>%",
      "  dplyr::left_join(biom$metadata, by = \".sample\")" )


    # Factor levels on .sample control x-axis ordering.
    #________________________________________________________
    .hc   <- stats::hclust(d = stats::dist(x = t(mtx), method = dist), method = clust)
    .lvls <- .hc$labels[.hc$order]

    .data_code <- c(
      .data_code,
      "",
      "# Arrange similar samples near one another.",
      sprintf("hc   <- stats::hclust(stats::dist(t(mtx), method = %s), method = %s)", double_quote(dist), double_quote(clust)),
      "lvls <- hc$labels[hc$order]" )

    if (!is.null(order.by)) {

      .lvls <- biom$samples[do.call(order, c(unname(as.list(biom$metadata[order.by])), list(match(biom$samples, .lvls))))]

      .data_code <- c(
        .data_code,
        sprintf("lvls <- biom$samples[do.call(order, c(unname(as.list(biom$metadata[%s])), list(match(biom$samples, lvls))))]", deparse_vec(order.by)) )
    }

    .ggdata[['.sample']] <- factor(.ggdata[['.sample']], levels = .lvls)
    .data_code <- c(.data_code, sprintf("%s[[\".sample\"]] <- factor(%s[[\".sample\"]], levels = lvls)", .dv, .dv))

    remove(".hc", ".lvls")


    .xcol  <- '.sample'
    .ycol  <- '.abundance'
    .xmode <- 'factor'

    stat.by    <- ".taxa"
    color.by   <- ".taxa"
    pattern.by <- ".taxa"
  })


  #________________________________________________________
  # Automatic colors: the "Other" taxa group is light gray.
  #________________________________________________________
  if (isTRUE(params$colors) && !(isFALSE(params$other) || is.null(params$other))) {
    lvls  <- levels(params$.ggdata[['.taxa']])
    other <- if (is_scalar_character(params$other)) params$other else "Other"
    if (eq(tail(lvls, 1), other))
      params$colors <- c(get_n_colors(length(lvls) - 1), "#BBBBBB")
    remove("lvls", "other")
  }



  #________________________________________________________
  # Provenance-tracked ggplot2 functions.
  #________________________________________________________
  .label_number  <- P('scales::label_number')
  .cut_si        <- P('scales::cut_si')



  #________________________________________________________
  # Initialize the `layers` object.
  #________________________________________________________
  init_layers(params, do_init = "stack")



  #________________________________________________________
  # Display a label on the x-axis other than sample name.
  #________________________________________________________
  if (is.null(params$label.by)) {

    set_layer(params, 'xaxis', labels = NULL, breaks = NULL)

  } else if (!eq(params$label.by, ".sample")) {

    # Convert `.sample` to "<int>-<label>" format for easy parsing.
    with(params, {

      levels(.ggdata[['.sample']]) <- local({
        lvls <- levels(.ggdata[['.sample']])
        labs <- pull(biom, label.by)
        paste(seq_along(lvls), labs[lvls], sep = "-")
      })

      .data_code <- c(
        .data_code,
        sprintf("levels(%s[[\".sample\"]]) <- paste(seq_along(lvls), pull(biom, %s)[lvls], sep = \"-\")", .dv, double_quote(label.by)) )
    })

    set_layer(params, 'xaxis', labels = ~ sub("^\\d+\\-", "", .))
  }


  #________________________________________________________
  # Caption about labeling and ordering
  #________________________________________________________
  set_layer(
    params = params,
    layer  = 'labs',
    'x'       = "Sample",
    'caption' = with(params, {

      if ((is_null(label.by) || eq(label.by, '.sample')) && is_null(order.by)) {
        sprintf("Arranged by %s distance and '%s' clustering.", dist, clust)

      } else if ((is_null(label.by) || eq(label.by, '.sample'))) {
        sprintf("Ordered by %s.", paste(order.by, collapse = ", "))

      } else if (is_null(order.by)) {
        sprintf("Labeled by %s and arranged by\n%s distance and\n'%s' clustering.", label.by, dist, clust)

      } else if (eq(label.by, order.by)) {
        sprintf("Labeled and ordered by %s.", label.by)

      } else {
        sprintf("Labeled by %s and ordered by %s.", label.by, paste(order.by, collapse = ", "))
      }
    }))


  #________________________________________________________
  # If rarefied, display relative abundance in %
  #________________________________________________________
  depth <- unique(sample_sums(params$biom))
  if (length(depth) == 1) {

    if (depth != 1) {

      params$.ggdata[[params$.ycol]] <- params$.ggdata[[params$.ycol]] / depth

      params$.data_code <- c(
        params$.data_code,
        "",
        "# Counts to relative abundances.",
        sprintf("%s[[\".abundance\"]] <- %s[[\".abundance\"]] / %s", params$.dv, params$.dv, deparse_vec(depth)) )
    }

    set_layer(
      params = params,
      layer  = 'labs',
      'y' = "Relative Abundance" )

    set_layer(
      params = params,
      layer  = 'yaxis',
      'expand' = c(0,0),
      'limits' = c(0,1),
      'labels' = as.cmd(scales::percent) )


  } else {

    set_layer(
      params = params,
      layer  = 'labs',
      'y' = "Raw Abundance" )

    set_layer(
      params = params,
      layer  = 'yaxis',
      'expand' = c(0, 0, 0.02, 0),
      'breaks' = base::pretty(params$.ggdata[[params$.ycol]]),
      'labels' = .label_number(scale_cut = .cut_si("")) )
  }



  #________________________________________________________
  # aes() parameters
  #________________________________________________________
  set_layer(
    params = params,
    layer  = 'stack',
    'mapping|x' = params$.xcol,
    'mapping|y' = params$.ycol )



  #________________________________________________________
  # Non-aes parameters
  #________________________________________________________
  legend_title <- if (eq(params$.dots$labs.title, params$rank)) "" else params$rank
  if (!(is.null(params$colors) || isFALSE(params$colors))) {
    set_layer(params, 'labs', color = legend_title)
    set_layer(params, 'labs', fill  = legend_title)

  } else if (!(is.null(params$patterns) || isFALSE(params$patterns))) {
    set_layer(params, 'labs', fill  = legend_title)
  }
  remove("legend_title")



  #________________________________________________________
  # Create the plot and add each layer with its arguments.
  # Also attaches a human-readable version of the plot command.
  #________________________________________________________
  p <- params %>%
    plot_facets() %>%
    plot_build()

  return (p)
}
