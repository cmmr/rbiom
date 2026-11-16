#' Display taxa shared between groups as an UpSet plot.
#'
#' Draws the result of [taxa_shared()]: a bar for each combination of groups
#' showing how many taxa belong to exactly that combination, above a matrix
#' of dots that says which groups the combination is made of. Unlike a Venn
#' diagram, it stays readable with more than three groups.
#'
#' @inherit documentation_default
#' @inherit taxa_shared params details
#'
#' @family taxa_abundance
#' @family visualization
#'
#' @param stat.by   Dataset field with the groups to compare. Must be
#'        categorical, with at least two groups. Samples with `NA` in this
#'        field are left out. Required.
#'
#' @param sets   The maximum number of combinations (bars) to show. The
#'        largest are kept. Default: `15`
#'
#' @param totals   Show a bar for each group with the total number of taxa
#'        that belong to it. Default: `TRUE`
#'
#' @param ...   Parameters for the underlying ggplot2 functions, prefixed
#'        with the name of the part of the figure they apply to: `bar.`
#'        (the bars for each combination; [ggplot2::geom_col()]), `label.`
#'        (the counts above them; [ggplot2::geom_text()]), `dot.` and
#'        `empty.` (the filled and unfilled dots; [ggplot2::geom_point()]),
#'        `line.` (the lines joining dots; [ggplot2::geom_line()]), `total.`
#'        (the bars for each group; [ggplot2::geom_col()]), and
#'        `theme_rbiom.` ([theme_rbiom()], for all panels). For example,
#'        `bar.fill = "steelblue"`. An argument that no part accepts is an
#'        error.
#'
#' @return A `ggplot2` plot (a patchwork of panels). The [taxa_shared()]
#'         table is available as `$data` and `$stats`, and the code that
#'         rebuilds the figure as `$code`.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(hmp50)
#'
#'     taxa_upset(biom, stat.by = "Body Site")
#'
#'     # Core genera: present in at least half of a body site's samples.
#'     p <- taxa_upset(biom, stat.by = "Body Site", prev = 0.5, bar.fill = "steelblue")
#'     p
#'
#'     # The taxa behind the bars.
#'     p$stats
#'

taxa_upset <- function (
    biom, stat.by, rank = -1, detect = 0, prev = 0,
    unc = "singly", lineage = FALSE, sets = 15, totals = TRUE, ... ) {

  if (missing(stat.by))
    cli_abort("`stat.by` is required: the metadata field with the groups to compare.")

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('taxa_upset', params)
  schema_check('taxa_upset', params)


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('taxa_upset', params)
  p          <- get_cache_value(cache_file)

  if (is.null(p)) {

    #________________________________________________________
    # The table that the figure draws.
    #________________________________________________________
    tbl <- taxa_shared(
      biom    = biom,
      stat.by = stat.by,
      rank    = rank,
      detect  = detect,
      prev    = prev,
      unc     = unc,
      lineage = lineage )

    groups <- setdiff(names(attr(tbl, 'sets', exact = TRUE)), c('.set', '.n', '.n_groups'))

    if (length(groups) < 2)
      cli_abort(c(
        'x' = "An UpSet plot compares at least two groups, but {.field {stat.by}} has {length(groups)}: {.val {groups}}.",
        'i' = "Use {.fn taxa_shared} to list the taxa of a single group." ))

    if (nrow(tbl) == 0)
      cli_abort("No taxa belong to any group with `detect = {detect}` and `prev = {prev}`.")

    params$.fn <- 'taxa_upset()'
    p <- taxa_upset_build(tbl, groups, params)

    set_cache_value(cache_file, p)
  }

  return (plot_prov(p, cmd, biom))
}



#________________________________________________________
# The patchwork figure for a taxa_shared() table, with the code that
# rebuilds it from `data` (that table). The data-prep code below mirrors
# the computation line for line; the ggplot2 code is recorded by the P()
# wrappers as the layers are created.
#________________________________________________________
taxa_upset_build <- function (data, groups, params) {

  .set <- .group <- .n <- NULL # for CRAN check only

  n_sets     <- params$sets
  show_total <- isTRUE(params$totals)


  #________________________________________________________
  # Tables for the bars, the dot matrix, and the group totals.
  #________________________________________________________
  sets <- head(attr(data, "sets"), n_sets)
  sets[[".set"]] <- factor(as.character(sets[[".set"]]), levels = as.character(sets[[".set"]]))

  dots <- tibble::tibble(
    .set   = rep(sets[[".set"]], times = length(groups)),
    .group = factor(rep(groups, each = nrow(sets)), levels = rev(groups)),
    .in    = unlist(lapply(groups, function (g) sets[[g]])) )

  member <- dots[dots[[".in"]], ]

  totals <- tibble::tibble(
    .group = factor(groups, levels = rev(groups)),
    .n     = unname(vapply(groups, function (g) sum(data[[g]]), numeric(1))) )

  code <- c(
    "# Groups, from the first matrix row down.",
    sprintf("groups <- %s", deparse_vec(groups, indent = 4)),
    "",
    sprintf("# The %s largest sets: taxa belonging to exactly that combination of groups.", deparse_vec(n_sets)),
    sprintf("sets <- head(attr(data, \"sets\"), %s)", deparse_vec(n_sets)),
    "sets[[\".set\"]] <- factor(as.character(sets[[\".set\"]]), levels = as.character(sets[[\".set\"]]))",
    "",
    "# One dot per set and group; the filled ones are the set's groups.",
    "dots <- tibble::tibble(",
    "  .set   = rep(sets[[\".set\"]], times = length(groups)),",
    "  .group = factor(rep(groups, each = nrow(sets)), levels = rev(groups)),",
    "  .in    = unlist(lapply(groups, function (g) sets[[g]])) )",
    "",
    "member <- dots[dots[[\".in\"]], ]",
    "" )

  if (show_total)
    code <- c(
      code,
      "# Taxa belonging to each group.",
      "totals <- tibble::tibble(",
      "  .group = factor(groups, levels = rev(groups)),",
      "  .n     = unname(vapply(groups, function (g) sum(data[[g]]), numeric(1))) )",
      "" )


  #________________________________________________________
  # Provenance-tracked ggplot2 functions.
  #________________________________________________________
  .ggplot        <- P('ggplot2::ggplot')
  .geom_col      <- P('ggplot2::geom_col')
  .geom_text     <- P('ggplot2::geom_text')
  .geom_point    <- P('ggplot2::geom_point')
  .geom_line     <- P('ggplot2::geom_line')
  .labs          <- P('ggplot2::labs')
  .scale_y       <- P('ggplot2::scale_y_continuous')
  .scale_x_rev   <- P('ggplot2::scale_x_reverse')
  .scale_y_disc  <- P('ggplot2::scale_y_discrete')
  .theme         <- P('ggplot2::theme')
  .element_blank <- P('ggplot2::element_blank')
  .theme_rbiom   <- P('rbiom::theme_rbiom')


  #________________________________________________________
  # A layer's arguments: defaults, then the user's prefixed `...`.
  #________________________________________________________
  used <- character()

  layer_args <- function (prefix, ..., .formals = NULL) {

    args  <- list(...)
    regex <- paste0("^", prefix, "\\.")
    nms   <- names(params$.dots)

    for (i in intersect(nms, .formals)) {
      args[i] <- list(params$.dots[[i]])
      used   <<- c(used, i)
    }

    for (i in grep(regex, nms, value = TRUE)) {
      args[sub(regex, "", i)] <- list(params$.dots[[i]])
      used <<- c(used, i)
    }

    return (args)
  }

  # A plot and its code from a list of P() results.
  assemble <- function (var, parts) {
    p <- parts[[1]]
    for (part in parts[-1]) p <- p + part
    cmds <- vapply(parts, attr, character(1), 'display', exact = TRUE)
    list(plot = p, code = sprintf("%s <- %s", var, paste(cmds, collapse = " +\n  ")))
  }

  # Counts: no fractional axis breaks.
  whole <- aa(
    function (x) unique(floor(pretty(x))),
    display = "function (x) unique(floor(pretty(x)))" )

  theme_args <- layer_args('theme_rbiom', .formals = setdiff(formalArgs(theme_rbiom), 'frame'))


  #________________________________________________________
  # Panel 1: the number of taxa in each set.
  #________________________________________________________
  bars <- assemble("p1", list(
    .ggplot(data = structure(sets, display = "sets"), mapping = aes(x = .set, y = .n)),
    do.call(.geom_col,  layer_args('bar',   fill = "grey20", width = 0.6)),
    do.call(.geom_text, layer_args('label', mapping = aes(label = .n), vjust = -0.5, size = 3)),
    .scale_y(breaks = whole, expand = c(0, 0, 0.15, 0)),
    .labs(x = NULL, y = "Number of Taxa"),
    do.call(.theme_rbiom, theme_args),
    .theme(
      axis.text.x  = .element_blank(),
      axis.ticks.x = .element_blank() )))


  #________________________________________________________
  # Panel 2: the number of taxa in each group.
  #________________________________________________________
  if (show_total)
    side <- assemble("p2", list(
      .ggplot(data = structure(totals, display = "totals"), mapping = aes(x = .n, y = .group)),
      do.call(.geom_col, layer_args('total', fill = "grey20", width = 0.6)),
      .scale_x_rev(breaks = whole, expand = c(0.05, 0, 0, 0)),
      .scale_y_disc(position = "right"),
      .labs(x = "Taxa per Group", y = NULL),
      do.call(.theme_rbiom, theme_args),
      .theme(
        axis.text.y  = .element_blank(),
        axis.ticks.y = .element_blank() )))


  #________________________________________________________
  # Panel 3: which groups make up each set. This is the figure's
  # outer plot, so its data is the taxa_shared() table.
  #________________________________________________________
  matrix_theme <- theme_args
  matrix_theme[['frame']] <- "none"

  mtx <- assemble(if (show_total) "p3" else "p2", list(
    .ggplot(data = structure(data, display = "data")),
    do.call(.geom_point, layer_args(
      'empty',
      data    = structure(dots, display = "dots"),
      mapping = aes(x = .set, y = .group),
      color   = "grey85",
      size    = 3 )),
    do.call(.geom_line, layer_args(
      'line',
      data      = structure(member, display = "member"),
      mapping   = aes(x = .set, y = .group, group = .set),
      linewidth = 0.8 )),
    do.call(.geom_point, layer_args(
      'dot',
      data    = structure(member, display = "member"),
      mapping = aes(x = .set, y = .group),
      size    = 3 )),
    .labs(x = NULL, y = NULL),
    do.call(.theme_rbiom, matrix_theme),
    .theme(
      axis.text.x  = .element_blank(),
      axis.ticks.x = .element_blank(),
      axis.ticks.y = .element_blank() )))


  #________________________________________________________
  # Every `...` argument must have reached a layer.
  #________________________________________________________
  unused <- setdiff(names(params$.dots), used)

  if (length(unused) > 0)
    cli_abort(c(
      'x' = "{.val {unused}} {?is/are} not {?an argument/arguments} of {params$.fn} or of any of its layers.",
      'i' = "Prefix an argument with a layer name to pass it to that layer's ggplot2 function, e.g. {.code bar.fill = \"steelblue\"}. Available prefixes: {.val {paste0(c('bar', 'label', 'dot', 'empty', 'line', 'total', 'theme_rbiom'), '.')}}." ))


  #________________________________________________________
  # Bars above the matrix; group totals to its left.
  #________________________________________________________
  heights <- c(10, length(groups) + 2)

  if (show_total) {

    p <- patchwork::wrap_plots(
      A       = bars$plot,
      B       = side$plot,
      C       = mtx$plot,
      design  = "#A\nBC",
      widths  = c(1, 3),
      heights = heights )

    code <- c(
      code, bars$code, "", side$code, "", mtx$code, "",
      sprintf(
        fmt = "patchwork::wrap_plots(A = p1, B = p2, C = p3, design = \"#A\\nBC\", widths = c(1, 3), heights = %s)",
        deparse_vec(heights) ))

  } else {

    p <- patchwork::wrap_plots(
      A       = bars$plot,
      B       = mtx$plot,
      design  = "A\nB",
      heights = heights )

    code <- c(
      code, bars$code, "", mtx$code, "",
      sprintf(
        fmt = "patchwork::wrap_plots(A = p1, B = p2, design = \"A\\nB\", heights = %s)",
        deparse_vec(heights) ))
  }

  p$data <- data

  attr(p, 'stats')    <- data
  attr(p, 'code')     <- add_class(paste(code, collapse = "\n"), 'rbiom_code')
  attr(p, 'data_cmd') <- attr(data, 'cmd', exact = TRUE)

  p %<>% add_class('rbiom_plot')

  return (p)
}
