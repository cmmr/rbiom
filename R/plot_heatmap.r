#' Create a heatmap with tracks and dendrograms from any matrix.
#'
#' @inherit documentation_heatmap
#' @inherit documentation_default
#'
#' @family visualization
#'
#' @param mtx   A numeric `matrix` with named rows and columns.
#'
#' @param tracks   List of track definitions. See details below.
#'        Default: `NULL`.
#'
#' @param title   Plot title. Default: `NULL`.
#'
#' @param label_size   The font size to use for the row and column labels. You
#'        can supply a numeric vector of length two to control row label sizes
#'        and column label sizes separately, for example
#'        `c(rows = 20, cols = 8)`, or simply `c(20, 8)`.
#'        Default: `NULL`, which uses the theme's tick label size, shrunk to
#'        `pmax(5, 100 / dim(mtx))` when there are too many rows or columns
#'        for it.
#'
#' @param ...   Additional arguments to pass on to [ggplot2::theme()].
#'        Arguments prefixed with `labs.` go to [ggplot2::labs()] instead,
#'        for example `labs.subtitle = "Plot subtitle"`, and `base_size`,
#'        `base_family`, `base_line_size`, and `base_rect_size` go to
#'        [theme_rbiom()].
#'
#'
#' @section Track Definitions:
#'
#' One or more colored tracks can be placed on the left and/or top of the
#' heatmap grid to visualize associated metadata values.
#'
#' \preformatted{## Categorical ----------------------------
#' cat_vals <- sample(c("Male", "Female"), 10, replace = TRUE)
#' tracks   <- list('Sex' = cat_vals)
#' tracks   <- list('Sex' = list(values = cat_vals, colors = "bright"))
#' tracks   <- list('Sex' = list(
#'   values = cat_vals,
#'   colors = c('Male' = "blue", 'Female' = "red")) )
#'
#' ## Numeric --------------------------------
#' num_vals <- sample(25:40, 10, replace = TRUE)
#' tracks   <- list('Age' = num_vals)
#' tracks   <- list('Age' = list(values = num_vals, colors = "greens"))
#' tracks   <- list('Age' = list(values = num_vals, range = c(0,50)))
#' tracks   <- list('Age' = list(
#'   label  = "Age (Years)",
#'   values = num_vals,
#'   colors = c("azure", "darkblue", "darkorchid") ))
#'
#' ## Multiple Tracks ------------------------
#' tracks <- list('Sex' = cat_vals, 'Age' = num_vals)
#' tracks <- list(
#'   list(label = "Sex", values = cat_vals, colors = "bright"),
#'   list(label = "Age", values = num_vals, colors = "greens") )
#'
#'
#' mtx           <- matrix(sample(1:50), ncol = 10)
#' dimnames(mtx) <- list(letters[1:5], LETTERS[1:10])
#' plot_heatmap(mtx = mtx, tracks = tracks)
#' }
#'
#' The following entries in the track definitions are understood:
#'
#' \describe{
#'   \item{\code{values} - }{ The metadata values. When unnamed, order must match matrix. }
#'   \item{\code{range} - }{ The c(min,max) to use for scale values. }
#'   \item{\code{label} - }{ Label for this track. Defaults to the name of this list element. }
#'   \item{\code{side} - }{ Options are \code{"top"} (default) or \code{"left"}. }
#'   \item{\code{colors} - }{ A pre-defined palette name or custom set of colors to map to. }
#'   \item{\code{na.color} - }{ The color to use for \code{NA} values. }
#'   \item{\code{bins} - }{ Bin a gradient into this many bins/steps. }
#'   \item{\code{guide} - }{ A list of arguments for guide_colorbar() or guide_legend(). }
#' }
#'
#' The default categorical colors are colorblind-friendly for up to 12
#' groups: Okabe and Ito's palette (`"okabe"`) for up to 8, and Paul Tol's
#' (`"tol"`) for 9-12.
#'
#' The available categorical palette names are: `"okabe"`, `"carto"`, `"r4"`,
#' `"polychrome"`, `"tol"`, `"bright"`, `"light"`, `"muted"`, `"vibrant"`,
#' `"tableau"`, `"classic"`, `"alphabet"`, `"tableau20"`, `"kelly"`, and
#' `"fishy"`.
#'
#' The available numeric palette names are: `"reds"`, `"oranges"`,
#' `"greens"`, `"purples"`, `"grays"`, `"acton"`, `"bamako"`, `"batlow"`,
#' `"bilbao"`, `"buda"`, `"davos"`, `"devon"`, `"hawaii"`, `"imola"`,
#' `"lajolla"`, `"lapaz"`, `"nuuk"`, `"oslo"`, `"tokyo"`, `"turku"`, `"bam"`,
#' `"berlin"`, `"broc"`, `"cork"`, `"lisbon"`, `"roma"`, `"tofino"`,
#' `"vanimo"`, and `"vik"`.
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     set.seed(123)
#'     mtx <- matrix(runif(5*8), nrow = 5, dimnames = list(LETTERS[1:5], letters[1:8]))
#'
#'     plot_heatmap(mtx)
#'     plot_heatmap(mtx, grid="oranges")
#'     plot_heatmap(mtx, grid=list(colors = "oranges", label = "Litres"))
#'
#'     # Tracks go on top (one value per column) unless side = "left"
#'     # (one value per row).
#'     tracks <- list(
#'       'Number' = sample(1:ncol(mtx)),
#'       'Person' = list(
#'         values = factor(sample(c("Alice", "Bob"), ncol(mtx), TRUE)),
#'         colors = c('Alice' = "purple", 'Bob' = "darkcyan") ),
#'       'State' = list(
#'         side   = "left",
#'         values = sample(c("TX", "OR", "WA"), nrow(mtx), TRUE),
#'         colors = "bright" )
#'     )
#'
#'     plot_heatmap(mtx, tracks=tracks)
#'
plot_heatmap <- function (
    mtx, grid = list(label = "Grid Value", colors = "imola"), tracks = NULL,
    label = TRUE, label_size = NULL, rescale = "none",
    trees = TRUE, clust = "complete", dist = "euclidean",
    asp = 1, tree_height = 10, track_height = 10,
    legend = "right", title = NULL, xlab.angle = "auto", ... ) {

  params <- eval_envir(environment(), ...)

  cmd <- local({
    p <- as.list(params)
    p$mtx <- aa(NA, display = "mtx")
    if (inherits(p$dist,  'dist'))   p$dist  <- aa(NA, display = "dist")
    if (inherits(p$clust, 'hclust')) p$clust <- aa(NA, display = "clust")
    current_cmd('plot_heatmap', p)
  })


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- local({
    key <- as.list(params, all.names = TRUE)
    key$mtx <- strip_prov(mtx)
    get_cache_file('plot_heatmap', key)
  })

  p <- get_cache_value(cache_file)

  if (is.null(p)) {
    params$.fn <- 'plot_heatmap()'
    p <- set_cache_value(cache_file, plot_heatmap_build(params))
  }

  attr(p, 'cmd')     <- cmd
  attr(p, 'history') <- attr(mtx, 'history', exact = TRUE)

  return (p)
}




#________________________________________________________
# plot_heatmap() without the cache / provenance wrapper.
#
# Inputs read from `params`: the plot_heatmap() arguments, plus
# optionally
#   .mtx_code   - code lines that create `mtx` (default: the matrix's own
#                 `cmd` attribute, else a literal matrix).
#   .dist_code  - code for the `dist` object when one was given.
#   .data_var / .plot_var - see plot_build().
#________________________________________________________
plot_heatmap_build <- function (params) {

  stopifnot(is_bare_environment(params))

  dv <- params$.data_var %||% "data"
  pv <- params$.plot_var

  mtx          <- params$mtx
  grid         <- params$grid
  tracks       <- params$tracks
  label        <- params$label
  label_size   <- params$label_size
  rescale      <- params$rescale
  trees        <- params$trees
  clust        <- params$clust
  dist         <- params$dist
  asp          <- params$asp
  tree_height  <- params$tree_height
  track_height <- params$track_height
  legend       <- params$legend
  title        <- params$title
  xlab.angle   <- params$xlab.angle
  dots         <- params$.dots


  #________________________________________________________
  # Sanity Checks.
  #________________________________________________________
  if (!is.matrix(mtx) || !is.numeric(mtx))
    cli_abort("`mtx` must be a numeric matrix, not {.type {mtx}}.")
  if (is.null(rownames(mtx)) || is.null(colnames(mtx)))
    cli_abort("`mtx` must have row and column names.")
  if (nrow(mtx) < 1 || ncol(mtx) < 1)
    cli_abort("`mtx` must have at least one row and one column.")

  stopifnot(
    (is_logical(label) && length(label) <= 2) ||
    is_string(label, c("rows", "cols", "both", "bottom", "right", "none")) )
  stopifnot(
    (is_logical(trees) && length(trees) <= 2) ||
    is_string(trees, c("rows", "cols", "both", "left", "top", "none")) )

  if (!any(is_na(clust), is_false(clust), inherits(clust, 'hclust')))
    validate_clust(max = 2, na_ok = TRUE)

  if (!inherits(dist, 'dist'))
    validate_dist(max = 2)

  stopifnot(is_null(label_size)   || (is.numeric(label_size)   && length(label_size)   <= 2) )
  stopifnot(is_null(tree_height)  || (is.numeric(tree_height)  && length(tree_height)  <= 2))
  stopifnot(is_null(track_height) || (is.numeric(track_height) && length(track_height) <= 2))
  stopifnot(is.numeric(asp) && length(asp) == 1 && !is_na(asp))
  stopifnot(is_string(rescale, c("none", "rows", "cols")))
  stopifnot(is_string(legend, c("right", "bottom")))
  stopifnot(is_string(as.character(xlab.angle), c("auto", "0", "30", "90")))

  if (!(is.null(title) || (is_scalar_character(title) && !is.na(title))))
    cli_abort("`title` must be NULL or a character string, not {.type {title}}: {title}.")


  #________________________________________________________
  # `...` arguments: labs.* go to labs(), base_* to theme_rbiom(),
  # and the rest to theme().
  #________________________________________________________
  labs_args  <- list()
  rbiom_args <- list(frame = "none")
  theme_args <- list()

  for (i in names(dots)) {
    if (startsWith(i, "labs.")) { labs_args[[sub("^labs\\.", "", i)]] <- dots[[i]] }
    else if (startsWith(i, "base_")) { rbiom_args[[i]] <- dots[[i]] }
    else                        { theme_args[[i]] <- dots[[i]] }
  }

  if (length(bad <- setdiff(names(rbiom_args), formalArgs(theme_rbiom))) > 0)
    cli_abort("{.val {bad}} {?is/are} not {?an argument/arguments} of {.fn theme_rbiom}.")

  if (length(bad <- setdiff(names(theme_args), formalArgs(ggplot2::theme))) > 0)
    cli_abort(c(
      'x' = "{.val {bad}} {?is/are} not {?an argument/arguments} of {.fn ggplot2::theme}.",
      'i' = "Prefix plot labels with `labs.`, e.g. `labs.subtitle = \"...\"`." ))

  if (!is.null(title)) labs_args[['title']] <- title


  #________________________________________________________
  # Provenance-tracked ggplot2 functions.
  #________________________________________________________
  .ggplot             <- P('ggplot2::ggplot')
  .coord_fixed        <- P('ggplot2::coord_fixed')
  .labs               <- P('ggplot2::labs')
  .geom_raster        <- P('ggplot2::geom_raster')
  .scale_x_discrete   <- P('ggplot2::scale_x_discrete')
  .scale_y_discrete   <- P('ggplot2::scale_y_discrete')
  .annotate           <- P('ggplot2::annotate')
  .geom_segment       <- P('ggplot2::geom_segment')
  .geom_tile          <- P('ggplot2::geom_tile')
  .element_text       <- P('ggplot2::element_text')
  .theme              <- P('ggplot2::theme')
  .theme_rbiom        <- P('rbiom::theme_rbiom')
  .new_scale_fill     <- P('ggnewscale::new_scale_fill')
  .element_rect       <- P('ggplot2::element_rect')
  .margin             <- P('ggplot2::margin')
  .scale_x_continuous <- P('ggplot2::scale_x_continuous')
  .scale_y_continuous <- P('ggplot2::scale_y_continuous')
  .sec_axis           <- P('ggplot2::sec_axis')


  #________________________________________________________
  # Setting legend = "bottom" tiggers a bunch of presets
  #________________________________________________________
  if (eq(legend, "bottom")) {
    theme_args[['legend.box']]        %<>% if.null("horizontal")
    theme_args[['legend.position']]   %<>% if.null("bottom")
    theme_args[['legend.background']] %<>% if.null(.element_rect(color = "black"))
    theme_args[['legend.margin']]     %<>% if.null(.margin(6,12,6,12,"points"))
    theme_args[['legend.box.margin']] %<>% if.null(.margin(6,6,6,6,"points"))
  }


  #________________________________________________________
  # Allow different specs for rows and cols.
  # For example, dist = c("euclidean", "canberra")
  #________________________________________________________

  # for CRAN check only
  label_rows <- label_cols <- NULL
  trees_rows <- trees_cols <- NULL
  clust_rows <- clust_cols <- NULL
  dist_rows  <- dist_cols  <- NULL

  for (i in c("label", "trees", "clust", "dist")) {
    val <- get(i)

    # Allow specifying label = "Rows", trees = "both", etc
    if (i %in% c("label", "trees") && is.character(val) && length(val) == 1) {
      val <- tolower(val)
      val <- if (val == "both")                                  { TRUE
      } else if (val %in% c("rows", "left", "right"))            { c(TRUE, FALSE)
      } else if (val %in% c("cols", "columns", "top", "bottom")) { c(FALSE, TRUE)
      } else                                                     { FALSE }
    }

    if (length(val) == 2 && !inherits(val, c('dist', 'hclust'))) {

      if (!is_null(names(val)))
        val <- val[rev(order(names(val)))]

      assign(sprintf("%s_rows", i), val[[1]])
      assign(sprintf("%s_cols", i), val[[2]])
    } else {
      assign(sprintf("%s_rows", i), val)
      assign(sprintf("%s_cols", i), val)
    }
    remove("val")
  }


  #________________________________________________________
  # The code that creates `mtx`.
  #________________________________________________________
  code <- params$.mtx_code %||% local({

    mtx_cmd <- attr(mtx, 'cmd', exact = TRUE)

    if (!is.null(mtx_cmd))
      return (sprintf("mtx <- %s", mtx_cmd))

    c("# The matrix as given to plot_heatmap().",
      sprintf("mtx <- matrix(%s, nrow = %d, dimnames = list(%s, %s))",
        deparse_vec(as.vector(mtx), indent = 2), nrow(mtx), deparse_vec(rownames(mtx), indent = 2), deparse_vec(colnames(mtx), indent = 2) ))
  })

  mtx <- matrix(as.vector(mtx), nrow = nrow(mtx), dimnames = dimnames(mtx))


  #________________________________________________________
  # Combine heights + aspect ratio to get top/left heights
  #________________________________________________________
  stopifnot(length(tree_height) > 0 && length(track_height) > 0)

  asp <- asp * ncol(mtx) / nrow(mtx)

  tree_height_top   <- tree_height[[1]]  / 100 * nrow(mtx)
  track_height_top  <- track_height[[1]] / 100 * nrow(mtx)
  tree_height_left  <- if (length(tree_height)  == 1) tree_height_top  * asp else tree_height[[2]]  / 100 * ncol(mtx)
  track_height_left <- if (length(track_height) == 1) track_height_top * asp else track_height[[2]] / 100 * ncol(mtx)
  remove("tree_height", "track_height")


  #________________________________________________________
  # Heuristics for setting sizes and heights of elements
  #________________________________________________________
  # NA = the theme's size for tick labels (9.35 pt).
  if (is_null(label_size))         label_size <- ifelse(100 / dim(mtx) >= 9, NA, pmax(5, 100 / dim(mtx)))
  if (!is_null(names(label_size))) label_size <- label_size[rev(order(names(label_size)))]
  label_size_x <- signif(tail(label_size, 1), digits = 3)
  label_size_y <- signif(head(label_size, 1), digits = 3)

  # "auto" angle for the column labels, by how much room they need.
  if (eq(as.character(xlab.angle), "auto")) {
    width      <- sum(nchar(colnames(mtx))) + 2 * ncol(mtx)
    xlab.angle <- if (width <= 50) "0" else if (width <= 100) "30" else "90"
    remove("width")
  }

  if (is.null(tracks) || length(tracks) == 0) {
    n_left <- n_top <- 0
  } else {
    n_left <- sum(sapply(tracks, function (x) is.list(x) && eq(x[['side']], 'left')))
    n_top  <- length(tracks) - n_left
  }


  #________________________________________________________
  # Rescale matrix's row or column values
  #________________________________________________________
  if (eq(rescale, "rows")) {
    mtx  <- t(apply(mtx, 1L, scales::rescale))
    code <- c(code, "mtx <- t(apply(mtx, 1L, scales::rescale))")

  } else if (eq(rescale, "cols")) {
    mtx  <- apply(mtx, 2L, scales::rescale)
    code <- c(code, "mtx <- apply(mtx, 2L, scales::rescale)")
  }


  #________________________________________________________
  # Cluster matrix's rows/cols and generate data.frames for dendrograms
  #________________________________________________________

  # The code for a dist / hclust object given by the user.
  obj_code <- function (x, what) {
    if (!is.null(cc <- attr(x, 'cmd', exact = TRUE)))        return (list(expr = cc,   pre = NULL))
    if (!is.null(cc <- params[[paste0('.', what, '_code')]])) return (list(expr = cc,   pre = NULL))
    list(expr = what, pre = sprintf("# %s <- %s   # not reproducible from code", what, obj_repr(x)))
  }

  cluster <- function (dm, hc, margin) {

    # margin 1 = rows, 2 = cols
    mtx_expr <- if (margin == 1) "mtx" else "t(mtx)"
    var      <- if (margin == 1) "hc_rows" else "hc_cols"
    lines    <- NULL

    if (!inherits(hc, 'hclust')) {

      if (!inherits(dm, 'dist')) {
        dm_expr <- sprintf("stats::dist(%s, method = %s)", mtx_expr, double_quote(dm))
        dm      <- stats::dist(if (margin == 1) mtx else t(mtx), method = dm)
      } else {
        oc      <- obj_code(dm, 'dist')
        dm_expr <- oc$expr
        lines   <- c(lines, oc$pre)
      }

      lines <- c(lines, sprintf("%s <- stats::hclust(%s, method = %s)", var, dm_expr, double_quote(hc)))
      hc    <- stats::hclust(dm, method = hc)

    } else {
      oc    <- obj_code(hc, 'clust')
      lines <- c(lines, oc$pre, sprintf("%s <- %s", var, oc$expr))
    }

    list(hc = hc, lines = lines, var = var)
  }


  if (clust_wanted(clust_cols) && ncol(mtx) > 1) {

    cl   <- cluster(dist_cols, clust_cols, margin = 2)
    hc   <- cl$hc
    mtx  <- mtx[, hc[['order']], drop = FALSE]
    code <- c(code, cl$lines, "mtx <- mtx[, hc_cols$order, drop = FALSE]")

    if (isTRUE(trees_cols)) {

      bounds <- 0.5 + nrow(mtx) + c(0, tree_height_top)
      if (isTRUE(n_top > 0))
        bounds <- bounds + (nrow(mtx) / 50) + track_height_top * n_top

      trees_cols <- dendro(hc, bounds, side = "top")
    }
    remove("cl", "hc")

  } else {
    trees_cols <- FALSE
  }

  if (clust_wanted(clust_rows) && nrow(mtx) > 1) {

    cl   <- cluster(dist_rows, clust_rows, margin = 1)
    hc   <- cl$hc
    mtx  <- mtx[hc[['order']], , drop = FALSE]
    code <- c(code, cl$lines, "mtx <- mtx[hc_rows$order, , drop = FALSE]")

    if (isTRUE(trees_rows)) {

      bounds <- 0.5 - c(0, tree_height_left)
      if (isTRUE(n_left > 0))
        bounds <- bounds - (ncol(mtx) / 50) - track_height_left * n_left

      trees_rows <- dendro(hc, bounds, side = "left")
    }
    remove("cl", "hc")

  } else {
    trees_rows <- FALSE
  }

  all_trees <- rbind(
    if (!isFALSE(trees_cols)) trees_cols else NULL,
    if (!isFALSE(trees_rows)) trees_rows else NULL )

  remove("trees_cols", "trees_rows")
  remove("clust", "clust_rows", "clust_cols")
  remove("dist", "dist_rows", "dist_cols")



  #________________________________________________________
  # Convert all track definitions to long form.
  #________________________________________________________
  top_tracks  <- list()
  left_tracks <- list()
  for (i in seq_along(tracks)) {

    track <- tracks[[i]]
    if (!is_list(track)) track <- list(values = track)


    # Defaults and sanity checks.
    #________________________________________________________
    track[['label']] %<>% if.null(names(tracks)[[i]])
    if (!is_scalar_character(track[['label']])) track[['label']] <- paste("Track", i)
    if (!is_string(track[['side']], c("left", "top"))) track[['side']] <- "top"

    if (is_null(track[['values']])) cli_abort("`values` needed for track {.val {track[['label']]}}.")


    # Ensure that values is the same length as rows or cols.
    #________________________________________________________
    track[['values']] <- local({

      values     <- track[['values']]
      side_len   <- if (track[['side']] == "top") ncol(mtx)     else nrow(mtx)
      side_names <- if (track[['side']] == "top") colnames(mtx) else rownames(mtx)

      if (is_null(names(values))) {
        if (!isTRUE(length(values) == side_len))
          cli_abort("Length mismatch for track {.val {track[['label']]}}: expected {side_len} values, got {length(values)}.")

      } else {
        if (length(missing <- setdiff(side_names, names(values))) > 0)
          cli_abort("Track {.val {track[['label']]}} is missing values for {length(missing)} row/col name{?s}: {.val {head(missing)}}.")
        values <- as.vector(values[side_names])
      }

      if (!is.numeric(values))
        values %<>% as.factor()

      return (values)
    })


    # Define the colors for this track.
    #________________________________________________________
    track[['colors']] <- local({

      colors <- track[['colors']]

      if (is.numeric(track[['values']])) {

        colors %<>% if.null(c('reds', 'oranges', 'greens', 'purples')[(i %% 4) + 1])
        if (is_palette(colors)) colors <- get_palette(colors)

      } else {

        keys <- levels(track[['values']])

        if (is_null(names(colors))) {
          colors <- get_n_colors(length(keys), colors)

        } else {
          if (length(missing <- setdiff(keys, names(colors))) > 0)
            cli_abort("Track {.val {track[['label']]}} is missing colors for {length(missing)} value{?s}: {.val {head(missing)}}.")
          colors <- as.vector(colors[keys])
        }
      }

      return (colors)
    })


    # Add to either top or left collection.
    #________________________________________________________
    if (eq(track[['side']], "top")) { top_tracks  %<>% c(list(track))
    } else                          { left_tracks %<>% c(list(track)) }

  }



  #________________________________________________________
  # Generate data.frames for annotation tracks
  #________________________________________________________
  if (isTRUE(n_top > 0))
    top_tracks %<>% tracks_df(
      bounds = nrow(mtx) + 0.5 + (nrow(mtx) / 50) + c(0, track_height_top * n_top),
      side   = "top" )

  if (isTRUE(n_left > 0))
    left_tracks %<>% tracks_df(
      bounds = 0.5 - (ncol(mtx) / 50) - c(0, track_height_left * n_left),
      side   = "left" )

  all_tracks <- c(top_tracks, left_tracks)

  remove("tree_height_left", "tree_height_top", "track_height_left", "track_height_top")


  #________________________________________________________
  # The data for ggplot(): the matrix in long form.
  #________________________________________________________
  x <- y <- fill <- xend <- yend <- NULL # for CRAN check only

  df <- tibble(
    x    = factor(rep(colnames(mtx), each = nrow(mtx)), levels = colnames(mtx)),
    y    = factor(rep(rownames(mtx),        ncol(mtx)), levels = rownames(mtx)),
    fill = as.vector(mtx) )

  code <- c(
    code,
    "",
    sprintf("%s <- tibble::tibble(", dv),
    "  x    = factor(rep(colnames(mtx), each = nrow(mtx)), levels = colnames(mtx)),",
    "  y    = factor(rep(rownames(mtx),        ncol(mtx)), levels = rownames(mtx)),",
    "  fill = as.vector(mtx) )" )


  if (!is_null(all_trees)) {
    all_trees <- as_tibble(all_trees)
    code <- c(code, "", derived_literal_code("trees", all_trees, comment = "# Dendrogram segments (from hc_rows / hc_cols) in the grid's coordinates."))
  }

  for (i in seq_along(all_tracks)) {

    track <- all_tracks[[i]]

    values_code <- track[['values_code']] %||% deparse_vec(track[['values']], indent = 2)

    code <- c(code, "", sprintf(
      fmt = "%s <- tibble::tibble(fill = %s, %s = seq_along(fill))",
      track[['id']], values_code, if (eq(track[['side']], "top")) "x" else "y" ))
  }



  #________________________________________________________
  # Assemble the layers for the ggplot
  #________________________________________________________
  gglayers <- list()
  gglayers %<>% ggpush(.ggplot(structure(df, display = dv)))
  gglayers %<>% ggpush(.coord_fixed(ratio = asp))

  if (length(labs_args) > 0)
    gglayers %<>% ggpush(do.call(.labs, c(labs_args, .indent = 4)))


  #________________________________________________________
  # In order to label any tracks, switch to continuous
  #________________________________________________________
  if (length(all_tracks) == 0) {

    gglayers %<>% ggpush(.geom_raster(aes(x = x, y = y, fill = fill)))

    if (isFALSE(label_cols))
      gglayers %<>% ggpush(.scale_x_discrete(breaks = NULL))

    if (isTRUE(label_rows)) {
      gglayers %<>% ggpush(.scale_y_discrete(position = 'right'))
    } else {
      gglayers %<>% ggpush(.scale_y_discrete(breaks = NULL))
    }

  } else {

    gglayers %<>% ggpush(.geom_raster(aes(x = as.numeric(x), y = as.numeric(y), fill = fill)))


    #________________________________________________________
    # Column labels along the bottom; left tracks are labeled
    # along the top with a secondary axis.
    #________________________________________________________
    x_args <- list(expand = c(0.01, 0), .indent = 4)

    if (isTRUE(label_cols)) {
      x_args[['breaks']] <- structure(seq_len(ncol(mtx)), display = "seq_len(ncol(mtx))")
      x_args[['labels']] <- structure(colnames(mtx),      display = "colnames(mtx)")
    } else {
      x_args['breaks'] <- list(NULL)
    }

    if (isTRUE(n_left > 0))
      x_args[['sec.axis']] <- .sec_axis(
        transform = ~ .,
        breaks = sapply(left_tracks, function (x) x[['label_at']]),
        labels = sapply(left_tracks, function (x) x[['label']]) )


    #________________________________________________________
    # Row labels along the right; top tracks are labeled
    # along the left with a secondary axis.
    #________________________________________________________
    y_args <- list(expand = c(0.01, 0), position = "right", .indent = 4)

    if (isTRUE(label_rows)) {
      y_args[['breaks']] <- structure(seq_len(nrow(mtx)), display = "seq_len(nrow(mtx))")
      y_args[['labels']] <- structure(rownames(mtx),      display = "rownames(mtx)")
    } else {
      y_args['breaks'] <- list(NULL)
    }

    if (isTRUE(n_top > 0))
      y_args[['sec.axis']] <- .sec_axis(
        transform = ~ .,
        breaks = sapply(top_tracks, function (x) x[['label_at']]),
        labels = sapply(top_tracks, function (x) x[['label']]) )


    gglayers %<>% ggpush(do.call(.scale_x_continuous, x_args))
    gglayers %<>% ggpush(do.call(.scale_y_continuous, y_args))

    remove("x_args", "y_args")
  }



  #________________________________________________________
  # Outline around the main grid
  #________________________________________________________
  gglayers %<>% ggpush(.annotate(
    geom      = "rect",
    xmin      = 0.5,
    xmax      = ncol(mtx) + 0.5,
    ymin      = 0.5,
    ymax      = nrow(mtx) + 0.5,
    linewidth = 0.2,
    color     = "black",
    fill      = NA ))



  #________________________________________________________
  # Apply color gradient to the main grid
  #________________________________________________________
  if (is_null(grid))             grid <- list()
  if (is_scalar_character(grid)) grid <- list(colors = grid)
  if (!is_list(grid))            grid <- as.list(grid)

  grid[['label']]  %<>% if.null("Grid Value")
  grid[['colors']] %<>% if.null("imola")
  if (is_palette(grid[['colors']]))
    grid[['colors']] <- get_palette(grid[['colors']])

  layer <- plot_heatmap_layer(grid, length(all_tracks) + 1)

  layer[['args']][['.indent']] <- 4
  gglayers %<>% ggpush(do.call(layer[['fun']], layer[['args']]))

  remove("layer")



  if (!is_null(all_trees)) {

    gglayers %<>% ggpush(.geom_segment(
      data    = structure(all_trees, display = "trees"),
      mapping = aes(x = x, y = y, xend = xend, yend = yend),
      .indent = 4 ))
  }


  for (i in seq_along(all_tracks)) {

    track <- all_tracks[[i]]
    layer <- plot_heatmap_layer(track, i)

    gglayers %<>% ggpush(.new_scale_fill())
    gglayers %<>% ggpush(.geom_tile(data = structure(track[['data']], display = track[['id']]), mapping = track[['mapping']], .indent = 4))
    gglayers %<>% ggpush(do.call(layer[['fun']], c(layer[['args']], .indent = 4)))
    gglayers %<>% ggpush(do.call(.annotate, c(list(geom = "rect", color = "black", fill = NA, linewidth = 0.2), track[['outline']])))
  }


  #________________________________________________________
  # Theming. Track labels use the theme's tick label size; grid
  # labels use `label_size` (NA for the theme's size).
  #________________________________________________________
  args <- theme_args

  x_text <- function (size = NA) {
    el <- switch(
      EXPR = paste0("x", tolower(xlab.angle)),
      x0   = list(),
      x90  = list(angle = -90, vjust = 0.3, hjust = 0),
             list(angle = -30, vjust = 1.0, hjust = 0) )
    if (!is.na(size)) el[['size']] <- size
    if (length(el) == 0) NULL else do.call(.element_text, el)
  }

  if (!hasName(args, 'axis.text.y.right') && isTRUE(label_rows) && !is.na(label_size_y))
    args[['axis.text.y.right']] <- .element_text(size = label_size_y)

  if (!hasName(args, 'axis.text.x.bottom') && isTRUE(label_cols))
    args[['axis.text.x.bottom']] <- x_text(label_size_x)

  if (!hasName(args, 'axis.text.x.top') && isTRUE(n_left > 0))
    args[['axis.text.x.top']] <- x_text()


  gglayers %<>% ggpush(do.call(.theme_rbiom, rbiom_args))
  if (length(args) > 0)
    gglayers %<>% ggpush(do.call(.theme, c(args, .indent = 4)))

  remove("args")


  #________________________________________________________
  # Build the plot and its code.
  #________________________________________________________
  p <- ggbuild(gglayers)
  p$data <- df

  plot_code <- attr(p, 'code', exact = TRUE)
  if (!is.null(pv)) plot_code <- sprintf("%s <- %s", pv, plot_code)

  params$.plot_attrs %<>% if.null(list())

  return (plot_finalize(p, params, c(code, "", plot_code), df))
}




#________________________________________________________
# Is clustering requested? (`FALSE` and `NA` disable it; a method name
# or an hclust object enable it.)
#________________________________________________________
clust_wanted <- function (x) {
  if (inherits(x, 'hclust')) return (TRUE)
  !(isFALSE(x) || (is.atomic(x) && length(x) == 1 && is.na(x)))
}



#________________________________________________________
# Convert palette args to scale_fill_* function
#________________________________________________________
plot_heatmap_layer <- function (args, guide_order = 0) {

  colors  <- args[['colors']]
  guide   <- as.list(args[['guide']])
  binned  <- !is_null(args[['bins']])
  numeric <- is_null(args[['values']]) || is.numeric(args[['values']])

  .theme          <- P('ggplot2::theme')
  .guide_colorbar <- P('ggplot2::guide_colorbar')
  .guide_legend   <- P('ggplot2::guide_legend')


  # fn_args gets passed on to scale_*, e.g. scale_fill_gradient
  fn_args <- list()
  if (!is_null(args[['label']]))    fn_args[['name']]     <- args[['label']]
  if (!is_null(args[['range']]))    fn_args[['limits']]   <- args[['range']]
  if (!is_null(args[['bins']]))     fn_args[['n.breaks']] <- args[['bins']] + 1
  if (!is_null(args[['na.color']])) fn_args[['na.value']] <- args[['na.color']]

  fn <- NULL

  if (numeric) {

    # Create color gradient from >= 2 color codes
    if (length(colors) == 2) {
      fn <- if (binned) "scale_fill_steps" else "scale_fill_gradient"
      fn_args[['low']]  <- colors[[1]]
      fn_args[['high']] <- colors[[2]]

    } else {
      fn <- if (binned) "scale_fill_stepsn" else "scale_fill_gradientn"
      fn_args[['colours']]  <- colors
    }

    guide[['.indent']]   %<>% if.null(6)
    guide[['order']]     %<>% if.null(guide_order)
    guide[['direction']] %<>% if.null('horizontal')
    guide[['theme']]     %<>% if.null(.theme(legend.title.position = "top"))
    fn_args[['guide']] <- do.call(.guide_colorbar, guide)

  } else {

    # Create a discrete color scale
    fn <- "scale_fill_manual"
    fn_args[['values']] <- colors

    guide[['order']] %<>% if.null(guide_order)
    guide[['ncol']]  %<>% if.null(2)
    fn_args[['guide']] <- do.call(.guide_legend, guide)
  }

  fun <- P(paste0('ggplot2::', fn))


  return (list(fn = fn, fun = fun, args = fn_args))
}



#________________________________________________________
# Convert user's `tracks` into a named list of lists, with the sample
# metadata values and the code that fetches them.
#________________________________________________________
biom_tracks <- function (params, side = NULL) {

  stopifnot(is_bare_environment(params))

  tracks <- params$tracks
  df     <- params$biom$metadata

  if (is.null(tracks) || is.character(tracks))
    tracks <- as.list(tracks)

  if (!is.list(tracks))
    cli_abort("Invalid `tracks` specification.")

  keys <- names(tracks)
  for (i in seq_along(tracks)) {

    key <- keys[[i]] %||% ""
    val <- tracks[[i]]

    if (is.character(val) && !is.na(val)) {
      if (key == "" && length(val) == 1) {
        key <- val
        val <- list()
      } else {
        val <- list(colors = val)
      }
    }

    if (!is.list(val))
      cli_abort("Invalid `tracks` specification.")


    if (is.null(names(val$colors))) {

      validate_df_field('key', evar = 'tracks', drop_na = FALSE)

    } else {

      # allow abbreviations in color mapping, e.g. sex = c(m = "blue", f = "red")
      validate_df_field('key', evar = 'tracks', col_type = "cat", drop_na = FALSE)
      nms <- names(val$colors)
      validate_var_choices('nms', choices = levels(df[[key]]), evar = key, max = Inf)
      names(val$colors) <- nms
    }

    if (!is.null(side)) val$side <- side

    values <- setNames(df[[key]], df$.sample)

    names(tracks)[[i]]       <- key
    tracks[[i]]              <- val
    tracks[[i]]$values       <- values
    tracks[[i]]$values_code  <- sprintf(
      fmt = "%s(unname(pull(biom, %s)[%s]))",
      if (is.numeric(values)) "as.numeric" else "as.factor",
      double_quote(key),
      if (eq(val$side, 'left')) "rownames(mtx)" else "colnames(mtx)" )
  }

  params$tracks <- tracks

  return (invisible(params))
}
