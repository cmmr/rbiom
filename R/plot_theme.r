#' Publication-ready ggplot2 theme.
#'
#' The theme that every rbiom plot uses. It follows common journal figure
#' guidelines: a sans-serif font, black text, open x/y axis lines with
#' outward ticks, no gridlines or background shading, and facet labels
#' without boxes. Titles and captions align with the left edge of the figure,
#' and secondary text (subtitles and captions) is gray.
#'
#' rbiom's plot functions accept these arguments in their `...`, for example
#' `adiv_boxplot(biom, base_size = 8, grid = "y")` (heatmaps accept the
#' `base_*` ones). Given there, `base_size` also scales the plot's
#' annotations (p-values, taxa labels) to match.
#'
#' The sizes of all text elements are relative to `base_size`, so the text
#' of an existing plot can also be resized in one step with
#' `p + ggplot2::theme(text = ggplot2::element_text(size = 8))`. The same
#' works for `family`.
#'
#' @param base_size   Font size in points for axis titles, legend titles, and
#'        facet labels. Tick labels, legend labels, and captions are 85% or
#'        80% of this size; titles are 115%. Most journals want 6-8 pt text at
#'        the final printed size: for a figure saved at its final width
#'        (e.g. 7 inches for two columns), set `base_size = 8`.
#'        Default: `11`
#'
#' @param base_family   Font family. The default, `""`, is the graphics
#'        device's sans-serif font: Helvetica for [grDevices::pdf()], Arial
#'        on Windows. Default: `""`
#'
#' @param frame   How to outline the plotting area. Options are:
#'        \describe{
#'            \item{`"axes"` - }{ x and y axis lines with outward ticks. }
#'            \item{`"box"` - }{ A box around each panel with outward ticks. }
#'            \item{`"none"` - }{ No axis lines, ticks, or axis titles, as in
#'                                [ggplot2::theme_void()] (for heatmaps). }
#'        }
#'        Default: `"axes"`
#'
#' @param grid   Light gray major gridlines: `FALSE` (none), `TRUE` (both
#'        axes), `"x"` (vertical lines), or `"y"` (horizontal lines).
#'        Default: `FALSE`
#'
#' @param base_line_size   Width of axis lines and ticks, in mm (as for
#'        ggplot2's `linewidth`). Default: `base_size / 22`
#'
#' @param base_rect_size   Width of the panel box, in mm. Default:
#'        `base_size / 22`
#'
#' @return A complete [ggplot2::theme()] object.
#'
#' @family visualization
#' @export
#' @examples
#'     library(rbiom)
#'     library(ggplot2)
#'
#'     # The theme of every rbiom plot.
#'     biom <- rarefy(hmp50)
#'     adiv_boxplot(biom, x = "Body Site", stat.by = "Body Site")
#'
#'     # Make a ggplot2 figure match.
#'     ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_rbiom()
#'
#'     # 7 pt text for a single journal column (~3.5 inches wide).
#'     adiv_corrplot(biom, x = "Age", layers = "tp", base_size = 7)
#'

theme_rbiom <- function (
    base_size = 11, base_family = "", frame = "axes", grid = FALSE,
    base_line_size = base_size / 22, base_rect_size = base_size / 22 ) {

  if (!(is.numeric(base_size) && length(base_size) == 1 && isTRUE(base_size > 0)))
    cli_abort("`base_size` must be a positive number, not {.val {base_size}}.")

  if (!(is_scalar_character(base_family) && !is.na(base_family)))
    cli_abort("`base_family` must be a font family name, not {.type {base_family}}.")

  if (!is_string(frame, c("axes", "box", "none")))
    cli_abort("`frame` must be one of 'axes', 'box', or 'none', not {.val {frame}}.")

  if (!(isTRUE(grid) || isFALSE(grid) || is_string(grid, c("x", "y"))))
    cli_abort("`grid` must be TRUE, FALSE, 'x', or 'y', not {.val {grid}}.")


  half <- base_size / 2
  ink  <- "black"
  soft <- "grey30"   # subtitles and captions
  rule <- "grey88"   # gridlines

  text  <- ggplot2::element_text
  line  <- ggplot2::element_line
  rect  <- ggplot2::element_rect
  blank <- ggplot2::element_blank
  rel   <- ggplot2::rel

  gridline <- line(colour = rule, linewidth = rel(0.5))


  t <- ggplot2::theme_bw(
    base_size      = base_size,
    base_family    = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size ) %+replace% ggplot2::theme(

    # Black text; secondary text in gray.
    axis.text    = text(size = rel(0.85), colour = ink),
    axis.title   = text(colour = ink),
    legend.title = text(hjust = 0, colour = ink),
    legend.text  = text(size = rel(0.85), colour = ink),
    strip.text   = text(
      colour = ink,
      hjust  = 0,
      margin = ggplot2::margin(half / 2, half / 2, half / 2, 0) ),

    # Titles and captions start at the figure's left edge.
    plot.title            = text(size = rel(1.15), face = "bold", hjust = 0, vjust = 1, margin = ggplot2::margin(b = half)),
    plot.subtitle         = text(size = rel(0.9), colour = soft, hjust = 0, vjust = 1, margin = ggplot2::margin(b = half)),
    plot.caption          = text(size = rel(0.8), colour = soft, hjust = 0, vjust = 1, margin = ggplot2::margin(t = half)),
    plot.tag              = text(size = rel(1.2), face = "bold", hjust = 0.5, vjust = 0.5),
    plot.title.position   = "plot",
    plot.caption.position = "plot",

    # Open axes, outward ticks, no shading or gridlines. (The gridlines
    # are blanked individually so that `grid` can restore them.)
    axis.line          = line(colour = ink, lineend = "square"),
    axis.ticks         = line(colour = ink),
    panel.border       = blank(),
    panel.background   = blank(),
    panel.grid         = gridline,
    panel.grid.major.x = blank(),
    panel.grid.major.y = blank(),
    panel.grid.minor   = blank(),
    panel.spacing      = grid::unit(base_size, "pt"),
    strip.background   = blank(),
    plot.background    = rect(fill = "white", colour = NA),

    # Room on the right for the last x-axis tick label to overhang.
    plot.margin = ggplot2::margin(half, base_size, half, half),

    # Legends: no key backgrounds, aligned with the top of the panels.
    legend.background    = blank(),
    legend.key           = blank(),
    legend.key.size      = grid::unit(1.5 * base_size, "pt"),
    legend.justification = "top",

    complete = TRUE )


  if (frame == "box")
    t <- t %+replace% ggplot2::theme(
      axis.line    = blank(),
      panel.border = rect(fill = NA, colour = ink),
      complete     = TRUE )

  if (frame == "none")
    t <- t %+replace% ggplot2::theme(
      axis.line  = blank(),
      axis.ticks = blank(),
      axis.title = blank(),
      complete   = TRUE )

  if (isTRUE(grid) || eq(grid, "x"))
    t <- t %+replace% ggplot2::theme(panel.grid.major.x = line(), complete = TRUE)

  if (isTRUE(grid) || eq(grid, "y"))
    t <- t %+replace% ggplot2::theme(panel.grid.major.y = line(), complete = TRUE)

  return (t)
}



#________________________________________________________
# Font size (in mm, as for geom_text) for annotations such as p-values
# and taxa labels: 3 (8.5 pt) at the default `base_size` of 11, scaled
# with the `base_size` given to theme_rbiom() through a plot's `...`.
#________________________________________________________
annot_size <- function (params) {
  bs <- params$layers[['theme_rbiom']][['base_size']]
  if (!(is.numeric(bs) && length(bs) == 1 && isTRUE(bs > 0))) bs <- 11
  signif(3 * bs / 11, 3)
}
