
#________________________________________________________
# Helpers for plot_heatmap(): a list of provenance-tracked ggplot2
# layers that is combined into a plot with `+`, dendrogram segments
# from an hclust object, and the geometry of the annotation tracks.
#________________________________________________________



#________________________________________________________
# Add a layer to a list of layers.
#________________________________________________________

ggpush <- function (gglayers, gglayer) {
  gglayers[[length(gglayers) + 1]] <- gglayer
  return (gglayers)
}


#________________________________________________________
# Combine a list of logged commands into a plot.
#________________________________________________________

ggbuild <- function (gglayers) {

  p   <- NULL
  cmd <- NULL


  for (i in seq_along(gglayers)) {

    gglayer <- gglayers[[i]]

    # In case this layer was built by init_layers / set_layer
    if (!is_null(fun <- attr(gglayer, 'function', exact = TRUE)))
      gglayer <- do.call(fun, c(gglayer, '.indent' = 4))

    if (is_null(p)) {
      p   <- gglayer
      cmd <- attr(gglayer, 'display')
    } else {
      p   <- p + gglayer
      cmd <- sprintf("%s +\n  %s", cmd, attr(gglayer, 'display'))
    }
  }

  attr(p, 'display') <- NULL
  attr(p, 'code')    <- add_class(cmd, 'rbiom_code')

  p$plot_env <- emptyenv()

  return (p)
}



#________________________________________________________
# Create the dendrogram based on hclust output
#________________________________________________________
dendro <- function (hc, bounds=c(0, 1), side = "top") {

  side <- match.arg(side, choices = c("top", "right", "bottom", "left"))


  #________________________________________________________
  # Allow user to control absolute positioning
  #________________________________________________________
  hc[['height']] <- scales::rescale(
    x    = hc[['height']],
    to   = bounds,
    from = c(0, max(hc[['height']])))


  #________________________________________________________
  # geom_segments
  #________________________________________________________
  func <- function (i, prev_ht = NULL) {

    if (i < 0) {
      x <- which(hc[['order']] == abs(i))
      return (data.frame(x = x, y = prev_ht, xend = x, yend = bounds[[1]]))
    }

    ht   <- hc[['height']][[i]]
    df1  <- func(hc[['merge']][i,1], ht)
    df2  <- func(hc[['merge']][i,2], ht)
    x    <- (df1[1,'x'] + df1[1,'xend']) / 2
    xend <- (df2[1,'x'] + df2[1,'xend']) / 2

    df <- if (is_null(prev_ht)) {
      data.frame(x = x, y = ht, xend = xend, yend = ht)
    } else {
      center <- (x + xend) / 2
      data.frame(
        x    = c(x, center),
        y    = c(ht, prev_ht),
        xend = c(xend, center),
        yend = c(ht, ht))
    }

    return (rbind(df, rbind(df1, df2)))
  }

  # A single leaf has no merges and no dendrogram.
  if (nrow(hc[['merge']]) == 0)
    return (data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric()))

  df <- func(nrow(hc[['merge']]))

  if (side %in% c("bottom", "right")) {
    df[['yend']] <- 1 - df[['yend']]
    df[['y']]    <- 1 - df[['y']]
  }
  if (side %in% c("left", "right"))
    df <- with(df, data.frame(x=y, y=x, xend=yend, yend=xend))

  rownames(df) <- NULL

  return (df)
}



tracks_df <- function (tracks, bounds=c(0,1), side="top") {



  if (length(tracks) == 0) return (tracks)

  side <- match.arg(side, choices = c("top", "left"))


  #________________________________________________________
  # compute the center position of the short edge for each track
  #________________________________________________________
  bounds_w <- abs(diff(bounds)) / length(tracks)
  bounds   <- rev(seq_along(tracks) - 0.5) * bounds_w + min(bounds)


  x <- y <- fill <- NULL # for CRAN check only

  for (i in seq_along(tracks)) {

    values <- tracks[[i]][['values']]
    sf     <- floor(log10(length(values))) + 3 # Sig. fig. digits

    if (isTRUE(side == "top")) {

      tracks[[i]][['data']] <- tibble(
        fill = values,
        x    = seq_along(values) )

      tracks[[i]][['mapping']] <- aes(
        x      = x,
        y      = !!signif(bounds[i], sf),
        fill   = fill,
        height = !!signif(bounds_w, sf),
        width  = 1 )

      tracks[[i]][['outline']] <- list(
        xmin = 0.5,
        xmax = length(values) + 0.5,
        ymin = signif(bounds[i] - bounds_w / 2, sf),
        ymax = signif(bounds[i] + bounds_w / 2, sf) )

    } else {

      tracks[[i]][['data']] <- tibble(
        fill = values,
        y    = seq_along(values) )

      tracks[[i]][['mapping']] <- aes(
        x      = !!signif(bounds[i], sf),
        y      = y,
        fill   = fill,
        height = 1,
        width  = !!signif(bounds_w, sf) )

      tracks[[i]][['outline']] <- list(
        xmin = signif(bounds[i] - bounds_w / 2, sf),
        xmax = signif(bounds[i] + bounds_w / 2, sf),
        ymin = 0.5,
        ymax = length(values) + 0.5 )
    }

    tracks[[i]][['label_at']] <- signif(bounds[i], sf)
    tracks[[i]][['id']]       <- paste0(side, "_track_", i)
  }

  return (tracks)
}
