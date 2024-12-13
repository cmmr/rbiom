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
#' All built-in color palettes are colorblind-friendly. See 
#' [Mapping Metadata to Aesthetics](https://cmmr.github.io/rbiom/articles/aes.html#discrete-palettes) 
#' for images of the palettes.
#' 
#' Categorical palette names: \code{"okabe"}, \code{"carto"}, \code{"r4"}, 
#' \code{"polychrome"}, \code{"tol"}, \code{"bright"}, \code{"light"}, 
#' \code{"muted"}, \code{"vibrant"}, \code{"tableau"}, \code{"classic"}, 
#' \code{"alphabet"}, \code{"tableau20"}, \code{"kelly"}, and \code{"fishy"}.
#' 
#' Numeric palette names: \code{"reds"}, \code{"oranges"}, \code{"greens"}, 
#' \code{"purples"}, \code{"grays"}, \code{"acton"}, \code{"bamako"}, 
#' \code{"batlow"}, \code{"bilbao"}, \code{"buda"}, \code{"davos"}, 
#' \code{"devon"}, \code{"grayC"}, \code{"hawaii"}, \code{"imola"}, 
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
#'     set.seed(123)
#'     mtx <- matrix(runif(5*8), nrow = 5, dimnames = list(LETTERS[1:5], letters[1:8]))
#'     
#'     plot_heatmap(mtx)
#'     plot_heatmap(mtx, grid="oranges")
#'     plot_heatmap(mtx, grid=list(colors = "oranges", label = "Some %", bins = 5))
#'     
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
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('plot_heatmap', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
    
  theme_args <- params$.dots
  
  
  #________________________________________________________
  # Sanity Checks.
  #________________________________________________________
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
  
  stopifnot(is_null(label_size)   || (is_double(label_size)   && length(label_size)   <= 2) )
  stopifnot(is_null(tree_height)  || (is_double(tree_height)  && length(tree_height)  <= 2))
  stopifnot(is_null(track_height) || (is_double(track_height) && length(track_height) <= 2))
  stopifnot(is_scalar_double(asp) && !is_na(asp))
  stopifnot(is_string(rescale, c("none", "rows", "cols")))
  stopifnot(is_string(legend, c("right", "bottom")))
  stopifnot(is_string(as.character(xlab.angle), c("auto", "0", "30", "90")))
  
  
  #________________________________________________________
  # Setting legend = "bottom" tiggers a bunch of presets
  #________________________________________________________
  if (eq(legend, "bottom")) {
    legend <- list(nrow = 2, title.position = "top")
    theme_args[['legend.box']]        %<>% if.null("horizontal")
    theme_args[['legend.position']]   %<>% if.null("bottom")
    theme_args[['legend.background']] %<>% if.null(element_rect(color = "black"))
    theme_args[['legend.margin']]     %<>% if.null(margin(6,12,6,12,"points"))
    theme_args[['legend.box.margin']] %<>% if.null(margin(6,6,6,6,"points"))
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
    
    if (length(val) == 2) {
      
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
  # Combine heights + aspect ratio to get top/left heights
  #________________________________________________________
  stopifnot(is.numeric(tree_height) && is.numeric(track_height) && is.numeric(asp))
  stopifnot(length(tree_height) > 0 && length(track_height) > 0 && length(asp) == 1)
  
  asp <- asp * ncol(mtx) / nrow(mtx)
  
  tree_height_top   <- tree_height[[1]]  / 100 * nrow(mtx)
  track_height_top  <- track_height[[1]] / 100 * nrow(mtx)
  tree_height_left  <- if (length(tree_height_top) == 1) tree_height_top  * asp else tree_height[[2]]  / 100 * ncol(mtx)
  track_height_left <- if (length(track_height)    == 1) track_height_top * asp else track_height[[2]] / 100 * ncol(mtx)
  remove("tree_height", "track_height")
  
  
  
  
  #________________________________________________________
  # Heuristics for setting sizes and heights of elements
  #________________________________________________________
  if (is_null(label_size))         label_size <- pmax(5, pmin(14, 100 / dim(mtx)))
  if (!is_null(names(label_size))) label_size <- label_size[rev(order(names(label_size)))]
  label_size_x <- tail(label_size, 1) %>% signif(digits = 3)
  label_size_y <- head(label_size, 1)
  
  if (is.null(tracks) || length(tracks) == 0) {
    n_left <- n_top <- 0
  } else {
    n_left <- sum(sapply(tracks, function (x) is.list(x) && eq(x[['side']], 'left')))
    n_top  <- length(tracks) - n_left
  }
  if (isTRUE(n_left > 0)) label_size_x <- .c(.rep(10, n_left), .rep(label_size_x, ncol(mtx)))
  if (isTRUE(n_top  > 0)) label_size_y <- .c(.rep(10, n_top),  .rep(label_size_y, nrow(mtx)))
  label_size_x <- unit(label_size_x, "points")
  label_size_y <- unit(label_size_y, "points")
  
  
  
  #________________________________________________________
  # Rescale matrix's row or column values
  #________________________________________________________
  if (eq(rescale, "rows")) {
    mtx <- t(apply(mtx, 1L, scales::rescale))
    
  } else if (eq(rescale, "cols")) {
    mtx <- apply(mtx, 2L, scales::rescale)
  }
  
  
  #________________________________________________________
  # Cluster matrix's rows/cols and generate data.frames for dendrograms
  #________________________________________________________
  if (!isFALSE(clust_cols) && !is.na(clust_cols)) {

    dm <- dist_cols
    hc <- clust_cols
    if (!inherits(dm, 'dist') && !inherits(hc, 'hclust')) dm <- stats::dist(t(mtx), method = dm)
    if (!inherits(hc, 'hclust'))                    hc <- stats::hclust(dm,   method = hc)
    mtx <- mtx[,hc[['order']]]

    if (isTRUE(trees_cols)) {

      bounds <- 0.5 + nrow(mtx) + c(0, tree_height_top)
      if (isTRUE(n_top > 0))
        bounds <- bounds + (nrow(mtx) / 50) + track_height_top * n_top

      trees_cols <- dendro(hc, bounds, side = "top")
    }
    remove("dm", "hc")

  } else {
    trees_cols <- FALSE
  }
  
  if (!isFALSE(clust_rows) && !is.na(clust_rows)) {

    dm <- dist_rows
    hc <- clust_rows
    if (!inherits(dm, 'dist') && !inherits(hc, 'hclust')) dm <- stats::dist(mtx,  method = dm)
    if (!inherits(hc, 'hclust'))                    hc <- stats::hclust(dm, method = hc)

    hc  <- rev(hc) # Move interesting taxa to top rows of plot.
    mtx <- mtx[hc[['order']],]

    if (isTRUE(trees_rows)) {

      bounds <- 0.5 - c(0, tree_height_left)
      if (isTRUE(n_left > 0))
        bounds <- bounds - (ncol(mtx) / 50) - track_height_top * n_left

      trees_rows <- dendro(hc, bounds, side = "left")
    }
    remove("dm", "hc")

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
    
    if (is_null(track[['values']])) stop("`values` needed for track ", track[['label']])
    
    
    # Ensure that values is the same length as rows or cols.
    #________________________________________________________
    track[['values']] <- local({
      
      values     <- track[['values']]
      side_len   <- if (track[['side']] == "top") ncol(mtx)     else nrow(mtx)
      side_names <- if (track[['side']] == "top") colnames(mtx) else rownames(mtx)
      
      if (is_null(names(values)) || is_null(side_names)) {
        if (!isTRUE(length(values) == side_len))
          stop("Length mismatch for ", track[['label']], ". Expected ", side_len, " values, got ", length(values))
        
      } else {
        if (length(missing <- setdiff(side_names, names(values))) > 0)
          stop("Missing ", length(missing), " row/col names: ", paste(collapse = ", ", head(missing)))
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
            stop("Missing ", length(missing), " color mappings: ", paste(collapse = ", ", head(missing)))
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
  
  df <- data.frame(
      x    = factor(rep(colnames(mtx), each = nrow(mtx)), levels=colnames(mtx)),
      y    = factor(rep(rownames(mtx),        ncol(mtx)), levels=rownames(mtx)),
      fill = as.vector(mtx) ) %>%
    as_rbiom_tbl()
  
  
  
  #________________________________________________________
  # Assemble the layers for the ggplot
  #________________________________________________________
  gglayers <- list()
  gglayers %<>% ggpush(ggplot(df))
  gglayers %<>% ggpush(coord_fixed(ratio=asp))
  
  if (!is.null(title))
    gglayers %<>% ggpush(labs(title=title))
  
  
  #________________________________________________________
  # In order to label any tracks, switch to continuous
  #________________________________________________________
  if (length(all_tracks) == 0) {
    
    x <- y <- fill <- NULL # for CRAN check only
    gglayers %<>% ggpush(geom_raster(aes(x=x, y=y, fill=fill)))
    remove("x", "y", "fill")
    
    if (isFALSE(label_cols))
      gglayers %<>% ggpush(scale_x_discrete(breaks = NULL))
    
    if (isTRUE(label_rows)) {
      gglayers %<>% ggpush(scale_y_discrete(position = 'right'))
    } else {
      gglayers %<>% ggpush(scale_y_discrete(breaks = NULL))
    }
    
  } else {
    
    x <- y <- fill <- NULL # for CRAN check only
    gglayers %<>% ggpush(geom_raster(aes(x=as.numeric(x), y=as.numeric(y), fill=fill)))
    remove("x", "y", "fill")
    
    
    #________________________________________________________
    # scale_x_continuous labels and breaks
    #________________________________________________________
    
    xlabels <- NULL
    xbreaks <- NULL
    
    if (isTRUE(n_left > 0)) {
      
      xlabels <- sapply(left_tracks, function (x) x[['label']])
      xbreaks <- sapply(left_tracks, function (x) x[['label_at']]) %>% 
        signif(floor(log10(nrow(mtx))) + 3)
      
      if (isTRUE(label_cols)) {
        
        xlabels <- structure(
          .Data   = c(xlabels, levels(df[['x']])),
          display = sprintf(fmt = "c(%s, levels(data[['x']]))", 
            paste(collapse = ", ", single_quote(xlabels)) ))
        
        xbreaks <- structure(
          .Data   = c(xbreaks, 1:nlevels(df[['x']])),
          display = sprintf(fmt = "c(%s, 1:%i)", 
            paste(collapse = ", ", xbreaks), nlevels(df[['x']]) ))
      }
      
    } else if (isTRUE(label_cols)) {
      
      xlabels <- structure(
        .Data   = levels(df[['x']]),
        display = "levels(data[['x']])" )
      
      xbreaks <- structure(
        .Data   = 1:nlevels(df[['x']]),
        display = sprintf("1:%i", nlevels(df[['x']])) )
    }
    
    
    #________________________________________________________
    # scale_y_continuous labels and breaks
    #________________________________________________________
    
    ylabels <- NULL
    ybreaks <- NULL
    
    if (isTRUE(n_top > 0)) {
      
      ylabels <- sapply(top_tracks, function (x) x[['label']])
      ybreaks <- sapply(top_tracks, function (x) x[['label_at']]) %>% 
        signif(floor(log10(ncol(mtx))) + 3)
      
      if (isTRUE(label_rows)) {
        
        ylabels <- structure(
          .Data   = c(ylabels, levels(df[['y']])),
          display = sprintf(fmt = "c(%s, levels(data[['y']]))", 
            paste(collapse = ", ", single_quote(ylabels)) ))
        
        ybreaks <- structure(
          .Data   = c(ybreaks, 1:length(levels(df[['y']]))),
          display = sprintf( fmt = "c(%s, 1:%i)", 
            paste(collapse = ", ", ybreaks), length(levels(df[['y']])) ))
      }
      
    } else if (isTRUE(label_rows)) {
      
      ylabels <- structure(
        .Data   = levels(df[['y']]),
        display = "levels(data[['y']])" )
      
      ybreaks <- structure(
        .Data   = 1:nlevels(df[['y']]),
        display = sprintf("1:%i", nlevels(df[['y']])) )
    }
    
    
    
    #________________________________________________________
    # Call scale_y_continuous and scale_y_continuous
    #________________________________________________________
    
    x_args <- list(expand = c(0.01,0), .indent = 4)
    y_args <- list(expand = c(0.01,0), .indent = 4)
    if (!is_null(xlabels)) x_args %<>% c(list(labels = xlabels, breaks = xbreaks))
    if (!is_null(ylabels)) y_args %<>% c(list(labels = ylabels, breaks = ybreaks, position = "right"))
    
    gglayers %<>% ggpush(do.call(scale_x_continuous, x_args))
    gglayers %<>% ggpush(do.call(scale_y_continuous, y_args))
    
    remove("xlabels", "xbreaks", "ylabels", "ybreaks", "x_args", "y_args")
    
  }
  
  
  
  #________________________________________________________
  # Outline around the main grid
  #________________________________________________________
  gglayers %<>% ggpush(annotate(
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
    
    attr(gglayers[[1]][['data']], 'trees') <- as_rbiom_tbl(all_trees)
    
    xend <- yend <- NULL # for CRAN check only
    
    gglayers %<>% ggpush(geom_segment( 
      data    = ~ attr(., 'trees', exact = TRUE),
      mapping = aes(x=x, y=y, xend=xend, yend=yend),
      .indent = 4 ))
  }
  
  
  for (i in seq_along(all_tracks)) {
    
    track <- all_tracks[[i]]
    layer <- plot_heatmap_layer(track, i)
    
    attr(gglayers[[1]][['data']], track[['id']]) <- as_rbiom_tbl(track[['data']])
    data <- as.formula(sprintf("~attr(., '%s', exact = TRUE)", track[['id']]))
    
    gglayers %<>% ggpush(new_scale_fill())
    gglayers %<>% ggpush(geom_tile(data = data, mapping = track[['mapping']], .indent = 4))
    gglayers %<>% ggpush(do.call(layer[['fun']], c(layer[['args']], .indent = 4)))
    gglayers %<>% ggpush(do.call(annotate, c(list(geom = "rect", color = "black", fill = NA, linewidth = 0.2), track[['outline']])))
  }
  
  
  #________________________________________________________
  # Theming
  #________________________________________________________
  args <- theme_args %||% list()
  args <- within(args, {
    
    if (!exists('axis.text.y') && (isTRUE(label_rows) || isTRUE(n_top > 0)))
      axis.text.y <- suppressWarnings(element_text(size=label_size_y, hjust=0))
    
    if (!exists('axis.text.x') && (isTRUE(label_cols) || isTRUE(n_left > 0)))
      axis.text.x <- suppressWarnings(switch(
        EXPR = paste0("x", tolower(xlab.angle)),
        x0   = element_text(size=label_size_x),
        x90  = element_text(size=label_size_x, angle=-90, vjust=0.3, hjust=0),
               element_text(size=label_size_x, angle=-30, vjust=1.0, hjust=0) ))
  })
  
  
  gglayers %<>% ggpush(theme_void())
  gglayers %<>% ggpush(do.call(theme, c(args, .indent = 4)))
  
  remove("args")
  p <- ggbuild(gglayers)
  p$plot_env <- emptyenv()
  p %<>% add_class('rbiom_plot')
  
  
  attr(p, 'cmd') <- current_cmd('plot_heatmap')
  set_cache_value(cache_file, p)
  
  return (p)
}




#________________________________________________________
# Convert palette args to scale_fill_* function
#________________________________________________________
plot_heatmap_layer <- function (args, guide_order = 0) {
  
  colors  <- args[['colors']]
  guide   <- as.list(args[['guide']])
  binned  <- !is_null(args[['bins']])
  numeric <- is_null(args[['values']]) || is.numeric(args[['values']])
  
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
    
    guide[['.indent']]        %<>% if.null(6)
    guide[['order']]          %<>% if.null(guide_order)
    guide[['direction']]      %<>% if.null('horizontal')
    guide[['title.position']] %<>% if.null('top')
    guide[['barheight']]      %<>% if.null(unit(21.5, "points"))
    guide[['label.theme']]    %<>% if.null(element_text(size = unit(8, "points")))
    fn_args[['guide']] <- do.call(guide_colorbar, guide)
  
  } else {
    
    # Create a discrete color scale
    fn <- "scale_fill_manual"
    fn_args[['values']] <- colors
    
    guide[['order']] %<>% if.null(guide_order)
    guide[['ncol']]  %<>% if.null(2)
    fn_args[['guide']] <- do.call(guide_legend, guide)
  }
  
  
  
  return (list(fn = fn, fun = eval(parse(text=fn)), args = fn_args))
}



#________________________________________________________
# Convert user's `tracks` into a named list of lists.
#________________________________________________________
biom_tracks <- function (params) {
  
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
    
    
    names(tracks)[[i]] <- key
    tracks[[i]]        <- val
    tracks[[i]]$values <- setNames(df[[key]], df$.sample)
  }
  
  params$tracks        <- tracks
  params$biom$metadata <- df
  
  
  return (invisible(params))
}
