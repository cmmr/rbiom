#' Create a heatmap with tracks and dendrograms from any matrix.
#' 
#' @name plot_heatmap
#' 
#' @family plotting
#' 
#' @param mtx   A numeric \code{matrix} with named rows and columns.
#' 
#' @param grid   Color palette name, or a list with entries for \code{label}, 
#'        \code{colors}, \code{range}, \code{bins}, \code{na.color}, and/or 
#'        \code{guide}. See the Track Definitions section for details.
#'        Default: \code{list(label = "Grid Value", colors = "imola")}.
#' 
#' @param tracks   List of track definitions. See details below.
#'        Default: \code{NULL}.
#'        
#' @param label   Label the matrix rows and columns. You can supply a list
#'        or logical vector of length two to control row labels and column 
#'        labels separately, for example 
#'        \code{label = c(rows = TRUE, cols = FALSE)}, or simply 
#'        \code{label = c(T, F)}. Other valid options are \code{"rows"},
#'        \code{"cols"}, \code{"both"}, \code{"bottom"}, \code{"right"},
#'        and \code{"none"}.
#'        Default: \code{TRUE}.
#'        
#' @param label_size   The font size to use for the row and column labels. You 
#'        can supply a numeric vector of length two to control row label sizes 
#'        and column label sizes separately, for example 
#'        \code{c(rows = 20, cols = 8)}, or simply \code{c(20, 8)}.
#'        Default: \code{NULL}, which computes: 
#'        \code{pmax(8, pmin(20, 100 / dim(mtx)))}.
#'        
#' @param rescale   Rescale rows or columns to all have a common min/max.
#'        Options: \code{"none"}, \code{"rows"}, or \code{"cols"}.
#'        Default: \code{"none"}.
#'        
#' @param trees  Draw a dendrogram for rows (left) and columns (top). You can 
#'        supply a list or logical vector of length two to control the row tree 
#'        and column tree separately, for example 
#'        \code{trees = c(rows = T, cols = F)}, or simply \code{trees = c(T, F)}. 
#'        Other valid options are \code{"rows"}, \code{"cols"}, \code{"both"}, 
#'        \code{"left"}, \code{"top"}, and \code{"none"}.
#'        Default: \code{TRUE}.
#'        
#' @param clust   Clustering algorithm for reordering the rows and columns by 
#'        similarity. You can supply a list or character vector of length two to 
#'        control the row and column clustering separately, for example 
#'        \code{clust = c(rows = "complete", cols = NA)}, or simply 
#'        \code{clust = c("complete", NA)}.
#'        Default: \code{"complete"}.
#'        
#'        Options are:
#'        \itemize{
#'          \item{\code{FALSE} or \code{NA} - }{ Disable reordering. }
#'          \item{An \link[stats]{hclust} object} { }
#'          \item{A method name - }{ \code{"ward.D"}, 
#'            \code{"ward.D2"}, \code{"single"}, \code{"complete"}, 
#'            \code{"average"}, \code{"mcquitty"}, \code{"median"}, or 
#'            \code{"centroid"}. }
#'        }
#'        
#' @param dist   Distance algorithm to use when reordering the rows and columns 
#'        by similarity. You can supply a list or character vector of length
#'        two to control the row and column clustering separately, for example 
#'        \code{dist = c(rows = "euclidean", cols = "maximum")}, or simply 
#'        \code{dist = c("euclidean", "maximum")}.
#'        Default: \code{"euclidean"}.
#'        
#'        Options are:
#'        \itemize{
#'          \item{A \link[stats]{dist} object} { }
#'          \item{A method name - }{ \code{"euclidean"}, 
#'            \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, 
#'            \code{"binary"}, or \code{"minkowski"}. }
#'        }
#'        
#' @param tree_height,track_height   The height of the dendrogram or annotation
#'        tracks in multiples (or fractions) of the smaller dimension of the
#'        grid cell size. Use a numeric vector of length two to assign
#'        \code{c(left, top)} independently. 
#'        Default: \code{NULL}, which computes:
#'        \code{tree_height = sqrt(min(dim(mtx))), track_height = tree_height / 4}.
#'        
#' @param ratio   Height/width ratio for entire grid. 
#'        Default: \code{1} (square).
#'        
#' @param legend   Where to place the legend. Options are: \code{"right"} or
#'        \code{"bottom"}. Default: \code{"right"}.
#' 
#' @param xlab.angle   Angle of the labels at the bottom of the plot. 
#'        Options are \code{"auto"}, \code{0}, \code{30}, and \code{90}. 
#'        Default: \code{"auto"}.
#'        
#' @param ...   Additional arguments to pass on to ggplot2::theme().
#'        For example, \code{labs.title = "Plot Title"}.
#'        
#' @return A \code{ggplot2} plot. The constructed ggplot command will be
#'         attached as \code{attr(,'cmd')}, and the underlying computed data
#'         as \code{$data}.
#' 
#' @section Track Definitions:
#' 
#' One or more colored tracks can be placed on the left and/or top of the 
#' heatmap grid to visualize associated metadata values.
#'          
#' \preformatted{  ## Categorical ----------------------------
#' cat_vals = sample(c("Male", "Female"), 10, replace = TRUE)
#' tracks   = list('Sex' = cat_vals)
#' tracks   = list('Sex' = list('values' = cat_vals, 'colors' = "bright"))
#' tracks   = list('Sex' = list(
#'   'values' = cat_vals, 
#'   'colors' = c('Male' = "blue", 'Female' = "red")) )
#' 
#' ## Numeric --------------------------------
#' num_vals = sample(25:40, 10, replace = TRUE)
#' tracks   = list('Age' = num_vals)
#' tracks   = list('Age' = list('values' = num_vals, 'colors' = "greens"))
#' tracks   = list('Age' = list('values' = num_vals, 'range' = c(0,50)))
#' tracks   = list('Age' = list(
#'   'label'  = "Age (Years)",
#'   'values' = num_vals, 
#'   'colors' = c("azure", "darkblue", "darkorchid") ))
#' 
#' ## Multiple Tracks ------------------------
#' tracks = list('Sex' = cat_vals, 'Age' = num_vals)
#' tracks = list(
#'   list('label' = "Sex", values' = cat_vals, 'colors' = "bright"),
#'   list('label' = "Age", values' = num_vals, 'colors' = "greens") )
#'   
#' plot_heatmap(matrix(sample(1:50), ncol=10), tracks = tracks)
#' }
#' 
#' The following entries in the track definitions are understood: 
#' 
#' \itemize{
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
#' All built-in color palettes are colorblind-friendly.
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
    tree_height  = NULL, track_height = NULL, ratio=1, 
    legend = "right", xlab.angle = "auto", ... ) {
    
  theme_args <- list(...)
  
  with_cache("plot_heatmap", environment(), NULL, local({
    
    
    #________________________________________________________
    # Sanity Checks.
    #________________________________________________________
    stopifnot(
      (is_logical(label) && length(label) <= 2) ||
      is_string(label, c("rows", "cols", "both", "bottom", "right", "none")) )
    stopifnot(
      (is_logical(trees) && length(trees) <= 2) ||
      is_string(trees, c("rows", "cols", "both", "left", "top", "none")) )
    
    if (!any(is_na(clust), is_false(clust), is(clust, 'hclust')))
      clust %<>% validate_arg(arg = 'clust', n = 1:2, allow_na = TRUE)
    
    if (!is(dist, 'dist'))
      dist %<>% validate_arg(arg = 'dist', n = 1:2)
    
    stopifnot(is_null(label_size)   || (is_double(label_size)   && length(label_size)   <= 2) )
    stopifnot(is_null(tree_height)  || (is_double(tree_height)  && length(tree_height)  <= 2))
    stopifnot(is_null(track_height) || (is_double(track_height) && length(track_height) <= 2))
    stopifnot(is_scalar_double(ratio))
    stopifnot(is_string(rescale, c("none", "rows", "cols")))
    stopifnot(is_string(legend, c("right", "bottom")))
    stopifnot(is_string(as.character(xlab.angle), c("auto", "0", "30", "90")))
    
    
    #________________________________________________________
    # Setting legend = "bottom" tiggers a bunch of presets
    #________________________________________________________
    if (identical(legend, "bottom")) {
      legend <- list(nrow = 2, title.position = "top")
      
      theme_args <- within(theme_args, {
        legend.box        %||=% "horizontal"
        legend.position   %||=% "bottom"
        legend.background %||=% element_rect(color = "black")
        legend.margin     %||=% margin(6,12,6,12,"points")
        legend.box.margin %||=% margin(6,6,6,6,"points")
      })
    }
    
    
    #________________________________________________________
    # Allow different specs for rows and cols.
    # For example, dist = c("euclidean", "canberra")
    #________________________________________________________
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
      if (identical(track[['side']], "top")) { top_tracks  %<>% c(list(track))
      } else                                 { left_tracks %<>% c(list(track)) }
      
    }
    n_top  <- length(top_tracks)
    n_left <- length(left_tracks)
    
    
    
    #________________________________________________________
    # Heuristics for setting sizes and heights of elements
    #________________________________________________________
    if (is_null(tree_height))  tree_height  <- min(dim(mtx)) * 0.15
    if (is_null(track_height)) track_height <- tree_height / 4
    
    ratio <- ratio * ncol(mtx) / nrow(mtx)
    
    if (is_null(label_size))         label_size <- pmax(5, pmin(14, 100 / dim(mtx)))
    if (!is_null(names(label_size))) label_size <- label_size[rev(order(names(label_size)))]
    label_size_x <- tail(label_size, 1)
    label_size_y <- head(label_size, 1)
    if (isTRUE(n_left > 0)) label_size_x <- .c(.rep(10, n_left), .rep(label_size_x, ncol(mtx)))
    if (isTRUE(n_top  > 0)) label_size_y <- .c(.rep(10, n_top),  .rep(label_size_y, nrow(mtx)))
    label_size_x <- unit(label_size_x, "points")
    label_size_y <- unit(label_size_y, "points")
    
    
    
    
    #________________________________________________________
    # Use specified height(s) + ratio to compute top & left heights
    #________________________________________________________
    stopifnot(is.numeric(tree_height) && is.numeric(track_height) && is.numeric(ratio))
    stopifnot(length(tree_height) > 0 && length(track_height) > 0 && length(ratio) == 1)
    tree_height_top   <- tail(tree_height,  1) * ifelse(ratio < 1, ratio, 1)
    tree_height_left  <- head(tree_height,  1) * ifelse(ratio > 1, ratio, 1)
    track_height_top  <- tail(track_height, 1) * ifelse(ratio < 1, ratio, 1)
    track_height_left <- head(track_height, 1) * ifelse(ratio > 1, ratio, 1)
    remove("tree_height", "track_height")
    
    
    #________________________________________________________
    # Rescale matrix's row or column values
    #________________________________________________________
    if (identical(rescale, "rows")) {
      mtx <- t(apply(mtx, 1L, scales::rescale))
      
    } else if (identical(rescale, "cols")) {
      mtx <- apply(mtx, 2L, scales::rescale)
    }
    
    
    #________________________________________________________
    # Cluster matrix's rows/cols and generate data.frames for dendrograms
    #________________________________________________________
    if (!isFALSE(clust_cols) && !is.na(clust_cols)) {
      
      dm <- dist_cols
      hc <- clust_cols
      if (!is(dm, 'dist') && !is(hc, 'hclust')) dm <- stats::dist(t(mtx), method = dm)
      if (!is(hc, 'hclust'))                    hc <- stats::hclust(dm,   method = hc)
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
      if (!is(dm, 'dist') && !is(hc, 'hclust')) dm <- stats::dist(mtx,  method = dm)
      if (!is(hc, 'hclust'))                    hc <- stats::hclust(dm, method = hc)
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
      fill = as.vector(mtx) )
    
    
    
    #________________________________________________________
    # Assemble the layers for the ggplot
    #________________________________________________________
    gglayers <- list()
    gglayers %<>% ggpush(ggplot(df))
    gglayers %<>% ggpush(coord_fixed(ratio=ratio))
    
    
    #________________________________________________________
    # In order to label any tracks, switch to continuous
    #________________________________________________________
    if (length(all_tracks) == 0) {
      gglayers %<>% ggpush(geom_raster(aes(x=x, y=y, fill=fill)))
      
      if (isFALSE(label_cols))
        gglayers %<>% ggpush(scale_x_discrete(breaks = NULL))
      
      if (isTRUE(label_rows)) {
        gglayers %<>% ggpush(scale_y_discrete(position = 'right'))
      } else {
        gglayers %<>% ggpush(scale_y_discrete(breaks = NULL))
      }
      
    } else {
      gglayers %<>% ggpush(geom_raster(aes(x=as.numeric(x), y=as.numeric(y), fill=fill)))
      
      
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
            .Data   = c(xbreaks, 1:length(levels(df[['x']]))),
            display = sprintf( fmt = "c(%s, 1:%i)", 
              paste(collapse = ", ", xbreaks), length(levels(df[['x']])) ))
        }
        
      } else if (isTRUE(label_cols)) {
        
        xlabels <- structure(
          .Data   = levels(df[['x']]),
          display = "levels(data[['x']])" )
        
        xbreaks <- structure(
          .Data   = 1:length(levels(df[['x']])),
          display = sprintf("1:%i", length(levels(df[['x']]))) )
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
          .Data   = 1:length(levels(df[['x']])),
          display = sprintf("1:%i", length(levels(df[['x']]))) )
      }
      
      
      
      #________________________________________________________
      # Outline around the main grid
      #________________________________________________________
      gglayers %<>% ggpush(geom_rect(
        mapping = aes(
          xmin = 0.5, 
          xmax = !!(ncol(mtx) + 0.5), 
          ymin = 0.5, 
          ymax = !!(nrow(mtx) + 0.5) ), 
        size  = 0.2,
        color = "black", 
        fill  = NA ))
      
      
      
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
      
      attr(gglayers[[1]][['data']], 'trees') <- all_trees
      
      gglayers %<>% ggpush(geom_segment( 
        data    = ~ attr(., 'trees', exact = TRUE),
        mapping = aes(x=x, y=y, xend=xend, yend=yend),
        .indent = 4 ))
    }
    
    
    for (i in seq_along(all_tracks)) {
      
      track <- all_tracks[[i]]
      layer <- plot_heatmap_layer(track, i)
      
      attr(gglayers[[1]][['data']], track[['id']]) <- track[['data']]
      data <- as.formula(sprintf("~attr(., '%s', exact = TRUE)", track[['id']]))
      
      gglayers %<>% ggpush(new_scale_fill())
      gglayers %<>% ggpush(geom_tile(data = data, mapping = track[['mapping']], .indent = 4))
      gglayers %<>% ggpush(do.call(layer[['fun']], c(layer[['args']], .indent = 4)))
      gglayers %<>% ggpush(geom_rect(mapping = track[['outline']], color = "black", fill = NA, size = 0.2 ))
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
    
    return (p)
    
  }))
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
  if (!is_null(args[['label']]))      fn_args[['name']]       <- args[['label']]
  if (!is_null(args[['range']]))      fn_args[['limits']]     <- args[['range']]
  if (!is_null(args[['bins']]))       fn_args[['n.breaks']]   <- args[['bins']] + 1
  if (!is_null(args[['na.color']]))   fn_args[['na.value']]   <- args[['na.color']]
  
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


# #________________________________________________________
# # Choose a suitable palette for a set of values.
# #________________________________________________________
# palette_autoselect <- function (vals, i) {
#   
#   if (is.numeric(vals)) { # Continuous variable
#     
#     opts <- list(
#       lo = paste0("grDevices::", c("Oranges", "Purples", "Greens", "Reds", "Blues", "Grays")),
#       hi = paste0("grDevices::", c("Oranges", "Purples", "Greens", "Reds", "Blues", "Grays")),
#       bi = paste0("khroma::",    c("cork", "vik", "bam")),
#       un = paste0("viridis::",   c("turbo", "viridis", "magma", "plasma", "inferno", "cividis", "mako", "bilbao")) )
#     
#     n <- length(vals)
#     x <- scales::rescale(vals, to = c(-1, 1))
#     
#     opt <- if (sum(x > 0) < n * 2/3)        { "lo"
#     } else if (sum(x > 0) > n * 2/3)        { "hi"
#     } else if (sum(abs(x) > 0.5) > n * 2/3) { "bi"
#     } else {                                  "un" }
#     
#     pal <- opts[[opt]][[(i - 1) %% length(opts[[opt]]) + 1]]
#     pal <- palette_colorcodes(pal)
#     if (opt == "lo") pal <- rev(pal)
#     
#   } else { # Categorical variable
#     
#     opts <- paste0("khroma::", c("muted", "light", "bright", "pale", "dark"))
#     pal  <- opts[[(i - 1) %% length(opts) + 1]]
#   }
#   
#   return (pal)
# }


# #________________________________________________________
# # Convert user-provided palette name to vector of colors.
# #________________________________________________________
# palette_colorcodes <- function (pal, n = 10, direction = 1) {
#   
#   if (length(pal) != 1) return (pal) # Common passthrough case
#   pal <- as.vector(pal)
#   
#   cols <- c("package", "palette", "type", "mode")
#   opts <- rbind(
#     data.frame(paletteer::palettes_c_names, mode = "c")[,cols], 
#     data.frame(paletteer::palettes_d_names, mode = "d")[,cols], 
#     data.frame(paletteer::palettes_dynamic_names, mode = "y")[,cols] )
#   opts[['id']]   <- paste(sep="::", opts[['package']], opts[['palette']])
#   opts[['type']] <- c(s="seq", d="div", q="qual")[substr(opts[['type']], 1, 1)]
#   
#   # Favor matching to RColorBrewer, khroma, or viridis
#   opts <- opts[order(!opts[['package']] %in% c("RColorBrewer", "khroma", "viridis")),]
#   
#   # Reverse color order when palette name is prefixed by '-'
#   if (identical(substr(pal, 1, 1), "-")) {
#     direction <- -1
#     pal <- substr(pal, 2, nchar(pal))
#   }
#   
#   # Match either "RColorBrewer::Oranges" or "oranges"
#   pal   <- sub("okabeito", "okabe_ito", tolower(pal))
#   parts <- strsplit(pal, "::", fixed = TRUE)[[1]]
#   row   <- if (length(parts) > 1) {
#     which(
#       tolower(opts[['package']]) == parts[[1]] & 
#       tolower(opts[['palette']]) == parts[[2]] )
#   } else {
#     which(tolower(opts[['palette']]) == parts[[1]])
#   }
#   
#   if (length(row) == 0)
#     stop("There is no paletteer palette named ", pal)
#   
#   # Pull the hex color values out of paletteer
#   pal <- as.list(opts[head(row, 1),])
#   if (pal[['mode']] == "c") {
#     res <- paletteer::paletteer_c(pal[['id']], n)
#     
#   } else if (pal[['mode']] == "y") {
#     res <- paletteer::paletteer_dynamic(pal[['id']], n)
#     
#   } else if (pal[['mode']] == "d") {
#     res <- paletteer::paletteer_d(pal[['id']])
#     
#     
#     if (identical(pal[['id']], "khroma::okabe_ito"))
#       res <- res[c(2:8, 1)] # put black last
#     
#     
#     if (!missing(n)) {
#       
#       if (n > length(res)) {
#         res <- rep_len(res, n)
#         
#       } else if (n < length(res)) {
#         if (pal[['type']] == "qual") {
#           res <- res[1:n]
#         } else {
#           res <- res[as.integer(seq(from=1, to=length(res), length.out=n))]
#         }
#       }
#       
#     }
#   }
#   
#   if (identical(direction, -1))
#     res <- rev(res)
#   
#   return (substr(as.vector(res), 1, 7))
# }



#________________________________________________________
# Create the dendrogram based on hclust output
#________________________________________________________
dendro <- function (hc, bounds=c(0, 1), side = "top") {
  
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
  fn <- function (i, prev_ht = NULL) {
    
    if (i < 0) {
      x <- which(hc[['order']] == abs(i))
      return (data.frame(x = x, y = prev_ht, xend = x, yend = bounds[[1]]))
    }
    
    ht   <- hc[['height']][[i]]
    df1  <- fn(hc[['merge']][i,1], ht)
    df2  <- fn(hc[['merge']][i,2], ht)
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
  
  df <- fn(nrow(hc[['merge']]))
  
  if (side %in% c("bottom", "right")) {
    df[['yend']] <- 1 - df[['yend']]
    df[['y']]    <- 1 - df[['y']]
  }
  if (side %in% c("left", "right"))
    df <- with(df, data.frame(x=y, y=x, xend=yend, yend=xend))
  
  return (df)
}



tracks_df <- function (tracks, bounds=c(0,1), side="top") {
  
  if (length(tracks) == 0) return (tracks)
  stopifnot(is_string(side, c("top", "left")))
  
  
  #________________________________________________________
  # compute the center position of the short edge for each track
  #________________________________________________________
  bounds_w <- abs(diff(bounds)) / length(tracks)
  bounds   <- rev(seq_along(tracks) - 0.5) * bounds_w + min(bounds)
  
  
  for (i in seq_along(tracks)) {
    
    values <- tracks[[i]][['values']]
    sf     <- floor(log10(length(values))) + 3 # Sig. fig. digits
    
    if (isTRUE(side == "top")) {
      
      tracks[[i]][['data']] <- data.frame(
        fill = values, 
        x    = seq_along(values) )
      
      tracks[[i]][['mapping']] <- aes(
        x      = x,
        y      = !!signif(bounds[i], sf),
        fill   = fill,
        height = !!signif(bounds_w, sf),
        width  = 1 )
      
      tracks[[i]][['outline']] <- aes(
        xmin = 0.5,
        xmax = !!(length(values) + 0.5),
        ymin = !!signif(bounds[i] - bounds_w / 2, sf),
        ymax = !!signif(bounds[i] + bounds_w / 2, sf) )
      
    } else {
      
      tracks[[i]][['data']] <- data.frame(
        fill = values, 
        y    = seq_along(values) )
      
      tracks[[i]][['mapping']] <- aes(
        x      = !!signif(bounds[i], sf),
        y      = y,
        fill   = fill,
        height = 1,
        width  = !!signif(bounds_w, sf) )
      
      tracks[[i]][['outline']] <- aes(
        xmin = !!signif(bounds[i] - bounds_w / 2, sf),
        xmax = !!signif(bounds[i] + bounds_w / 2, sf),
        ymin = 0.5,
        ymax = !!(length(values) + 0.5) )
    }
    
    tracks[[i]][['label_at']] <- signif(bounds[i], sf)
    tracks[[i]][['id']]       <- paste0(side, "_track_", i)
  }
  
  return (tracks)
}
