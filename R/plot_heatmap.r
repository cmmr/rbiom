#' Create a heatmap with tracks and dendrograms from any matrix.
#' 
#' @name plot_heatmap
#' 
#' @family plotting
#' 
#' @param mtx   A numeric \code{matrix} with named rows and columns.
#'        
#' @param colors   A named \code{list()} with names matching \code{".grid"},
#'        \code{colnames(top_tracks)}, and \code{colnames(left_tracks)}. Values
#'        in this list should be a \bold{palette} definition (see section below).
#'        Default: \code{"rocket"}.
#'        
#'        For simple cases, a character vector can be given instead of a list. 
#'        For example, \code{colors = "viridis"} is understood to mean the 
#'        color palette name for the main grid, and 
#'        \code{colors = c(Age = "magma", Male = "blue", Female = "red")}
#'        will apply the magma palette to the Age metadata track, and blue and
#'        red to the Sex metadata track. The default "rocket" palette will be 
#'        used for the main grid unless otherwise defined.
#'        
#'        To set a custom title in the color scale legend, provide a key 
#'        prefixed by \code{.} like this:
#'        
#'        \code{colors = list(.Abundance = list(palette = "rocket", name = "Grid Value")}.
#'        
#'        To reverse the color order in a predefined palette, prepend a hyphen
#'        to the name, for example \code{"-rocket"}.
#'        
#' @param top_tracks   A \code{data.frame()} with row names matching 
#'        \code{colnames(mtx)}. Each column will be its own track at the top
#'        of the plot. Default: \code{NULL}.
#'        
#' @param left_tracks   A \code{data.frame()} with row names matching 
#'        \code{rownames(mtx)}. Each column will be its own track at the left
#'        of the plot. Default: \code{NULL}.
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
#'        Options: \code{"rows"}, \code{"cols"}, or \code{NULL}.
#'        Default: \code{NULL} (no rescaling).
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
#'          \item{An \link[stats]{hclust} method name - }{ \code{"ward.D"}, 
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
#'        Default: \code{dist = "euclidean"}.
#'        
#'        Options are:
#'        \itemize{
#'          \item{A \link[stats]{dist} object} { }
#'          \item{A \link[stats]{dist} method name - }{ \code{"euclidean"}, 
#'            \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, 
#'            \code{"binary"}, or \code{"minkowski"}. }
#'        }
#'        
#' @param tree_height,track_height   The height of the dendrogram or annotation
#'        tracks in multiples (or fractions) of the smaller dimension of the
#'        grid cell size. Use a numeric vector of length two to assign
#'        \code{c(left, top)} independently. Default: \code{NULL}, which computes:
#'        \code{tree_height = sqrt(min(dim(mtx))), track_height = tree_height / 4}.
#'        
#' @param ratio   Height/width ratio for entire grid.
#'        Default: \code{1} (square).
#'        
#' @param legend   Where to place the legend. Options are: \code{"right"} or
#'        \code{"bottom"}. Default: \code{"right"}.
#'        
#' @return A \code{ggplot2} plot. The constructed ggplot command will be
#'         attached as \code{attr(,'cmd')}, and the underlying computed data
#'         as \code{$data}.
#' 
#' @section Palettes:
#'          A palette can be as simple as the name of a color set, for example:
#'          
#'          \code{colors = "oranges"}
#'          
#'          You can use any name defined in the \code{paletteer} R package, 
#'          which aggregates color sets from a multitude of other R
#'          packages including \code{viridis}, \code{RColorBrewer}, and
#'          \code{ggpubr}, to name a few.
#'          
#'          Gradients from the \bold{viridis} and \bold{khroma} packages have
#'          good support for colorblindness:   \code{"viridis"}, \code{"magma"}, 
#'          \code{"plasma"}, \code{"inferno"}, \code{"cividis"}, \code{"mako"}, 
#'          \code{"rocket"}, \code{"turbo"},   \code{"broc"},    \code{"cork"},
#'          \code{"vik"},    \code{"lisbon"},  \code{"tofino"},  \code{"berlin"},
#'          \code{"roma"},   \code{"turku"},   \code{"vanimo"},  \code{"batlow"},
#'          \code{"devon"},  \code{"lajolla"}, \code{"bamako"},  \code{"davos"},
#'          \code{"bilbao"}, \code{"nuuk"},    \code{"oslo"},    \code{"grayC"},
#'          \code{"hawaii"}, \code{"lapaz"},   \code{"tokyo"},   \code{"buda"},
#'          \code{"acton"},  \code{"bam"},  or \code{"imola"}.
#'          
#'          To assign metadata track colors in addition to grid colors, combine
#'          them using a list. Any missing definitions will be auto-assigned.
#'          
#'          \code{colors = list("oranges", Sex = "wsj_red_green", 'Body Site' = "bugs")}
#'          
#'          You can also create custom color sets. For categorical data,
#'          provide a named character vector of color names. For continuous
#'          data, provide an unnamed character vector defining a gradient.
#'          
#'          \code{colors = list(Sex = list(palette = c('Male' = "blue", 'Female' = "#FF3390")))}
#'          
#'          \code{colors = list(Age = list(palette = c("azure", "darkblue", "darkorchid")))}
#'          
#'          You can also use this nested list structure to pass additional 
#'          parameters to the underlying ggplot2 \link[ggplot2]{discrete_scale} 
#'          or \link[ggplot2]{continuous_scale}, for example: 
#'          
#'          \code{colors = list(.grid = list(palette = "turbo", n.breaks = 5, 
#'          limits = c(0,1), reverse = TRUE, name = "Relative Abundance", 
#'          na.value = "grey50")}
#'          
#'          \itemize{
#'            \item{\code{name} - } { The title for this color scale in the legend. }
#'            \item{\code{na.value} - }{ The color to use for \code{NA} values. }
#'            \item{\code{limits} - }{ The c(min,max) to use for scale values. }
#'            \item{\code{n.breaks} - }{ Bin a gradient into this many bins/steps. }
#'            \item{\code{reverse} - }{ Reverse the order of colors. }
#'          }
#'          
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     set.seed(123)
#'     mtx <- matrix(runif(5*8), nrow = 5)
#'     rownames(mtx) <- sample(LETTERS, nrow(mtx))
#'     colnames(mtx) <- sample(letters, ncol(mtx))
#'     
#'     plot_heatmap(mtx)
#'     plot_heatmap(mtx, colors="oranges")
#'     plot_heatmap(mtx, colors=list(palette = "oranges", name = "New Label", n.breaks = 5))
#'     
#'     top_tracks <- data.frame(
#'       row.names = colnames(mtx), 
#'       'Bool'    = sample(c(TRUE, FALSE), ncol(mtx), TRUE), 
#'       'Int'     = sample(1:20,           ncol(mtx), TRUE) )
#'     left_tracks <- data.frame(
#'       row.names = rownames(mtx),
#'       'City'    = sample(names(precip), nrow(mtx)) )
#'     plot_heatmap(mtx, top_tracks=top_tracks, left_tracks=left_tracks)
#'     
plot_heatmap <- function (
    mtx, colors = "rocket",
    top_tracks = NULL, left_tracks = NULL,
    label = TRUE, label_size = NULL, rescale = NULL,
    trees = TRUE, clust = "complete", dist = "euclidean",
    tree_height  = NULL, track_height = NULL, ratio=1, 
    legend = NULL, colorbar = NULL, ... ) {
  
  theme_args <- list(...)
  
  
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
    
    colorbar %||=% list()
    colorbar <- within(colorbar, {
      barheight      %||=% unit(21.5, "points")
      title.position %||=% "top"
      label.theme    %||=% element_text(size = unit(8, "points"))
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
      
      if (!is.null(names(val)))
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
  # Split colors into grid_colors and track_colors.
  #________________________________________________________
  
  grid_colors  <- list()
  track_colors <- list()
  
  if (is.list(colors)) {
    
    # A proper list of lists
    if (all(sapply(colors, is.list))) {
      
      track_keys <- c(colnames(left_tracks), colnames(top_tracks))
      
      for (i in seq_along(colors)) {
        key <- names(colors)[[i]]
        
        if (is.null(key))               { grid_colors         <- colors[[i]] 
        } else if (key == "")           { grid_colors         <- colors[[i]]
        } else if (key %in% track_keys) { track_colors[[key]] <- colors[[i]] 
        } else if (substr(key, 1, 1) == ".") {
          grid_colors <- colors[[i]]
          if (!is.list(grid_colors))          grid_colors <- list(palette = grid_colors)
          if (is.null(grid_colors[['name']])) grid_colors[['name']] <- sub("^\\.", "", key)
        } else {
          warning(
            "In plot_heatmap(), '", key, 
            "' is in `colors`, but not `left_tracks` or `top_tracks`" )
        }
      }
      
    } else if (all(!sapply(colors, is.list))) {
      # A less proper unlisted list: color = list(palette = "oranges", name = "Vals")
      grid_colors <- colors
      
    } else {
      # A list of lists and non-lists
      stop("In plot_heatmap, `colors` contains mixed content.")
    }
    
    remove(list = intersect(ls(), c("i", "key")))
    
  } else if (is.character(colors)) {
    
    if (is.null(names(colors))) {
      
      # colors = "magma"
      # colors = c("blue", "white", "red")
      #________________________________________________________
      
      grid_colors[['palette']] <- colors
      
    } else {
      
      
      keys <- names(colors)
      
      # colors = c("viridis", ...)
      #________________________________________________________
      grid_colors[['palette']] <- colors[which(keys == "")]
      
      # colors = c(.Abundance = "viridis", ...)
      #________________________________________________________
      if (length(i <- which(substr(names(colors), 1, 1) == ".")) > 0) {
        grid_colors[['palette']] <- colors[head(i, 1)]
        grid_colors[['name']]    <- sub("^\\.", "", names(colors)[head(i, 1)])
      }
      
      
      for (df in list(top_tracks, left_tracks)) {
        
        
        # colors = c(Age="mako", ...)
        #________________________________________________________
        for (x in intersect(keys, colnames(df)))
          track_colors[[x]] <- list(palette = colors[[x]])
        
        
        # colors = c(Male="blue", Female="red", ...)
        #________________________________________________________
        for (i in colnames(df))
          if (!is.numeric(df[[i]]))
            if (length(x <- intersect(keys, unique(df[[i]]))) > 0)
              track_colors[[i]] <- list(palette = colors[x])
      }
          
    }
    
    remove(list = intersect(ls(), c("keys", "x", "df", "i")))
    
  } else {
    stop("In plot_heatmap(), `colors` must be a list or character vector.")
  }
  
  
  if (length(grid_colors[['palette']]) == 0)
    grid_colors[['palette']] <- attr(colors, "grid_color", exact = TRUE) %||% "rocket"
  
  if (is.null(grid_colors[['name']]))
    grid_colors[['name']] <- attr(colors, "grid_name", exact = TRUE) %||% "Grid Value"
  
  # bdiv heatmaps should have the color scheme reverse compared to taxa heatmaps
  if (isTRUE(attr(colors, "grid_rev", exact = TRUE)))
    grid_colors[['palette']] <- local({
      pal <- grid_colors[['palette']]
      if (!(length(pal) == 1 && is.character(pal))) return (pal)
      if (!grepl(pattern = "^\\-", pal)) return (paste0("-", pal))
      return (substr(pal, 2, nchar(pal)))
    })
  
  remove("colors")
  
  
  #________________________________________________________
  # Subset out any missing colors
  #________________________________________________________
  for (i in seq_along(track_colors)) {
    
    track <- names(track_colors)[[i]]  # column name in top_tracks or left_tracks
    keys  <- names(track_colors[[i]])  # metadata values being mapped to colors
    
    if (is.null(keys)) next
    
    if (track %in% colnames(top_tracks) && !is.numeric(top_tracks[[track]]))
      if (length(x <- intersect(keys, unique(top_tracks[[track]]))) > 0) {
        top_tracks[[track]] <- factor(as.character(top_tracks[[track]]), levels=x)
        top_tracks          <- top_tracks[!is.na(top_tracks[[track]]),,drop=FALSE]
        mtx                 <- mtx[,rownames(top_tracks)]
      }
    
    if (track %in% colnames(left_tracks) && !is.numeric(left_tracks[[track]]))
      if (length(x <- intersect(keys, unique(left_tracks[[track]]))) > 0) {
        left_tracks[[track]] <- factor(as.character(left_tracks[[track]]), levels=x)
        left_tracks          <- left_tracks[!is.na(left_tracks[[track]]),,drop=FALSE]
        mtx                  <- mtx[rownames(left_tracks),]
      }
  }
  
  remove(list = intersect(ls(), c("i", "track", "keys", "x")))
  
  
  #________________________________________________________
  # Heuristics for setting sizes and heights of elements
  #________________________________________________________
  if (is.null(tree_height))  tree_height  <- min(dim(mtx)) * 0.15
  if (is.null(track_height)) track_height <- tree_height / 4
  
  ratio <- ratio * ncol(mtx) / nrow(mtx)
  
  if (is.null(label_size))         label_size <- pmax(5, pmin(14, 100 / dim(mtx)))
  if (!is.null(names(label_size))) label_size <- label_size[rev(order(names(label_size)))]
  label_size_x <- tail(label_size, 1)
  label_size_y <- head(label_size, 1)
  if (isTRUE(ncol(left_tracks) > 0)) label_size_x <- .c(.rep(10, ncol(left_tracks)), .rep(label_size_x, ncol(mtx)))
  if (isTRUE(ncol(top_tracks)  > 0)) label_size_y <- .c(.rep(10, ncol(top_tracks)),  .rep(label_size_y, nrow(mtx)))
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
      if (!is.null(top_tracks))
        bounds <- bounds + (nrow(mtx) / 50) + track_height_top * ncol(top_tracks)
      
      trees_cols <- dendro(hc, bounds, side = "top")
    }
    remove("dm", "hc")
  }
  
  if (!isFALSE(clust_rows) && !is.na(clust_rows)) {
    
    dm <- dist_rows
    hc <- clust_rows
    if (!is(dm, 'dist') && !is(hc, 'hclust')) dm <- stats::dist(mtx,  method = dm)
    if (!is(hc, 'hclust'))                    hc <- stats::hclust(dm, method = hc)
    mtx <- mtx[hc[['order']],]
    
    if (isTRUE(trees_rows)) {
      
      bounds <- 0.5 - c(0, tree_height_left)
      if (!is.null(left_tracks))
        bounds <- bounds - (ncol(mtx) / 50) - track_height_top * ncol(left_tracks)
      
      trees_rows <- dendro(hc, bounds, side = "left")
    }
    remove("dm", "hc")
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
  if (!is.null(top_tracks))
    top_tracks <- tracks(
      df     = top_tracks[colnames(mtx),,drop=FALSE],
      bounds = nrow(mtx) + 0.5 + (nrow(mtx) / 50) + c(0, track_height_top * ncol(top_tracks)),
      side   = "top" )
  
  if (!is.null(left_tracks))
    left_tracks <- tracks(
      df     = left_tracks[rownames(mtx),,drop=FALSE],
      bounds = 0.5 - (ncol(mtx) / 50) - c(0, track_height_left * ncol(left_tracks)),
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
    
    if (length(left_tracks) > 0) {
      
      xlabels <- sapply(left_tracks, function (x) x[['label']])
      xbreaks <- sapply(left_tracks, function (x) x[['label_at']]) %>% 
        signif(floor(log10(nrow(mtx))) + 3)
      
      if (isTRUE(label_cols)) {
        
        xlabels <- structure(
          .Data   = c(xlabels, levels(df[['x']])),
          display = sprintf(fmt = "c(%s, levels(data[['x']]))", 
            paste(collapse = ", ", glue::single_quote(xlabels)) ))
        
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
    
    if (length(top_tracks) > 0) {
      
      ylabels <- sapply(top_tracks, function (x) x[['label']])
      ybreaks <- sapply(top_tracks, function (x) x[['label_at']]) %>% 
        signif(floor(log10(ncol(mtx))) + 3)
      
      if (isTRUE(label_rows)) {
        
        ylabels <- structure(
          .Data   = c(ylabels, levels(df[['y']])),
          display = sprintf(fmt = "c(%s, levels(data[['y']]))", 
            paste(collapse = ", ", glue::single_quote(ylabels)) ))
        
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
    if (!is.null(xlabels)) x_args %<>% c(list(labels = xlabels, breaks = xbreaks))
    if (!is.null(ylabels)) y_args %<>% c(list(labels = ylabels, breaks = ybreaks, position = "right"))
    
    gglayers %<>% ggpush(do.call(scale_x_continuous, x_args))
    gglayers %<>% ggpush(do.call(scale_y_continuous, y_args))
    
    remove("xlabels", "xbreaks", "ylabels", "ybreaks", "x_args", "y_args")
    
  }
  
  
  
  #________________________________________________________
  # Apply color gradient to the main grid
  #________________________________________________________
  pal <- plot_heatmap_palette(pal = grid_colors, vals = df[['fill']])
  
  if (length(all_tracks) > 0)
    pal[['args']][['guide']] <- guide_colorbar(
      direction      = 'horizontal', 
      title.position = 'top', 
      barheight      = unit(21.5, "points"),
      label.theme    = element_text(size = unit(8, "points")),
      order          = length(all_tracks) + 1)
  
  pal[['args']][['.indent']] <- 4
  gglayers %<>% ggpush(do.call(pal[['fun']], pal[['args']]))
  
  remove("pal", "grid_colors")
  
  
  
  if (!is.null(all_trees)) {
    
    attr(gglayers[[1]][['data']], 'trees') <- all_trees
    
    gglayers %<>% ggpush(geom_segment( 
      data    = ~ attr(., 'trees', exact = TRUE),
      mapping = aes(x=x, y=y, xend=xend, yend=yend),
      .indent = 4 ))
  }
  
  
  for (i in seq_along(all_tracks)) {
    
    track <- all_tracks[[i]]
    
    label   <- track[['label']]
    data    <- track[['data']]
    mapping <- track[['mapping']]
    outline <- track[['outline']]
    
    pal <- local({
      pal <- track_colors
      if (is.list(pal)) pal <- pal[[label]]
      if (is.null(pal)) pal <- i
      plot_heatmap_palette(pal = pal, vals = data[['fill']], title = label)
    })
    
    if (is.numeric(data[['fill']])) {
      args <- colorbar %||% list(
        direction      = 'horizontal', 
        title.position = 'top',
        barheight      = unit(21.5, "points"),
        label.theme    = element_text(size = unit(8, "points")) )
      args[['order']] <- i
      pal[['args']][['guide']] <- do.call(guide_colorbar, args)
    }
    
    attr(gglayers[[1]][['data']], track[['id']]) <- data
    data <- as.formula(sprintf("~attr(., '%s', exact = TRUE)", track[['id']]))
    
    gglayers %<>% ggpush(new_scale_fill())
    gglayers %<>% ggpush(geom_tile(data = data, mapping = mapping, .indent = 4))
    gglayers %<>% ggpush(do.call(pal[['fun']], c(pal[['args']], .indent = 4)))
    
    if (!is.numeric(attr(gglayers[[1]][['data']], track[['id']])[['fill']])) {
      args <- legend %||% list(ncol = 2)
      args[['order']] <- i
      gglayers %<>% ggpush(guides(fill = do.call(guide_legend, args)))
    }
    
    #________________________________________________________
    # Outline around the track
    #________________________________________________________
    gglayers %<>% ggpush(geom_rect(
      mapping = outline, color = "black", fill = NA, size = 0.2 ))
  }
  
  
  #________________________________________________________
  # Theming
  #________________________________________________________
  args <- theme_args %||% list()
  args <- within(args, {
    
    if (isTRUE(label_rows) || length(top_tracks) > 0)
      axis.text.y %||=% suppressWarnings(
        element_text(size=label_size_y, hjust=0) )
    
    if (isTRUE(label_cols) || length(left_tracks) > 0)
      axis.text.x %||=% suppressWarnings(
        element_text(size=label_size_x, angle=-30, vjust=1, hjust=0) )
  })
  
  
  gglayers %<>% ggpush(theme_void())
  gglayers %<>% ggpush(do.call(theme, c(args, .indent = 4)))
  
  remove("args")
  
  p <- ggbuild(gglayers)
  p$plot_env <- emptyenv()
  
  return (p)
}




#________________________________________________________
# Convert palette args to scale_fill_* function
#________________________________________________________
plot_heatmap_palette <- function (pal, vals, title = NULL) {
  
  # Arguments to pass on to e.g. scale_fill_manual
  args <- if (is.list(pal)) pal else list()
  if (is.null(args[['name']])) args[['name']] <- title
  
  # The name of the palette, or the custom colors to use
  if (is.list(pal)) pal <- pal[['palette']]
  stopifnot(length(pal) > 0)
  
  # If pal is a number, auto-select a palette
  if (is.numeric(pal)) pal <- palette_autoselect(vals, i = pal)
  stopifnot(is.character(pal))
  
  
  fn <- NULL
  
  if (is.numeric(vals)) {
    
    binned <- !is.null(args[['n.breaks']])
    
    # Create color gradient from >= 2 color codes
    if (length(pal) == 2) {
      args %<>% c('low' = pal[[1]], 'high' = pal[[2]])
      fn <- if (binned) "scale_fill_steps" else "scale_fill_gradient"
    } else {
      codes <- palette_colorcodes(pal, direction=args[['direction']])
      fn    <- if (binned) "scale_fill_stepsn" else "scale_fill_gradientn"
      if (length(codes) > 10) {
        i <- scales::rescale(seq_along(codes), to = c(1, 12))
        codes <- codes[which(!duplicated(as.integer(i)))[2:11]]
      }
      
      args[['colours']] <- codes
    }
  
  } else {
    # Create a discrete color scale
    fn <- "scale_fill_manual"
    args[['values']] <- palette_colorcodes(pal, n = length(unique(vals)))
  }
  
  args[['palette']]   <- NULL
  args[['direction']] <- NULL
  
  return (list(fn = fn, fun = eval(parse(text=fn)), args = args))
}


#________________________________________________________
# Choose a suitable palette for a set of values.
#________________________________________________________
palette_autoselect <- function (vals, i) {
  
  if (is.numeric(vals)) { # Continuous variable
    
    opts <- list(
      lo = paste0("grDevices::", c("Oranges", "Purples", "Greens", "Reds", "Blues", "Grays")),
      hi = paste0("grDevices::", c("Oranges", "Purples", "Greens", "Reds", "Blues", "Grays")),
      bi = paste0("khroma::",    c("cork", "vik", "bam")),
      un = paste0("viridis::",   c("turbo", "viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket")) )
    
    n <- length(vals)
    x <- scales::rescale(vals, to = c(-1, 1))
    
    opt <- if (sum(x > 0) < n * 2/3)        { "lo"
    } else if (sum(x > 0) > n * 2/3)        { "hi"
    } else if (sum(abs(x) > 0.5) > n * 2/3) { "bi"
    } else {                                  "un" }
    
    pal <- opts[[opt]][[(i - 1) %% length(opts[[opt]]) + 1]]
    pal <- palette_colorcodes(pal)
    if (opt == "lo") pal <- rev(pal)
    
  } else { # Categorical variable
    
    opts <- paste0("khroma::", c("muted", "light", "bright", "pale", "dark"))
    pal  <- opts[[(i - 1) %% length(opts) + 1]]
  }
  
  return (pal)
}


#________________________________________________________
# Convert user-provided palette name to vector of colors.
#________________________________________________________
palette_colorcodes <- function (pal, n = 10, direction = 1) {
  
  if (length(pal) != 1) return (pal) # Common passthrough case
  pal <- as.vector(pal)
  
  cols <- c("package", "palette", "type", "mode")
  opts <- rbind(
    data.frame(paletteer::palettes_c_names, mode = "c")[,cols], 
    data.frame(paletteer::palettes_d_names, mode = "d")[,cols], 
    data.frame(paletteer::palettes_dynamic_names, mode = "y")[,cols] )
  opts[['id']]   <- paste(sep="::", opts[['package']], opts[['palette']])
  opts[['type']] <- c(s="seq", d="div", q="qual")[substr(opts[['type']], 1, 1)]
  
  # Favor matching to RColorBrewer, khroma, or viridis
  opts <- opts[order(!opts[['package']] %in% c("RColorBrewer", "khroma", "viridis")),]
  
  # Reverse color order when palette name is prefixed by '-'
  if (identical(substr(pal, 1, 1), "-")) {
    direction <- -1
    pal <- substr(pal, 2, nchar(pal))
  }
  
  # Match either "RColorBrewer::Oranges" or "oranges"
  pal   <- sub("okabeito", "okabe_ito", tolower(pal))
  parts <- strsplit(pal, "::", fixed = TRUE)[[1]]
  row   <- if (length(parts) > 1) {
    which(
      tolower(opts[['package']]) == parts[[1]] & 
      tolower(opts[['palette']]) == parts[[2]] )
  } else {
    which(tolower(opts[['palette']]) == parts[[1]])
  }
  
  if (length(row) == 0)
    stop("There is no paletteer palette named ", pal)
  
  # Pull the hex color values out of paletteer
  pal <- as.list(opts[head(row, 1),])
  if (pal[['mode']] == "c") {
    res <- paletteer::paletteer_c(pal[['id']], n)
    
  } else if (pal[['mode']] == "y") {
    res <- paletteer::paletteer_dynamic(pal[['id']], n)
    
  } else if (pal[['mode']] == "d") {
    res <- paletteer::paletteer_d(pal[['id']])
    
    
    if (identical(pal[['id']], "khroma::okabe_ito"))
      res <- res[c(2:8, 1)] # put black last
    
    
    if (!missing(n)) {
      
      if (n > length(res)) {
        res <- rep_len(res, n)
        
      } else if (n < length(res)) {
        if (pal[['type']] == "qual") {
          res <- res[1:n]
        } else {
          res <- res[as.integer(seq(from=1, to=length(res), length.out=n))]
        }
      }
      
    }
  }
  
  if (identical(direction, -1))
    res <- rev(res)
  
  return (substr(as.vector(res), 1, 7))
}



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
    
    df <- if (is.null(prev_ht)) {
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



tracks <- function (df, bounds=c(0,1), side="top", colors = NULL) {
  
  if (is.null(df))   return (NULL)
  if (ncol(df) == 0) return (NULL)
  
  
  #________________________________________________________
  # compute the center position of the short edge for each track
  #________________________________________________________
  n        <- ncol(df)
  bounds_w <- abs(diff(bounds)) / n
  bounds   <- (seq_len(n) - 0.5) * bounds_w + min(bounds)
  sf       <- floor(log10(nrow(df))) + 3 # Sig. fig. digits
  
  
  results <- rev(lapply(rev(seq_along(df)), function (i) {
    
    res <- if (isTRUE(side == "left")) {
      list(
        data    = data.frame(fill = df[[i]], y = seq_len(nrow(df))),
        mapping = aes(
          x      = !!signif(bounds[i], sf),
          y      = y,
          fill   = fill,
          height = 1,
          width  = !!signif(bounds_w, sf) ),
        outline = aes(
          xmin = !!signif(bounds[i] - bounds_w / 2, sf),
          xmax = !!signif(bounds[i] + bounds_w / 2, sf),
          ymin = 0.5,
          ymax = !!(nrow(df) + 0.5)
        ) )
      
    } else if (isTRUE(side == "top")) {
      list(
        data    = data.frame(fill = df[[i]], x = seq_len(nrow(df))),
        mapping = aes(
          x      = x,
          y      = !!signif(bounds[i], sf),
          fill   = fill,
          height = !!signif(bounds_w, sf),
          width  = 1 ),
        outline = aes(
          xmin = 0.5,
          xmax = !!(nrow(df) + 0.5),
          ymin = !!signif(bounds[i] - bounds_w / 2, sf),
          ymax = !!signif(bounds[i] + bounds_w / 2, sf) 
        ))
      
    } else {
      stop("Annotation track `side` must be either 'left' or 'top'.")
    }
    
    res[['label']]    <- names(df)[[i]]
    res[['label_at']] <- signif(bounds[i], sf)
    res[['side']]     <- side
    res[['id']]       <- paste0(side, "_track_", i)
    
    return (res)
  }))
  
  return (results)
}
