#' Visualize diversity or abundance as a boxplot, dotplot, etc.
#' 
#' @name plot
#' @param x     A BIOM object, as returned from \link{read.biom}.
#' @param formula  Definition of the y-axis metric and metadata column to plot along the x-axis. 
#'                 For example, \bold{Shannon ~ `Body Site`}. Y-axis options are:
#'     \describe{
#'         \item{Alpha Diversity Metrics (one or more)}{
#'           \bold{OTUs}, \bold{Shannon}, \bold{Chao1}, \bold{Simpson}, and/or 
#'           \bold{InvSimpson}. Use \bold{Diversity} to get all five metrics.
#'         }
#'         \item{Beta Diversity Metrics (one only)}{
#'           \bold{manhattan}, \bold{euclidean}, \bold{bray-curtis}, \bold{jaccard}, 
#'           or \bold{unifrac}. \bold{Distance} will use \bold{unifrac} if a phylogenetic 
#'           tree is present, or \bold{bray-curtis} otherwise. Use in combination with
#'           the \code{weighted} parameter. Metadata column names can be prefixed with
#'           \bold{==} or \bold{!=} to limit distance calculations to \emph{within} or
#'           \emph{between}, respectively, those categories. See examples below.
#'         }
#'         \item{Taxa Abundances (one only)}{
#'           \bold{Kingdom}, \bold{Phylum}, \bold{Class}, \bold{Order}, \bold{Family}, 
#'           \bold{Genus}, \bold{Species}, \bold{Strain}, or \bold{OTU}. Supported ranks 
#'           will vary by biom. Run \code{taxa.ranks(biom)} to see the available options. 
#'           Specifying \bold{Abundance} will default to the most precise rank possible.
#'         }
#'        }
#' @param layers   What kind of plot to create. Options are \bold{box}, \bold{violin}, \bold{dot}, 
#'                 \bold{strip}, \bold{crossbar}, \bold{errorbar}, \bold{linerange}, and/or 
#'                 \bold{pointrange}. Single letter abbreviations are also accepted. For instance,
#'                 \code{c("box", "dot")} is equivalent to \code{c("b", "d")} and \code{"bd"}.
#' @param color.by  Metadata column to color by. If that column is a \code{factor}, the ordering of levels
#'                  will be maintained in the plot.
#' @param facet.by  A character vector of one or two metadata columns to facet by. Ordering of \code{factor}
#'                  levels will be maintained in the plot.
#' @param stats  How to display significance brackets and p-values on the plot. Options are: 
#'     \describe{
#'       \item{\code{TRUE}}{ Auto-select \bold{pairwise} or \bold{group}. (Default) }
#'       \item{\code{FALSE}}{ Do not display significance values on the plot. }
#'       \item{\bold{"p"} or \bold{"pairwise"}}{ Compare all against all. }
#'       \item{\bold{"g"} or \bold{"group"}}{ Compute a single p-value per group. }
#'       \item{\emph{data.frame}}{ Data to use for statistical annotations. I.e. from \link{stats.table}. }
#'     }
#' @param vline       How to calculate min/max of the \bold{crossbar}, \bold{errorbar}, \bold{linerange},
#'                    and \bold{pointrange} layers. Options are \bold{ci} (confidence interval), \bold{sd} 
#'                    (standard deviation), \bold{se} (standard error), and \bold{mad} (median absolute 
#'                    deviation). You may optionally append a number to \bold{ci} to specify the confidence
#'                    level, for instance \code{vline = "ci95"} (the default), will calculate the 95%
#'                    confidence interval, whereas \code{vline = "ci99"} will give the 99% confidence
#'                    interval.
#' @param xlab.angle  How to rotate the tick labels on the x-axis. \bold{'auto'} (the default), 
#'                    automatically selects a rotation value. \bold{0}, \bold{30}, and \bold{90} sets the 
#'                    angle to horizontal, angled, and vertical, respectively.
#' @param ...  Parameters passed on to ggplot2 functions. Prefixing a parameter name with \bold{b.}, 
#'             \bold{v.}, \bold{d.}, or \bold{s.} forces that parameter to be passed to, and only to,
#'             \bold{geom_boxplot}, \bold{geom_violin}, \bold{geom_dotplot}, or \bold{geom_jitter},
#'             respectively. Otherwise, parameters are passed along by matching against formal arguments.
#' @return A \code{ggplot2} plot. The computed data points and statistics will be attached as
#'         \code{attr(p, 'data')} and \code{attr(p, 'stats')}, respectively.
#' @export
#' @seealso \code{\link{stats.table}}
#' @examples
#'     library(rbiom)
#'     
#'     infile <- system.file("extdata", "hmp50.bz2", package = "rbiom")
#'     biom <- read.biom(infile)
#'     
#'     plot(biom, Shannon ~ `Body Site`)
#'     plot(biom, Shannon ~ Sex, layers="vb", color.by="Body Site")
#'     
#'     plot(biom, Diversity ~ `Body Site`, layers="p", color.by="Sex", xlab.angle=30)
#'     
#'     #plot(biom, UniFrac ~ `==Body Site`)
#'     
#'     plot(biom, bray ~ nmds)
#'     
#'
plot.BIOM <- function (x, formula, layers = "box", color.by = NULL, facet.by = NULL, stats = TRUE, vline = "ci95", xlab.angle = 'auto', ...) {
  
  biom <- x
  dots <- list(...)
  
  
  #--------------------------------------------------------------
  # Parse formula and sanity checks
  #--------------------------------------------------------------
  if (!is(biom,    'BIOM'))    stop("Please provide a BIOM object.")
  if (!is(formula, 'formula')) stop("Please provide a valid formula.")
  if (length(formula) != 3)    stop("Please provide a valid formula.")
  
  x  <- all.vars(formula[[3]])
  y  <- all.vars(formula[[2]])   %>% validate_metrics(biom, .)
  md <- c(x, color.by, facet.by) %>% validate_md_cols(biom, .)
  
  if (length(y) < 1) stop("Please provide a valid formula.")
  
  
  #--------------------------------------------------------------
  # We're running an ordination
  #--------------------------------------------------------------
  if (attr(y, 'mode') == "bdiv" && tolower(x) %in% ord_metrics()) {
    weighted <- ifelse(is.null(dots[['weighted']]), TRUE, dots[['weighted']])
    k        <- ifelse(is.null(dots[['k']]), 2, dots[['k']])
    
    ord <- ordinate(biom, x, y, weighted, md=md, k=2)
    
    p <- ggplot(ord, aes(x=Axis.1, y=Axis.2)) + geom_point()
      
    return(p)
  }
  
  
  #--------------------------------------------------------------
  # Convert biom object to a data.frame
  #--------------------------------------------------------------
  args <- list(biom = biom, metric=y, md=md)
  if ('weighted' %in% names(dots)) args[['weighted']] <- dots[['weighted']]
  
  df <- do.call(distill, args)
  y  <- attr(df, 'response')
  
  if (all(y %in% bdiv_metrics())) {
    if(!is.null(x))        x        <- sub("^[!=]=", "", x)
    if(!is.null(color.by)) color.by <- sub("^[!=]=", "", color.by)
    if(!is.null(facet.by)) facet.by <- sub("^[!=]=", "", facet.by)
  }
  
  remove("args")
  
  
  
  #--------------------------------------------------------------
  # Automatically facet by Metric or Taxa
  #--------------------------------------------------------------
  if (!is.null(attr(df, 'facet'))) {
    facet.by <- head(unique(c(attr(df, 'facet'), facet.by)), 2)
    if (is.null(dots[['scales']])) dots[['scales']] <- "free_y"
  }
  
  
  #--------------------------------------------------------------
  # Example `layers` argument: "box", c("BOX", "dot"), "V", "vbs"
  #--------------------------------------------------------------
  layerlist <- c(
    "violin", "box", "dot", "strip", 
    "crossbar", "errorbar", "linerange", "pointrange" )
  if (length(layers) == 1) layers <- strsplit(layers, "")[[1]]
  layers <- layerlist[pmatch(tolower(layers), layerlist)]
  
  
  #--------------------------------------------------------------
  # Merge sets of arguments. Preference given in following order:
  #  1. Prefixed `plot` arguments: plot(`box.fill`="red")
  #  2. `plot` arguments matching formalArgs: plot(`fill`="red")
  #  3. `withdots` arguments: withdots("box", `fill`="red")
  #--------------------------------------------------------------
  full2short <- c("facet_grid" = "facet", "facet_wrap" = "facet")
  
  short2full <- c("violin" = "geom_violin",   "crossbar"   = "geom_crossbar", 
                  "box"    = "geom_boxplot",  "errorbar"   = "geom_errorbar",
                  "dot"    = "geom_dotplot",  "linerange"  = "geom_linerange",
                  "strip"  = "geom_jitter",   "pointrange" = "geom_pointrange")
  
  short2regx <- c("violin" = "^v(|iolin)\\.", "crossbar"   = "^c(|rossbar)\\.", 
                  "box"    = "^b(|ox)\\.",    "errorbar"   = "^e(|rrorbar)\\.",
                  "dot"    = "^d(|ot)\\.",    "linerange"  = "^l(|inerange)\\.",
                  "strip"  = "^s(|trip)\\.",  "pointrange" = "^p(|ointrange)\\.",
                  "theme"  = "^t(|heme)\\.",  "facet"      = "^f(|acet)\\.")
  
  alldots <- function (fn, ...) {
    
    full  <- ifelse(fn %in% names(short2full), short2full[[fn]], fn)
    short <- ifelse(fn %in% names(full2short), full2short[[fn]], fn)
    regx  <- short2regx[[short]]
    func  <- do.call(`::`, list("ggplot2", full))
    
    specific_dots        <- dots[grep(regx, names(dots))]
    names(specific_dots) <- sub(regx, "", names(specific_dots), perl = TRUE)
    
    generic_dots <- dots[intersect(names(dots), formalArgs(func))]
    generic_dots <- generic_dots[setdiff(names(generic_dots), names(specific_dots))]
    args         <- c(generic_dots, specific_dots)
    
    defaults <- list(...)
    args     <- c(args, defaults[setdiff(names(defaults), names(args))])
    
    do.call(func, args)
  }
  
  
  #--------------------------------------------------------------
  # Help multiple plots overlay well together
  #--------------------------------------------------------------
  dodge  <- ggplot2::position_dodge(width = 0.8)
  jdodge <- ggplot2::position_jitterdodge(dodge.width = 0.8, jitter.width = 0.05)
  
  # Violin arguments
  v_args <- list(fn = "violin", position = dodge, color = "black")
  
  # Boxplot arguments
  b_args <- list(fn = "box", position = dodge, color = "black", width = 0.7)
  b_args[['width']] <- ifelse("violin" %in% layers, 0.1, 0.7)
  if ("violin" %in% layers)               b_args[['fill']]          <- "white"
  if (any(c("dot", "strip") %in% layers)) b_args[['outlier.shape']] <- NA
  
  # Dotplot arguments
  d_args <- list(fn = "dot", position = dodge, binaxis = "y", stackdir = "center")
  d_args[['binwidth']] <- diff(range(df[[y]])) / 30
  if (any(c("violin", "box") %in% layers))
    d_args <- c(d_args, color = "black", fill = "black", dotsize = 0.7)
  
  # Stripchart arguments
  s_args <- list(fn = "strip",  position = jdodge)
  if (any(c("violin", "box") %in% layers))
    s_args <- c(s_args, color = "black")
  
  # Vertical line arguments
  if (any(c("crossbar", "errorbar", "linerange", "pointrange") %in% layers)) {
    
    if (substr(vline, 1, 2) == "ci") {
      cl <- ifelse(vline == "ci", "95", substr(vline, 3, nchar(vline)))
      cl <- as.numeric(cl) / 100
      vlineFn <- function (vals) {
        tt <- try(t.test(vals, conf.level = cl), silent = TRUE)
        if (!is(tt, "htest")) tt <- list(estimate = NA, conf.int = c(NA, NA))
        data.frame(y = unname(tt$estimate), ymin = tt$conf.int[1], ymax = tt$conf.int[2])
      }
    } else if (vline == "mad") {
      vlineFn <- function (vals) {
        med <- median(vals); dev <- mad(vals, med)
        data.frame(y = med, ymin = med - dev, ymax = med + dev)
      }
    } else if (vline == "sd") {
      vlineFn <- function (vals) {
        avg <- mean(vals); dev <- sd(vals)
        data.frame(y = avg, ymin = avg - dev, ymax = avg + dev)
      }
    } else if (vline == "se") {
      vlineFn <- function (vals) {
        avg <- mean(vals); dev <- sqrt(var(vals)/length(vals))
        data.frame(y = avg, ymin = avg - dev, ymax = avg + dev)
      }
    } else {
      stop("vline must be one of 'ci95', 'mad', 'sd', or 'se'.")
    }
    
    plyby <- NULL
    for (i in c(x, color.by, facet.by))
      plyby <- c(plyr::as.quoted(as.name(i)), plyby)
    args <- list(
      data     = plyr::ddply(df, plyby, function (v) { vlineFn(v[[y]]) }),
      mapping  = aes(y=y, ymin=ymin, ymax=ymax),
      position = dodge
    )
    
    if (any(c("violin", "box") %in% layers))
      args <- c(args, color = "black")
    
    c_args <- c(list(fn = "crossbar"),   args, width = 0.5, fill = NA)
    e_args <- c(list(fn = "errorbar"),   args, width = 0.5)
    l_args <- c(list(fn = "linerange"),  args)
    p_args <- c(list(fn = "pointrange"), args)
    remove("args")
  }
  
  
  # Facet arguments
  f_args <- list(fn = "facet_wrap", facets = facet.by)
  if (length(facet.by) > 1) f_args[['fn']] <- "facet_grid"
  
  
  
  #--------------------------------------------------------------
  # Set up the global aesthetics
  #--------------------------------------------------------------
  aes.args <- list(x = as.name(x), y = as.name(y))
  if (!is.null(color.by)) {
    aes.args <- c(
      aes.args, 
      group = c(as.name(x), as.name(color.by)),
      color = as.name(color.by),
      fill  = as.name(color.by)
      )
  }
  
  
  #--------------------------------------------------------------
  # Create the plot and add each layer with its args
  #--------------------------------------------------------------
  p <- ggplot(df, do.call(aes_string, aes.args, quote = TRUE))

  for (layer in intersect(layerlist, layers)) {
    args <- get(paste0(substr(layer, 1, 1), "_args"))
  
    if (is.null(facet.by) || layer != "dot") {
      p <- p + do.call("alldots", args)
      
    } else { # faceted dotplots need extra help with binwidth
      for (i in unique(df[[facet.by]])) {
        args[['data']]     <- df[df[[facet.by]] == i,,drop=F]
        args[['binwidth']] <- diff(range(args[['data']][[y]])) / 30
        p <- p + do.call("alldots", args)
      }
    }
    
  }
  
  if (!is.null(facet.by)) p <- p + do.call("alldots", f_args)
  
  
  #--------------------------------------------------------------
  # Arguments to pass to ggplot2::theme()
  #--------------------------------------------------------------
  t_args <- list(fn = "theme")
  
  # x-axis Text Angle
  if (xlab.angle == 'auto') {
    if (!is.null(x))
      if (sum(nchar(unique(df[[x]]))) > 50)
        xlab.angle <- 30
  }
  
  if (xlab.angle == 90) {
    t_args[['axis.text.x']] <- ggplot2::element_text(angle=-90, vjust=0.3, hjust=0)
    
  } else if (xlab.angle == 30) {
    t_args[['axis.text.x']] <- ggplot2::element_text(angle=-30, vjust=1, hjust=0)
    
    # Ensure long x-axis labels don't get truncated at the figure's edge
    rpad <- strwidth(tail(levels(df[[x]]), 1), units="inches")
    rpad <- rpad * 0.8660254 # sin((90-30) * pi / 180) / sin(90 * pi / 180)
    if (!is.null(color.by))
      rpad <- rpad - max(c(
        strwidth(color.by,               units="inches", cex=1.2) + .20,
        strwidth(levels(df[[color.by]]), units="inches")          + .52 ))
    t_args[['plot.margin']] <- ggplot2::unit(x=c(.1, max(.1, rpad), .1, .1), units="inches")
  }
  
  
  p <- p + do.call("alldots", t_args)
  
  # p <- p + scale_color_discrete(aes_string(color=as.name(color.by)))
  
  
  return (p)
}



