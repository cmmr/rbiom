#' Visualize the number of observations per sample.
#' 
#' @inherit documentation_default
#' @inherit documentation_plot_return return
#' 
#' @family rarefaction
#' @family visualization
#'        
#' @param rline   Where to draw a horizontal line on the plot, intended to show
#'        a particular rarefaction depth. Set to `TRUE` to show an 
#'        auto-selected rarefaction depth, `FALSE` to not show a line, or
#'        an integer for a custom position.
#'        Default: `TRUE`.
#' 
#' @param counts   Display the number of samples and reads remaining after
#'        rarefying to \code{rline} reads per sample. Default: `TRUE`.
#'        
#' @param labels   Show sample names under each bar. Default: `TRUE`.
#'        
#' @param trans   Y-axis transformation. Options are \code{"log10"} or 
#'        `NULL`.  Default: \code{"log10"}.
#'        
#' @param ...   Additional parameters to pass along to ggplot2 functions. 
#'        Prefix a parameter name with \code{r.} to ensure it gets 
#'        passed to (and only to) \link[ggplot2]{geom_hline}. For instance, 
#'        \code{r.color = "black"} ensures only the horizontal rarefaction line 
#'        has its color set to \code{"black"}.
#' 
#' 
#' @export
#' @examples
#'     library(rbiom)
#'     
#'     rare_stacked(hmp50)
#'     
#'     rare_stacked(hmp50, rline = 500, r.size = 2, r.linetype = "twodash")
#'     
#'     fig <- rare_stacked(hmp50, counts = FALSE)
#'     fig$code
#'     

rare_stacked <- function (
    biom, rline = TRUE, counts = TRUE, labels = TRUE, y.trans = "log10", ...) {
  
  biom <- as_rbiom(biom)
  
  
  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  params     <- slurp_env(..., .dots = TRUE)
  cache_file <- get_cache_file('rare_stacked', params)
  if (isTRUE(attr(cache_file, 'exists', exact = TRUE)))
    return (readRDS(cache_file))
  
  params <- list2env(params)
  
  
  #________________________________________________________
  # Sanity Checks
  #________________________________________________________
  with(params, {
    stopifnot(is_scalar_logical(rline) || is_scalar_integerish(rline))
    stopifnot(!is.na(rline))
  })
  
  
  #________________________________________________________
  # Split the counts into kept/dropped.
  #________________________________________________________
  with(params, {
    
    
    # Find the number of observations per sample.
    #________________________________________________________
    .ss <- sort(sample_sums(biom))
    
    
    # Default rarefaction depth.
    #________________________________________________________
    if (isTRUE(rline)) {
      rline <- (sum(.ss) * .1) / length(.ss)
      rline <- min(.ss[.ss >= rline])
    }
    
    
    if (isFALSE(rline)) {
      
      .ggdata <- data.frame(
        check.names = FALSE,
        '.x'    = factor(names(.ss), levels = names(.ss)), 
        '.xmin' = seq_along(.ss) - 0.4,
        '.xmax' = seq_along(.ss) + 0.4,
        '.ymin' = 1,
        '.ymax' = as.vector(.ss) )
      
    } else {
      .ggdata <- data.frame(
        check.names = FALSE,
        '.x'     = factor(names(.ss), levels = names(.ss)),
        '.group' = factor(rep(c("Excluded", "Retained"), each = length(.ss))),
        '.xmin'  = seq_along(.ss) - 0.4,
        '.xmax'  = seq_along(.ss) + 0.4,
        '.ymin'  = as.vector(c(ifelse(.ss  < rline, 1, rline), ifelse(.ss < rline, 0, 1))),
        '.ymax'  = as.vector(c(ifelse(.ss == rline, 0, .ss),   ifelse(.ss < rline, 0, rline))) )
    }
    
    .ggdata %<>% subset(.ymin > 0 & .ymax > 0) %>% as_rbiom_tbl()
    
    .xcol  <- ".sample"
    .ycol  <- ".ymax"
    .xmode <- "factor"
    
    if (is.null(.dots$labs.x)) .dots$labs.x <- "Sample"
    if (is.null(.dots$labs.y)) .dots$labs.y <- "Sequencing Depth"
  })
  
  
  
  #________________________________________________________
  # Initialize the `layers` object.
  #________________________________________________________
  init_layers(params)
  
  
  #________________________________________________________
  # Add default layer parameters.
  #________________________________________________________
  set_layer(params, 'rect',  mapping = .qw(xmin, xmax, ymin, ymax), color=NA)
  set_layer(params, 'yaxis', expand = c(0,0))
  set_layer(params, 'theme', panel.grid.major.x = element_blank())
  
  # if (params$trans == "log10") {
  #   set_layer(params, 'yaxis', c(loglabels(params$.ss), trans="log10"))
  # } else {
  #   set_layer(params, 'yaxis', labels = label_number(scale_cut = cut_si("")))
  # }
  
  if (isTRUE(labels)) {
    set_layer(params, 'point', mapping     = list(x = ".x"), y = 1, alpha = 0)
    set_layer(params, 'theme', axis.text.x = element_text(angle=-30, vjust=1, hjust=0) )
    set_layer(params, 'labs',  x = NULL)
  }
  
  if (!isFALSE(params$rline)) {
    
    rline <- params$rline
    ss    <- params$.ss
    
    set_layer(params, 'rect', 'mapping|fill' = ".group")
    set_layer(params, 'labs',  fill = "Reads")
    set_layer(params, 'hline', yintercept = rline, color = "red", linetype="dashed")
    
    if (isTRUE(params$counts)) {
      
      set_layer(params, 'labs', subtitle = local({
        
        samples_before = length(ss)
        samples_after  = sum(ss >= rline)
        reads_before   = sum(ss)
        reads_after    = rline * sum(ss >= rline)
        
        samples_retained = sprintf(
          fmt = "Samples Retained: %s / %s = %s%%",
          formatC(samples_after,  format="d", big.mark=","),
          formatC(samples_before, format="d", big.mark=","),
          round(100 * samples_after / samples_before, 1) )
        reads_retained = sprintf(
          fmt = "Reads Retained: %s / %s = %s%%",
          formatC(reads_after,  format="d", big.mark=","),
          formatC(reads_before, format="d", big.mark=","),
          round(100 * reads_after / reads_before, 1) )
        
        return (paste0(samples_retained, "\n", reads_retained))
      }))
      
      set_layer(params, 'theme', plot.subtitle = element_text(size=10)) 
    }
      
    remove("rline", "ss")
  }
  
  
  
  # Convert layer definitions into a plot.
  #________________________________________________________
  p <- plot_build(params)
  
  
  attr(p, 'cmd') <- current_cmd('rare_stacked')
  set_cache_value(cache_file, p)
  
  return (p)
}






