
# To-do: add metacoder for overlaying on a phylogenetic tree.

#' Visualize BIOM data with boxplots.
#' 
#' @name taxa_boxplot
#' 
#' @inherit adiv_boxplot params sections return
#' 
#' @family plotting
#' 
#' 
#' @param x   A categorical metadata column name to use for the x-axis. The 
#'        default, \code{NULL} puts taxa along the x-axis.
#'        
#' @param rank   What rank(s) of taxa to display. E.g. "Phylum", "Genus", etc. 
#'        Run \code{taxa_ranks()} to see all options for a given BIOM object. 
#'        The default, \code{NULL}, selects the lowest level.
#'        
#' @param taxa   Which taxa to display. An integer value will show the top n
#'        most abundant taxa. A value 0 <= n < 1 will show any taxa with that 
#'        mean abundance or greater (e.g. 0.1). A character vector of
#'        taxon names will show only those taxa. Default: \code{5}.
#'
#' @param flip   Transpose the axes, so that taxa are present as rows instead
#'        of columns. Default: \code{TRUE}
#'
#' @param stripe   Shade every other x position. Default: \emph{same as flip}
#'
#' @param p.top   Only display taxa with the most significant differences in 
#'        abundance. If \code{p.top} is >= 1, then the \code{p.top} most 
#'        significant taxa are displayed. If \code{p.top} is less than one, all 
#'        taxa with an adjusted p-value <= \code{p.top} are displayed. 
#'        Recommended to be used in combination with the \code{taxa} parameter 
#'        to set a lower bound on the mean abundance of considered taxa. 
#'        Default: \code{Inf}
#'        
#' @param y.trans   The transformation to apply to the y-axis. Visualizing 
#'        differences of both high- and low-abundance taxa is best done with
#'        a non-linear axis. Options are: 
#'        \itemize{
#'          \item{\code{"sqrt"} - }{ square-root transformation }
#'          \item{\code{"log1p"} - }{ log(y + 1) transformation }
#'          \item{\emph{anything else} - }{ no transformation }
#'        }
#'        These methods allow visualization of both high- and low-abundance
#'        taxa simultaneously, without complaint about 'zero' count
#'        observations. Default: \code{"sqrt"}
#'        
#' @param ...   Parameters are matched to formal arguments of ggplot2
#'        functions. Prefixing parameter names with a layer name ensures that
#'        a particular parameter is passed to, and only to, that layer. For
#'        instance, \code{dot.size = 2} or \code{d.size = 2} ensures only the 
#'        dotplot layer has its size set to \code{2}.
#' 
#' 
#' @export
#' @seealso \code{\link{stats_table}}
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50) 
#'     taxa_boxplot(biom, rank = c("Phylum", "Genus"))
#'     taxa_boxplot(biom, rank = "Genus", taxa = 3, layers = "ps", color.by = list("Body Site" = c('Saliva' = "blue", 'Stool' = "red")), flip = FALSE)
#'     
#'
taxa_boxplot <- function (
    biom, x = NULL, rank = NULL, taxa = 5, layers = "lsb",
    color.by = NULL, pattern.by = NULL, shape.by = NULL, facet.by = NULL, limit.by = NULL, 
    flip = TRUE, stripe = flip, p.top = Inf, p.adj = "fdr", p.label = TRUE, 
    ci = 95, xlab.angle = 'auto', y.trans = "sqrt", ...) {
  
  
  #________________________________________________________
  # Record the function call in a human-readable format.
  #________________________________________________________
  params <- c(as.list(environment()), list(...))
  params[['...']] <- NULL
  history <- attr(biom, 'history')
  history %<>% c(sprintf("taxa_boxplot(%s)", as.args(params, fun = taxa_boxplot)))
  remove(list = setdiff(ls(), c("params", "history")))
  
  
  #________________________________________________________
  # Sanity checks. x and *.by are checked by boxplot_build.
  #________________________________________________________
  params %<>% within({
    if (!is(biom, 'BIOM')) stop("Please provide a BIOM object.")
    rank %<>% validate_arg(biom, 'rank', n = c(1,Inf), default = tail(c('OTU', taxa_ranks(biom)), 1))
  })
  
  
  #________________________________________________________
  # initLayer ignores formalArgs, so copy into equivalent.
  #________________________________________________________
  if (isTRUE(tolower(params[['y.trans']]) %in% c("sqrt", "log1p")))
    params[['yaxis.trans']] <- tolower(params[['y.trans']])
  
  
  
  #________________________________________________________
  # Use the generalized boxplot function to make the plot
  #________________________________________________________
  p <- boxplot_build(params, taxa_boxplot, taxa_boxplot_data, taxa_boxplot_layers)
  
  attr(p, 'history') <- history
  
  
  return (p)
  
}


#______________________________________________________________
# Convert biom object to a data.frame
#______________________________________________________________
taxa_boxplot_data <- function (params) {
  
  biom  <- params[['biom']]
  ranks <- params[['rank']]
  taxa  <- params[['taxa']]
  
  # Convert abundance spec to taxa names
  if (is.numeric(taxa))
    taxa <- as.vector(sapply(ranks, function (rank) {
      means <- taxa_means(as_percent(biom), rank)
      if (taxa < 1) { return (names(which(means >= taxa)))
      } else        { return (head(names(means), taxa)) }
    }))
  
  if (is_rarefied(biom))
    biom <- as_percent(biom)
  
  ggdata <- plyr::ldply(ranks, function (rank) {
    
    df <- taxa_table(
      biom = biom,
      rank = rank,
      md   = TRUE,
      safe = TRUE )
    
    df <- df[df[['.taxa']] %in% taxa,,drop=FALSE]
    if (nrow(df) == 0) return (NULL)
    
    df[['.rank']] <- rank
    df[['.y']]    <- df[['.value']]
    return (df)
  })
  
  
  ggdata[['.taxa']] %<>% factor(levels = taxa)
  ggdata[['.rank']] %<>% factor(levels = ranks)
  params[['.group.by']] %<>% c(".taxa", .)
  
  
  # Default to taxa along x-axis instead of all together.
  #________________________________________________________
  if (identical(params[['x']], ".all")) {
    
    params[['x']] <- ".taxa"
    params[['.group.by']] %<>% setdiff(".all")
    
    if (length(ranks) > 1) params[['facet.by']] %<>% c(".rank", .)
    
  } else {
    params[['facet.by']] %<>% c(".taxa", .)
  }
  
  
  attr(ggdata, 'params') <- params
  attr(ggdata, 'xcol')   <- params[['x']]
  attr(ggdata, 'ycol')   <- ".y"
  
  return (ggdata)
}


#______________________________________________________________
# Make taxa-specific layer tweaks
#______________________________________________________________
taxa_boxplot_layers <- function (layers) {
  
  setLayer("labs", y = "Relative Abundance")
  setLayer("yaxis", labels = scales::percent)
  
  return (layers)
  
}
