
plot_ordination <- function (biom, x, y, layers = "box", color.by = NULL, shape.by = NULL, facet.by = NULL, colors = NULL, shapes = NULL, weighted = TRUE, k = 2, ...) {
  
  df <- ordinate(biom, x, y, weighted, k=2)
  
  gg <- ggplot(df, aes(x=Axis.1, y=Axis.2)) + geom_point()
  
  return(gg)
  
}

