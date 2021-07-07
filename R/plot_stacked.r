
plot_stacked <- function (biom, x, y, facet.by = NULL, colors = NULL, facets = NULL, ...) {
  
  params <- c(as.list(environment()), list(...))
  mode   <- paste(attr(y, 'mode', exact = TRUE), "~", attr(x, 'mode', exact = TRUE))
  
  # 
  # p <- ggplot(TaxaTable, aes(x=SampleID, y=Abundance, fill=Taxa)) +
  #   geom_bar(stat="identity") + 
  #   scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
  #   labs(x=ifelse(labelBy=="None", "", labelBy), y="Relative Abundance")
  # 
  # # Set the colors
  # p <- p + scale_fill_manual(values=taxaPalette, name='Taxa')
  # 
  # # Theme modifications
  # p <- p + GGSTYLE
  # p <- p + theme(
  #   plot.margin = unit(c(1,1,1,1), "lines"), 
  #   panel.grid  = element_blank()
  # )
  # 
  # # Set the x-axis labels
  # labels <- NULL
  # if (!is.null(labelBy)) labels <- setNames(as.character(TaxaTable[[labelBy]]), TaxaTable[['SampleID']])
  # p <- p + scale_x_discrete(labels=labels)
  # 
  # # Faceting
  # nrows <- 1
  # if (!is.null(facetBy)) {
  #   p <- p + facet_wrap(facetByTicked, ncol=facetsPerRow, scales="free_x")
  #   nrows <- ceiling(length(unique(TaxaTable[[facetBy]])) / facetsPerRow)
  # }
  # 
  # # x-axis Text Angle
  # x <- req(as.numeric(labelAngle) * -1)
  # if (x == -90)         { p <- p + theme(axis.text.x=element_text(angle=x, vjust=0.3, hjust=0)) }
  # if (x > -90 && x < 0) { p <- p + theme(axis.text.x=element_text(angle=x, vjust=1,   hjust=0)) }
  # 
  
}
