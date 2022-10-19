

#________________________________________________________
# Log ggplot command history, complete with arguments.
#________________________________________________________

ENV <- environment(NULL)

.onLoad <- function(libname, pkgname) {
  
  lapply(FUN = ggwrap, pkg="ggplot2",    {c(
    'ggplot', 'aes_string', 
    'coord_fixed', 'coord_flip', 'continuous_scale',
    'element_blank', 'element_rect', 'element_text', 
    'expansion', 'facet_grid', 'facet_wrap', 
    'geom_bar', 'geom_boxplot', 'geom_col',  
    'geom_hline', 'geom_vline', 'geom_segment',
    'geom_crossbar', 'geom_errorbar', 'geom_linerange', 'geom_point', 
    'geom_pointrange', 'geom_raster', 'geom_rect', 
    'geom_smooth', 'geom_text', 'geom_tile', 'geom_violin', 
    'guide_colorbar', 'guide_legend', 'guides', 
    'labs', 'margin', 
    'position_dodge', 'position_jitter', 'position_jitterdodge',
    'scale_color_manual', 
    'scale_fill_manual', 'scale_fill_gradient', 'scale_fill_gradientn', 
    'scale_fill_steps', 'scale_fill_stepsn', 
    'scale_shape_manual', 
    'scale_y_log10', 
    'scale_x_continuous', 'scale_x_discrete', 
    'scale_y_continuous', 'scale_y_discrete', 
    'stat_ellipse', 'stat_smooth', 
    'theme', 'theme_bw', 'theme_void' )})
  
  lapply(FUN = ggwrap, pkg="ggpattern",  {c(
    'geom_bar_pattern', 'geom_boxplot_pattern', 
    'geom_col_pattern', 'geom_crossbar_pattern', 
    'scale_pattern_color_manual', 'scale_pattern_fill_manual', 
    'scale_pattern_type_manual', 'geom_violin_pattern' )})
  
  lapply(FUN = ggwrap, pkg="ggtree",  {c(
    'ggtree', 'geom_tiplab', 'geom_cladelab', 
    'hexpand', 'vexpand' )})
  
  lapply(FUN = ggwrap, pkg="grid",       {c('arrow', 'unit')})
  lapply(FUN = ggwrap, pkg="ggtext",     {c('element_markdown')})
  lapply(FUN = ggwrap, pkg="ggdensity",  {c('geom_hdr', 'geom_hdr_lines')})
  lapply(FUN = ggwrap, pkg="ggbeeswarm", {c('geom_beeswarm', 'geom_quasirandom')})
  lapply(FUN = ggwrap, pkg="ggrepel",    {c('geom_label_repel')})
  lapply(FUN = ggwrap, pkg="ggnewscale", {c('new_scale_fill')})
  lapply(FUN = ggwrap, pkg="scales",     {c('alpha' )})
  
  lapply(FUN = basewrap, pkg="base", {c('c', 'rep')})
  
  
  # #________________________________________________________
  # # aes() uses non-standard evaluation. Quote as is.
  # #________________________________________________________
  # assign('aes', pos = ENV, function (..., .indent = 0) {
  #   res <- ggplot2::aes(...)
  #   cmd <- trimws(capture.output(match.call()))
  #   attr(res, 'display') <- paste(collapse = " ", cmd)
  #   return (res)
  # })
  
}
