

#________________________________________________________
# Log ggplot command history, complete with arguments.
#________________________________________________________

ENV <- environment(NULL)

.onLoad <- function(libname, pkgname) {
  
  lapply(FUN = cmd_wrap, pkg="ggplot2",    {c(
    'ggplot', 'aes_string', 
    'coord_fixed', 'coord_flip', 'continuous_scale',
    'element_blank', 'element_rect', 'element_text', 
    'facet_grid', 'facet_wrap', 'expansion', 'labs', 'margin', 
    'geom_bar', 'geom_col', 'geom_boxplot', 'geom_violin',  
    'geom_hline', 'geom_vline', 'geom_segment',
    'geom_crossbar', 'geom_errorbar', 'geom_linerange', 'geom_point', 
    'geom_pointrange', 'geom_raster', 'geom_rect', 'geom_tile', 
    'geom_smooth', 'geom_text', 'geom_label', 
    'guide_colorbar', 'guide_legend', 'guides', 
    'position_dodge', 'position_jitter', 'position_jitterdodge',
    'scale_color_manual', 'scale_color_continuous',
    'scale_fill_manual', 'scale_fill_gradient', 'scale_fill_gradientn', 
    'scale_fill_steps', 'scale_fill_stepsn', 
    'scale_shape_manual', 'scale_size',
    'scale_x_continuous', 'scale_x_discrete', 
    'scale_y_continuous', 'scale_y_discrete', 'scale_y_log10', 
    'stat_ellipse', 'stat_smooth', 
    'theme', 'theme_bw', 'theme_void' )})
  
  lapply(FUN = cmd_wrap, pkg="ggpattern",  {c(
    'geom_bar_pattern', 'geom_boxplot_pattern', 
    'geom_col_pattern', 'geom_crossbar_pattern', 
    'scale_pattern_color_manual', 'scale_pattern_fill_manual', 
    'scale_pattern_type_manual', 'geom_violin_pattern' )})
  
  lapply(FUN = cmd_wrap, pkg="ggtree",  {c(
    'ggtree', 'geom_tiplab', 'geom_cladelab', 
    'gheatmap', 'hexpand', 'vexpand' )})
  
  lapply(FUN = cmd_wrap, pkg="ggbeeswarm", {c('geom_beeswarm', 'geom_quasirandom')})
  lapply(FUN = cmd_wrap, pkg="ggdensity",  {c('geom_hdr', 'geom_hdr_lines')})
  # lapply(FUN = cmd_wrap, pkg="ggh4x",      {c('facetted_pos_scales')})
  lapply(FUN = cmd_wrap, pkg="ggnewscale", {c('new_scale_fill')})
  lapply(FUN = cmd_wrap, pkg="ggrepel",    {c('geom_label_repel')})
  lapply(FUN = cmd_wrap, pkg="ggtext",     {c('element_markdown')})
  lapply(FUN = cmd_wrap, pkg="grid",       {c('arrow', 'unit')})
  lapply(FUN = cmd_wrap, pkg="scales",     {c('alpha')})
  lapply(FUN = cmd_wrap, pkg="emmeans",    {c('emmeans', 'emtrends')})
  
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
