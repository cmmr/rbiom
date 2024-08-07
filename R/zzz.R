

#________________________________________________________
# Log ggplot commands, complete with arguments.
#________________________________________________________

ENV <- environment(NULL)

.onLoad <- function(libname, pkgname) {
  
  
  # Add imports of imports to our namespace (avoiding NOTE about 20+ imports).
  
  # 'ggplot2' -> 'scales' -> R6::R6Class
  # 'ggplot2' -> 'scales' -> labeling::extended
  # 'ggplot2' -> 'scales' -> lifecycle::deprecate_warn
  # 'ggplot2' ->  'cli', 'glue', 'rlang', 'tibble'
  #________________________________________________________
  for (fn in c("cli_text", "cli_abort", "cli_warn", "qty"))
    assign(x = fn, value = getFromNamespace(x = fn, ns = "cli"), pos = ENV)
  
  for (fn in c("glue", "single_quote", "double_quote"))
    assign(x = fn, value = getFromNamespace(x = fn, ns = "glue"), pos = ENV)
  
  for (fn in c("tibble", "as_tibble"))
    assign(x = fn, value = getFromNamespace(x = fn, ns = "tibble"), pos = ENV)
  
  for (fn in c(
   "%||%", ":=", ".data", "hash", "env_names", "env_has", 
   "is_na", "is_null", "is_bare_environment", "is_list", "is_formula", 
   "is_true", "is_false", "is_logical", "is_scalar_logical", 
   "is_character", "is_scalar_character", "is_string", 
   "is_integerish", "is_scalar_integerish", "is_double", "is_scalar_double" ))
    assign(x = fn, value = getFromNamespace(x = fn, ns = "rlang"), pos = ENV)
  
  

  
  
  lapply(FUN = cmd_wrap, pkg="ggplot2", {c(
    'ggplot', # `aes` is intentionally omitted here
    'annotate', 'annotation_custom', 
    'coord_fixed', 'coord_flip', 'continuous_scale',
    'element_blank', 'element_rect', 'element_text', 
    'facet_grid', 'facet_wrap', 'expansion', 'labs', 'margin', 
    'geom_bar', 'geom_col', 'geom_boxplot', 'geom_violin',  
    'geom_line', 'geom_hline', 'geom_vline', 'geom_segment', 'geom_ribbon',
    'geom_crossbar', 'geom_errorbar', 'geom_linerange', 'geom_point', 
    'geom_pointrange', 'geom_raster', 'geom_rect', 'geom_tile', 
    'geom_smooth', 'geom_text', 'geom_label', 
    'guide_colorbar', 'guide_legend', 'guides', 
    'position_dodge', 'position_jitter', 'position_jitterdodge',
    'scale_color_manual', 'scale_color_continuous', 'scale_color_gradientn',
    'scale_fill_manual', 'scale_fill_gradient', 'scale_fill_gradientn', 
    'scale_fill_steps', 'scale_fill_stepsn', 
    'scale_shape_manual', 'scale_size',
    'scale_x_continuous', 'scale_x_discrete', 
    'scale_y_continuous', 'scale_y_discrete', 'scale_y_log10', 
    'stat_ellipse', 'stat_smooth', 
    'theme', 'theme_bw', 'theme_void' )})

  # lapply(FUN = cmd_wrap, pkg="ggtree", {c(
  #   'ggtree', 'geom_tiplab', 'geom_cladelab', 
  #   'gheatmap', 'hexpand', 'vexpand' )})
  
  
  lapply(FUN = cmd_wrap, pkg="fillpattern", {c('fill_pattern', 'scale_fill_pattern')})
  lapply(FUN = cmd_wrap, pkg="generics",    {c('tidy', 'glance', 'augment')})
  lapply(FUN = cmd_wrap, pkg="ggbeeswarm",  {c('geom_beeswarm', 'geom_quasirandom')})
  lapply(FUN = cmd_wrap, pkg="ggdensity",   {c('geom_hdr', 'geom_hdr_lines')})
  lapply(FUN = cmd_wrap, pkg="ggh4x",       {c('facetted_pos_scales')})
  lapply(FUN = cmd_wrap, pkg="ggnewscale",  {c('new_scale_fill')})
  lapply(FUN = cmd_wrap, pkg="ggrepel",     {c('geom_text_repel', 'geom_label_repel')})
  lapply(FUN = cmd_wrap, pkg="ggtext",      {c('element_markdown', 'geom_textbox')})
  lapply(FUN = cmd_wrap, pkg="graphics",    {c('pairs')})
  lapply(FUN = cmd_wrap, pkg="grid",        {c('arrow', 'unit')})
  lapply(FUN = cmd_wrap, pkg="patchwork",   {c('inset_element')})
  lapply(FUN = cmd_wrap, pkg="scales",      {c('alpha', 'trans_new', 'label_number', 'cut_si')})
  
  lapply(FUN = basewrap, pkg="base", {c('c', 'rep', 'summary')})
  
  
  
  #________________________________________________________
  # Attach function names as attribute
  #________________________________________________________
  for (i in ls(ENV))
    if (is.function(ENV[[i]]))
      attr(ENV[[i]], 'fn') <- i
  
  
  #________________________________________________________
  # Trigger/catch the once-per-session ggbeeswarm warning.
  #________________________________________________________
  rlang::catch_cnd(
    ggplot2::ggsave(
      filename = nullfile(), 
      device = 'png', width = 7, height = 7, units = "in", 
      plot = ggplot2::ggplot() +
        ggbeeswarm::geom_beeswarm(
          mapping = ggplot2::aes(x = rep(LETTERS[1:3],50), y = 1:150), 
          method="center") ))
  
  
  
  
  #____________________________________________________________________
  # Extract attributes *WITH* exact matching by default.
  #____________________________________________________________________
  assign('attr', pos = ENV, function (x, which, exact = TRUE) {
    base::attr(x, which, exact)
  })
  
  
  
  #____________________________________________________________________
  # Empty the cache (mainly for during development)
  #____________________________________________________________________
  if (!is.null(x <- get_cache_dir()))
    unlink(x = dir(x, full.names = TRUE))
  
}
