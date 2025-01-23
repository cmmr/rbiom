

#________________________________________________________
# Log ggplot commands, complete with arguments.
#________________________________________________________

ENV <- environment()

.onLoad <- function(libname, pkgname) { # nocov start
  
  
  
  # Add imports of imports to our namespace.
  # Triggers NOTE about 20+ imports if they're in Imports.
  #
  # ggplot2 -> cli, glue, rlang, tibble
  # ggplot2 -> scales -> R6::R6Class
  # ggplot2 -> scales -> labeling::extended
  # ggplot2 -> scales -> lifecycle::deprecate_warn
  #________________________________________________________
  include <- function (pkg, ...) {
    
    require_package(pkg, reason = 'to run rbiom')
    
    for (i in c(...))
      assign(i, getFromNamespace(x = i, ns = pkg), ENV)
  }
  
  include("cli",     "cli_text", "cli_abort", "cli_warn", "qty")
  include("glue",    "glue", "single_quote", "double_quote")
  include("tibble",  "tibble", "as_tibble")
  include("rlang", c(
    "%||%", ":=", ".data", "hash", "env_names", "env_has", "is_empty",
    "is_na", "is_null", "is_bare_environment", "is_list", "is_formula", 
    "is_true", "is_false", "is_logical", "is_scalar_logical", 
    "is_character", "is_scalar_character", "is_string", 
    "is_integerish", "is_scalar_integerish", "is_double", "is_scalar_double" ))
  
  
  # rhdf5 is an undeclared Suggests (fails CRAN checks).
  #________________________________________________________
  rhdf5_funcs <- c(
    "h5createFile", "h5createGroup", "H5Dclose", "H5Dopen", "H5Fclose", "H5Fis_hdf5", 
    "H5Fopen", "h5ls", "h5readAttributes", "h5writeAttribute", "h5writeDataset" )
  
  if (nzchar(system.file(package = 'rhdf5'))) {
    
    for (i in rhdf5_funcs)
      assign(i, getFromNamespace(x = i, ns = 'rhdf5'), ENV)
    
  } else {
    
    f <- function (...)
      package_missing('rhdf5', 'to read/write HDF5 files.')
    
    for (i in rhdf5_funcs) assign(i, f, ENV)
  }
  
  
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
  lapply(FUN = cmd_wrap, pkg="ggbeeswarm",  {c('geom_beeswarm', 'geom_quasirandom')})
  lapply(FUN = cmd_wrap, pkg="ggdensity",   {c('geom_hdr', 'geom_hdr_lines')})
  lapply(FUN = cmd_wrap, pkg="ggnewscale",  {c('new_scale_fill')})
  lapply(FUN = cmd_wrap, pkg="ggrepel",     {c('geom_text_repel', 'geom_label_repel')})
  lapply(FUN = cmd_wrap, pkg="ggtext",      {c('element_markdown', 'geom_textbox')})
  lapply(FUN = cmd_wrap, pkg="graphics",    {c('pairs')})
  lapply(FUN = cmd_wrap, pkg="grid",        {c('arrow', 'unit')})
  lapply(FUN = cmd_wrap, pkg="patchwork",   {c('free', 'inset_element', 'wrap_plots')})
  lapply(FUN = cmd_wrap, pkg="scales",      {c('alpha', 'trans_new', 'label_number', 'cut_si')})
  
  lapply(FUN = basewrap, pkg="base", {c('c', 'rep')})
  
  
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
          method  = "center") ))
  
  
  #____________________________________________________________________
  # Extract attributes *WITH* exact matching by default.
  #____________________________________________________________________
  assign('attr', pos = ENV, function (x, which, exact = TRUE) {
    base::attr(x, which, exact)
  })
  
  
  #____________________________________________________________________
  # Empty the cache (mainly for during development)
  #____________________________________________________________________
  # if (!is.null(x <- get_cache_dir()))
  #   unlink(x = dir(x, full.names = TRUE))
  
} # nocov end
