

#________________________________________________________
# Log ggplot commands, complete with arguments.
#________________________________________________________

ENV <- environment()

.onLoad <- function(libname, pkgname) { # nocov start
  
  
  #________________________________________________________
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
    for (fn in c(...))
      assign(fn, getFromNamespace(x = fn, ns = pkg), ENV)
  }
  
  include("cli",     "cli_text", "cli_abort", "cli_warn", "qty")
  include("glue",    "glue", "single_quote", "double_quote")
  include("tibble",  "tibble", "as_tibble")
  include("rlang", c(
    "%||%", ":=", ".data", "hash", "env_bind", "is_empty", "env_has", 
    "env_names", "is_na", "is_null", "is_bare_environment", "is_list", 
    "is_formula", "is_true", "is_false", "is_logical", "is_scalar_logical", 
    "is_character", "is_scalar_character", "is_string", "is_integerish", 
    "is_scalar_integerish", "is_double", "is_scalar_double" ))
  
  
  
  #________________________________________________________
  # rhdf5 is an undeclared Suggests (fails CRAN checks).
  #________________________________________________________
  if (nzchar(system.file(package = 'rhdf5'))) {
    delayedAssign('h5createFile',     getFromNamespace('h5createFile',     'rhdf5'), ENV, ENV)
    delayedAssign('h5createGroup',    getFromNamespace('h5createGroup',    'rhdf5'), ENV, ENV)
    delayedAssign('H5Dclose',         getFromNamespace('H5Dclose',         'rhdf5'), ENV, ENV)
    delayedAssign('H5Dopen',          getFromNamespace('H5Dopen',          'rhdf5'), ENV, ENV)
    delayedAssign('H5Fclose',         getFromNamespace('H5Fclose',         'rhdf5'), ENV, ENV)
    delayedAssign('H5Fis_hdf5',       getFromNamespace('H5Fis_hdf5',       'rhdf5'), ENV, ENV)
    delayedAssign('H5Fopen',          getFromNamespace('H5Fopen',          'rhdf5'), ENV, ENV)
    delayedAssign('h5ls',             getFromNamespace('h5ls',             'rhdf5'), ENV, ENV)
    delayedAssign('h5readAttributes', getFromNamespace('h5readAttributes', 'rhdf5'), ENV, ENV)
    delayedAssign('h5writeAttribute', getFromNamespace('h5writeAttribute', 'rhdf5'), ENV, ENV)
    delayedAssign('h5writeDataset',   getFromNamespace('h5writeDataset',   'rhdf5'), ENV, ENV)
  }
  
  
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
