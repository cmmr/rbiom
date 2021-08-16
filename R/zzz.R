
# Memoise functions that do heavy-lifting

.onLoad <- function(libname, pkgname) {
  
  # Set rbiom.cache like this before calling library(rbiom)
  # options(rbiom.cache = cachem::cache_mem(max_size = 50 * 1024^2))
  # or equivalently, options(rbiom.cache = 50 * 1024^2)
  user_cm <- getOption('rbiom.cache', default = 50 * 1024^2)
  
  
  if (isTRUE(user_cm == 0)) {
    return (invisible(NULL))
    
  } else if (is.numeric(user_cm)) {
    cm <- cachem::cache_mem(max_size = user_cm)
    
  } else {
    cm <- user_cm
  }
  
  # memoise is required by DESCRIPTION, but double-checking can't hurt.
  # This lookup style is faster than rownames(installed.packages())
  stopifnot(nzchar(system.file(package = 'memoise')))
  
  alpha.div   <<- memoise::memoise(alpha.div,   cache = cm)
  apcoa       <<- memoise::memoise(apcoa,       cache = cm)
  bdply       <<- memoise::memoise(bdply,       cache = cm)
  beta.div    <<- memoise::memoise(beta.div,    cache = cm)
  distill     <<- memoise::memoise(distill,     cache = cm)
  ordinate    <<- memoise::memoise(ordinate,    cache = cm)
  plot.BIOM   <<- memoise::memoise(plot.BIOM,   cache = cm)
  rarefy      <<- memoise::memoise(rarefy,      cache = cm)
  select      <<- memoise::memoise(select,      cache = cm)
  stats.table <<- memoise::memoise(stats.table, cache = cm)
  subset.BIOM <<- memoise::memoise(subset.BIOM, cache = cm)
  subtree     <<- memoise::memoise(subtree,     cache = cm)
  taxa.rollup <<- memoise::memoise(taxa.rollup, cache = cm)
  unifrac     <<- memoise::memoise(unifrac,     cache = cm)
  
  return (invisible(NULL))
}
