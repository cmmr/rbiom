
# Memoise functions that do heavy-lifting

.onLoad <- function(libname, pkgname) {
  
  # Set rbiom.cache like this before calling library(rbiom)
  # options(rbiom.cache = cachem::cache_mem(max_size = 50 * 1024^2))
  
  cm <- getOption('rbiom.cache')
  
  if (!is(cm, "cachem"))
    return (invisible(NULL))
  
  # memoise covers all the necessary dependencies
  # This lookup style is faster than rownames(installed.packages())
  if (!nzchar(system.file(package = 'memoise'))) {
    cat("The 'memoise' R package is needed to use caching.", file = stderr())
    return (invisible(NULL))
  }

  alpha.div   <<- memoise::memoise(alpha.div,   cache = cm)
  apcoa       <<- memoise::memoise(apcoa,       cache = cm)
  beta.div    <<- memoise::memoise(beta.div,    cache = cm)
  rarefy      <<- memoise::memoise(rarefy,      cache = cm)
  select      <<- memoise::memoise(select,      cache = cm)
  subset.BIOM <<- memoise::memoise(subset.BIOM, cache = cm)
  subtree     <<- memoise::memoise(subtree,     cache = cm)
  taxa.rollup <<- memoise::memoise(taxa.rollup, cache = cm)
  unifrac     <<- memoise::memoise(unifrac,     cache = cm)
  
  return (invisible(NULL))
}