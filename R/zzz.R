
# Memoise functions that do heavy-lifting

.onLoad <- function(libname, pkgname) {
  
  # memoise covers all the necessary dependencies
  # This lookup style is faster than rownames(installed.packages())
  if (!nzchar(system.file(package = 'memoise')))
    return (invisible(NULL))
  
  cm <- getOption(
    x       = 'rbiom.cache', 
    default = try(cachem::cache_mem(max_size = 50 * 1024^2), silent=TRUE) )
  
  if (is(cm, "cachem")) {
    alpha.div   <<- memoise::memoise(alpha.div,   cache = cm)
    apcoa       <<- memoise::memoise(apcoa,       cache = cm)
    beta.div    <<- memoise::memoise(beta.div,    cache = cm)
    rarefy      <<- memoise::memoise(rarefy,      cache = cm)
    select      <<- memoise::memoise(select,      cache = cm)
    subset.BIOM <<- memoise::memoise(subset.BIOM, cache = cm)
    subtree     <<- memoise::memoise(subtree,     cache = cm)
    taxa.rollup <<- memoise::memoise(taxa.rollup, cache = cm)
    unifrac     <<- memoise::memoise(unifrac,     cache = cm)
  }
}