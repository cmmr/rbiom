
#________________________________________________________
# Use patchwork::plot_wrap(), adding `cmd` and `code`.
#________________________________________________________

plot_wrap <- function (plots, cmd, ...) {
    
  code <- sprintf("p <- %s\n\n", cmd)
  
  for (i in seq_along(plots)) {
    
    if (i == 1) { data_code <- "d1 <- p$data" }
    else        { data_code <- sprintf("d%i <- p$patches$plots[[%i]]$data", i, i - 1) }
    
    patch_code <- sub(
      x           = attr(plots[[i]], 'code', exact = TRUE), 
      pattern     = "ggplot(data)", 
      replacement = sprintf("p%i <- ggplot(d%i)", i, i),
      fixed       = TRUE )
    
    code <- paste0(code, data_code, "\n", patch_code, "\n\n")
  }
  
  plist <- I(paste0(collapse = ", ", "p", seq_along(plots)))
  code  <- paste0(code, sprintf("patchwork::wrap_plots(%s)", as.args(list(plist, ...))))
  
  
  p <- patchwork::wrap_plots(plots, ...) %>%
    add_class('rbiom_plot')
  
  attr(p, 'cmd')  <- cmd
  attr(p, 'code') <- add_class(code, 'rbiom_code')
  
  return (p)
}
