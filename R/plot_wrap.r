
#________________________________________________________
# Combine single plots into one figure with patchwork::wrap_plots().
#
# Each plot in `plots` must have been built in "panel mode", i.e. with
# `params$.data_var = "d<i>"` and `params$.plot_var = "p<i>"` for its
# position `i`, so that its code reads `p1 <- ggplot(d1, ...)` and its
# derived tables are computed from `d1`. The combined code is then:
#
#   d1 <- <command that produces panel 1's data>   (when there is one)
#   <panel 1's code>
#   d2 <- ...
#   <panel 2's code>
#   patchwork::wrap_plots(p1, p2, ...)
#
# The caller attaches `cmd` and the rbiom object's history afterwards
# with prov_attach(). patchwork stores the last plot as `p` itself and
# the others in `p$patches$plots`.
#________________________________________________________

plot_wrap <- function (plots, ...) {

  code <- character()

  for (i in seq_along(plots)) {

    data_cmd <- attr(plots[[i]], 'data_cmd', exact = TRUE)

    code <- c(
      code,
      if (!is.null(data_cmd)) sprintf("d%i <- %s", i, data_cmd),
      attr(plots[[i]], 'code', exact = TRUE),
      "" )
  }

  plist <- I(paste0("p", seq_along(plots), collapse = ", "))
  code  <- c(code, sprintf("patchwork::wrap_plots(%s)", as.args(list(plist, ...))))


  p <- patchwork::wrap_plots(plots, ...) %>%
    add_class('rbiom_plot')

  attr(p, 'code') <- add_class(paste(code, collapse = "\n"), 'rbiom_code')

  return (p)
}



#________________________________________________________
# Mark a params environment as building panel `i` of a multi-panel
# figure (see plot_wrap()).
#________________________________________________________
panel_params <- function (params, i) {
  params$.data_var <- paste0("d", i)
  params$.plot_var <- paste0("p", i)
  return (invisible(params))
}
