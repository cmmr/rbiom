
#-----------------------------------------------
# Enable provenance tracking on ggplot functions
#-----------------------------------------------

element_blank <- function (...) {
  res <- ggplot2::element_blank(...)
  attr(res, 'display') <- capture.output(match.call())
  return (res)
}
element_rect <- function (...) {
  res <- ggplot2::element_rect(...)
  attr(res, 'display') <- capture.output(match.call())
  return (res)
}
element_text <- function (...) {
  res <- ggplot2::element_text(...)
  attr(res, 'display') <- capture.output(match.call())
  return (res)
}

element_markdown <- function (...) {
  res <- ggtext::element_markdown(...)
  attr(res, 'display') <- paste0("ggtext::", capture.output(match.call()))
  return (res)
}

arrow <- function (...) {
  res <- grid::arrow(...)
  attr(res, 'display') <- capture.output(match.call())
  return (res)
}

expansion <- function (...) {
  res <- ggplot2::expansion(...)
  attr(res, 'display') <- capture.output(match.call())
  return (res)
}
