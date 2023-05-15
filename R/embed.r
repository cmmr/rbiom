#' Converts a data.frame to a download link.
#' 
#' @param x   An object coercable to a \code{data.frame}.
#' @param title   Text/HTML label for the link. Default: \code{"Download Data (CSV)"}.
#' @param filename   Default filename to download as. Default: \code{"data.csv"}.
#' @param ...   Additional arguments to pass to \code{write.csv()}.
#' @return An HTML string in which the object is encoded in base64.
#' @family embed
#' @export
#' @examples
#'     library(rbiom)
#'     p <- plot(hmp50, Bray ~ Sex)
#'     attr(p, 'stats') %>% embed_csv(row.names=FALSE) %>% cat("\n\n")
#' @md
#'

embed_csv <- function (x, label="Download Data (CSV)", filename="data.csv", ...) {
  
  as.data.frame(x) %>%
    write.csv(...) %>%
    capture.output() %>%
    paste(collapse = "\n") %>%
    openssl::base64_encode() %>%
    sprintf(
      fmt = '<a download=%s href="data:text/csv;base64,%s">%s</a>', 
      double_quote(filename), ., label) %>%
    structure(html = TRUE, class = c("html", "character"))
}


#' Wraps R code with Markdown syntax highlighting.
#' 
#' @param x   Some R code.
#' @return Markdown-compatible syntax highlighted code block.
#' @family embed
#' @export
#' @examples
#'     library(rbiom)
#'     p <- adiv_boxplot(hmp50, color.by = "Body Site")
#'     attr(p, 'cmd') %>% embed_code() %>% cat("\n\n")

embed_code <- function (x) {
  paste0("~~~~ {.R}\n", x, "\n~~~~\n\n")
}
