#' Converts a data.frame to a download link.
#' 
#' @noRd
#' 
#' @param df   An object coercable to a \code{data.frame}.
#' @param title   Text/HTML label for the link. Default: \code{"Download Data (CSV)"}.
#' @param filename   Default filename to download as. Default: \code{"data.csv"}.
#' @param ...   Additional arguments to pass to \code{write.csv()}.
#' @return An HTML string in which the object is encoded in base64.
#' @family embed
#' 
#' @examples
#'     library(rbiom)
#'     
#'     biom <- rarefy(hmp50)
#'     
#'     p <- bdiv_ord_plot(biom, stat.by="Sex")
#'     
#'     attr(p, 'stats') %>% embed_csv(row.names=FALSE) %>% cat("\n\n")
#'

embed_csv <- function (df, label="Download Data (CSV)", filename="data.csv", ...) {
  
  as.data.frame(df) %>%
    write.csv(...) %>%
    capture.output() %>%
    paste(collapse = "\n") %>%
    base64_enc() %>%
    sprintf(
      fmt = '<a download=%s href="data:text/csv;base64,%s">%s</a>', 
      double_quote(filename), ., label) %>%
    structure(html = TRUE, class = c("html", "character"))
}


#' Wraps R code with Markdown syntax highlighting.
#' 
#' @noRd
#' 
#' @param text   Some R code in a string.
#' @return Markdown-compatible syntax highlighted code block.
#' @family embed
#' 
#' @examples
#'     library(rbiom)
#'     p <- adiv_boxplot(hmp50, stat.by = "Body Site")
#'     attr(p, 'code') %>% embed_code() %>% cat("\n\n")

embed_code <- function (text) {
  paste0("~~~~ {.R}\n", text, "\n~~~~\n\n")
}
