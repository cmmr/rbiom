



#____________________________________________________________________
# Extract attributes *WITH* exact matching by default.
#____________________________________________________________________
attr <- function (x, which, exact = TRUE) {
  base::attr(x, which, exact)
}
