
# Colorblind-friendly color sets and gradients.
# Hardcoded here to avoid additional dependencies on khroma and/or other R packages.

PALETTES <- list(

  # Discrete color palettes from various sources.
  # 'okabe' is Okabe & Ito's palette (Wong 2011, Nature Methods 8:441),
  # ordered so that the first colors contrast most with each other and with
  # a white background; the yellow, which is hard to see on white, is last.
  'bright'     = c("4477AA", "EE6677", "228833", "CCBB44", "66CCEE", "AA3377", "BBBBBB"),
  'vibrant'    = c("EE7733", "0077BB", "33BBEE", "EE3377", "CC3311", "009988", "BBBBBB"),
  'okabe'      = c("0072B2", "D55E00", "009E73", "E69F00", "CC79A7", "56B4E9", "F0E442", "999999"),
  'r4'         = c("DF536B", "61D04F", "2297E6", "28E2E5", "CD0BBC", "F5C710", "9E9E9E", "000000"),
  'light'      = c("77AADD", "EE8866", "EEDD88", "FFAABB", "99DDFF", "44BB99", "BBCC33", "AAAA00"),
  'muted'      = c("CC6677", "332288", "DDCC77", "117733", "88CCEE", "882255", "44AA99", "999933", "AA4499"),
  'fishy'      = c("6388B4", "FFAE34", "EF6F6A", "8CC2CA", "C3BC3F", "55AD89", "BB7693", "BAA094", "767676"),
  'tableau'    = c("4E79A7", "F28E2B", "E15759", "76B7B2", "59A14F", "EDC948", "B07AA1", "FF9DA7", "9C755F", "BAB0AC"),
  'classic'    = c("1F77B4", "FF7F0E", "2CA02C", "D62728", "9467BD", "8C564B", "E377C2", "7F7F7F", "BCBD22", "17BECF"),
  'carto'      = c("88CCEE", "CC6677", "DDCC77", "117733", "332288", "AA4499", "44AA99", "999933", "882255", "661100", "6699CC", "888888"),
  'tol'        = c("332288", "6699CC", "88CCEE", "44AA99", "117733", "999933", "DDCC77", "661100", "CC6677", "AA4466", "882255", "AA4499"),
  'tableau20'  = c("1F77B4", "AEC7E8", "FF7F0E", "FFBB78", "2CA02C", "98DF8A", "D62728", "FF9896", "9467BD", "C5B0D5", "8C564B", "C49C94", "E377C2", "F7B6D2", "7F7F7F", "C7C7C7", "BCBD22", "DBDB8D", "17BECF", "9EDAE5"),
  'kelly'      = c("F3C300", "875692", "F38400", "A1CAF1", "BE0032", "C2B280", "848482", "008856", "E68FAC", "0067A5", "F99379", "604E97", "F6A600", "B3446C", "DCD300", "882D17", "8DB600", "654522", "E25822", "2B3D26", "222222"),
  'alphabet'   = c("AA0DFE", "3283FE", "85660D", "782AB6", "565656", "1C8356", "16FF32", "F7E1A0", "E2E2E2", "1CBE4F", "C4451C", "DEA0FD", "FE00FA", "325A9B", "FEAF16", "F8A19F", "90AD1C", "F6222E", "1CFFCE", "2ED9FF", "B10DA1", "C075A6", "FC1CBF", "B00068", "FBE426", "FA0087"),
  'polychrome' = c("5A5156", "E4E1E3", "F6222E", "FE00FA", "16FF32", "3283FE", "FEAF16", "B00068", "1CFFCE", "90AD1C", "2ED9FF", "DEA0FD", "AA0DFE", "F8A19F", "325A9B", "C4451C", "1C8356", "85660D", "B10DA1", "FBE426", "1CBE4F", "FA0087", "FC1CBF", "F7E1A0", "C075A6", "782AB6", "AAF400", "BDCDFF", "822E1C", "B5EFB5", "7ED7D1", "1C7F93", "D85FF7", "683B79", "66B0FF", "3B00FB"),

  # Sequential and diverging palettes from khroma R package. khroma::color("acton")(9)
  'acton'      = c("2D204C", "4F3D69", "765985", "9E6592", "C36C99", "D48CB1", "D3A5C4", "D9C5D8", "E5E5F0"),
  'bamako'     = c("003F4C", "144C40", "2A5933", "426A25", "5F7D13", "868D02", "B8A525", "E3C960", "FFE599"),
  'batlow'     = c("001959", "0E395F", "215A5F", "4C724C", "7F8033", "BF9038", "F39E71", "FDB3B4", "F9CCF9"),
  'bilbao'     = c("FFFFFF", "D2D2D2", "C0BAA5", "B3A779", "AB8B66", "A16C59", "974C4B", "722626", "4C0000"),
  'buda'       = c("B200B2", "B22D9B", "BA4C90", "C56986", "CC847E", "D39F77", "DABD70", "E3DA68", "FFFF66"),
  'davos'      = c("2C194C", "264477", "3469AF", "5984BF", "7C9AB4", "A1B3AA", "C7CD9E", "E6E6BF", "FFFFFF"),
  'devon'      = c("2B194C", "27386B", "26588D", "3C6FB9", "788BD9", "ACA6ED", "C7C3F3", "E3E0F8", "FFFFFF"),
  'hawaii'     = c("8C0172", "922C54", "964C3E", "9A6D27", "9B951B", "89BB48", "6BD38C", "66E7D2", "B2F2FD"),
  'imola'      = c("1933B2", "2547A7", "305D9D", "3E708D", "52857F", "70A377", "91C36D", "BDE566", "FFFF66"),
  'lajolla'    = c("FFFFCC", "F9E98D", "F2C559", "E99D52", "E0744F", "C24C4A", "873E39", "4B2B1A", "191900"),
  'lapaz'      = c("190C65", "213682", "2B5A99", "3E7CA5", "6598A4", "96A493", "C5AC8A", "F8D0B7", "FFF2F2"),
  'nuuk'       = c("04598C", "2A6380", "527684", "7F9092", "A4A699", "B7B790", "C4C482", "DDDD86", "FEFEB2"),
  'oslo'       = c("000000", "0C1D2F", "153C5F", "265C95", "4C7CC5", "7797CA", "9FAEC8", "CBCCD0", "FFFFFF"),
  'tokyo'      = c("190D33", "4D1F50", "7A466F", "8A6B7E", "908B86", "94AC90", "9FCE99", "D2F7B9", "FFFFD8"),
  'turku'      = c("000000", "262625", "4B473E", "706A52", "9A8C66", "C99F7E", "EBA597", "FEC4C3", "FFE5E5"),
  'bam'        = c("65024B", "A4428B", "D07EBB", "ECC5E1", "F6F1F0", "D7E7C0", "8CB464", "4B802E", "0D4C00"),
  'berlin'     = c("9EB0FF", "519FD2", "276785", "13303E", "180C09", "3F1200", "7A321B", "BB6C60", "FFACAC"),
  'broc'       = c("2C194C", "284B7D", "5D83A9", "A8BED1", "EDF1EC", "D8D8AB", "9E9E65", "5E5D2F", "262600"),
  'cork'       = c("2C194C", "284A7D", "5A82A8", "A5BAD0", "E5EDEC", "ADCFB1", "6BA870", "3F772D", "434C01"),
  'lisbon'     = c("E5E5FF", "90A6CC", "406A97", "15334F", "161919", "433F27", "7F784C", "C0B986", "FFFFD8"),
  'roma'       = c("7E1900", "A0621C", "C0A439", "E3E086", "D1ECC8", "79D4D9", "479AC5", "3165AD", "1A3399"),
  'tofino'     = c("DED8FF", "879DD8", "3E5E9A", "1D2C4A", "0C1413", "1B3B1F", "37723B", "7BB268", "DAE59A"),
  'vanimo'     = c("FFCDFD", "CD78BD", "923E80", "401B37", "1A1513", "293516", "517026", "7EAC45", "BEFDA5"),
  'vik'        = c("001260", "01437F", "2D7CA5", "91BCD1", "EAECE9", "D3BF95", "AF893D", "854B00", "601200"),

  # Basic sequential palettes from base R. RColorBrewer::brewer.pal(9, "Reds")
  'reds'       = c("FFF5F0", "FEE0D2", "FCBBA1", "FC9272", "FB6A4A", "EF3B2C", "CB181D", "A50F15", "67000D"),
  'oranges'    = c("FFF5EB", "FEE6CE", "FDD0A2", "FDAE6B", "FD8D3C", "F16913", "D94801", "A63603", "7F2704"),
  'greens'     = c("F7FCF5", "E5F5E0", "C7E9C0", "A1D99B", "74C476", "41AB5D", "238B45", "006D2C", "00441B"),
  'purples'    = c("FCFBFD", "EFEDF5", "DADAEB", "BCBDDC", "9E9AC8", "807DBA", "6A51A3", "54278F", "3F007D"),
  'grays'      = c("FFFFFF", "F0F0F0", "D9D9D9", "BDBDBD", "969696", "737373", "525252", "252525", "000000")
)



#________________________________________________________
# The kind of each palette in PALETTES.
#________________________________________________________
PALETTE_KINDS <- c(
  setNames(rep("discrete", 15), c(
    'bright', 'vibrant', 'okabe', 'r4', 'light', 'muted', 'fishy', 'tableau',
    'classic', 'carto', 'tol', 'tableau20', 'kelly', 'alphabet', 'polychrome' )),
  setNames(rep("sequential", 20), c(
    'acton', 'bamako', 'batlow', 'bilbao', 'buda', 'davos', 'devon', 'hawaii',
    'imola', 'lajolla', 'lapaz', 'nuuk', 'oslo', 'tokyo', 'turku',
    'reds', 'oranges', 'greens', 'purples', 'grays' )),
  setNames(rep("diverging", 9), c(
    'bam', 'berlin', 'broc', 'cork', 'lisbon', 'roma', 'tofino', 'vanimo', 'vik' ))
)



#' List the color palettes.
#'
#' The palettes that rbiom's `colors`, `grid`, and track `colors`
#' arguments accept by name. Prefix a name with `"-"` to reverse its
#' colors, e.g. `"-bilbao"`.
#'
#' The discrete palettes are for categories. `"okabe"`, `"tol"`,
#' `"bright"`, `"vibrant"`, `"muted"`, `"light"`, and `"carto"` are
#' designed for color vision deficiencies. The sequential and diverging
#' palettes are for numbers: heat map grids and numeric tracks.
#'
#' @param kind   Which kinds to list: any of `"discrete"`, `"sequential"`,
#'        and `"diverging"`. `NULL` lists them all. Default: `NULL`
#'
#' @return A tibble with one row per palette: `name`, `kind`, and
#'         `colors` (a list column of hex codes).
#'
#' @seealso [fn_schema()]
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     list_palettes()
#'
#'     list_palettes("diverging")$name
#'
list_palettes <- function (kind = NULL) {

  validate_var_choices('kind', c("discrete", "sequential", "diverging"), max = Inf, null_ok = TRUE)

  pals <- names(PALETTES)
  if (!is.null(kind)) pals <- pals[PALETTE_KINDS[pals] %in% kind]

  tibble::tibble(
    name   = pals,
    kind   = unname(PALETTE_KINDS[pals]),
    colors = unname(lapply(PALETTES[pals], function (x) paste0("#", x))) )
}



#________________________________________________________
# See if a given string is a palette name.
#________________________________________________________
is_palette <- function (palette) {

  if (!is_scalar_character(palette) || is.na(palette)) return (FALSE)

  palette <- tolower(palette)
  if (eq(substr(palette, 1, 1), "-"))
    palette <- substr(palette, 2, nchar(palette))

  if (!is_string(palette, names(PALETTES))) return (FALSE)

  return (TRUE)
}



#________________________________________________________
# Return the colors from a specific palette.
# Reverse the order of colors when prefixed with a '-'.
#________________________________________________________
get_palette <- function (palette) {

  stopifnot(is_scalar_character(palette) && !is_na(palette))
  palette <- tolower(palette)

  if (eq(substr(palette, 1, 1), "-")) {
    palette <- substr(palette, 2, nchar(palette))

    stopifnot(is_string(palette, names(PALETTES)))
    return (rev(paste0("#", PALETTES[[palette]])))

  } else {
    stopifnot(is_string(palette, names(PALETTES)))
    return (paste0("#", PALETTES[[palette]]))
  }
}



#________________________________________________________
# Return a vector of discrete colors.
#________________________________________________________
get_n_colors <- function (n, palette = NULL) {

  stopifnot(is_scalar_integerish(n) && !is_na(n))
  if (is.logical(palette)) palette <- NULL

  if (is_null(palette)) {
    if (n <= 8)         { palette <- "okabe"
    } else if (n <= 12) { palette <- "tol"
    } else if (n <= 20) { palette <- "tableau20"
    } else              { palette <- "alphabet" }
  }

  if (is_palette(palette))
    return (rep_len(get_palette(palette), n))

  return (rep_len(palette, n))
}



#________________________________________________________
# Return a vector of discrete shapes.
#________________________________________________________
get_n_shapes <- function (n, values = NULL) {

  stopifnot(is_scalar_integerish(n) && !is_na(n))
  if (is.logical(values)) values <- NULL

  # Solid shapes first, then outlines, then letters.
  if (is_null(values))
    values <- if (n <= 15) c(16, 17, 15, 18, 3, 4, 8, 1, 2, 0, 5, 6, 7, 9, 10) else c(65:90, 97:122)

  return (rep_len(values, n))
}



#________________________________________________________
# Return a vector of discrete patterns.
#________________________________________________________
get_n_patterns <- function (n, values = NULL) {

  stopifnot(is_scalar_integerish(n) && !is_na(n))
  if (is.logical(values)) values <- NULL

  if (is_null(values))
    values <- seq_len(n)

  return (rep_len(values, n))
}



#________________________________________________________
# Resolve a `colors` specification to a vector of colors.
#
#   pal  - A palette name (optionally prefixed with "-" to reverse it),
#          a vector of colors, or NULL to pick a discrete palette that
#          suits `n` colors. A named vector is returned as-is.
#   n    - How many colors to return (the palette is recycled).
#   keys - Names for the returned colors; a factor contributes its
#          levels, a character vector its unique values.
#
# Example calls/returns:
#
# color_palette(pal = "okabe") or color_palette(n = 8)
# c("#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#F0E442", "#999999")
#
# color_palette(pal = "r4", keys = c("Male", "Female"))
# c(Male = "#DF536B", Female = "#61D04F")
#
# color_palette(n = 3)
# c("#0072B2", "#D55E00", "#009E73")
#
# color_palette(pal = "-hawaii")
# c("#B2F2FD", "#66E7D2", "#6BD38C", "#89BB48", "#9B951B", "#9A6D27", "#964C3E", "#922C54", "#8C0172")
#________________________________________________________
color_palette <- function (pal = NULL, n = NULL, keys = NULL) {

  if (!is_null(names(pal))) {
    pal <- unlist(pal[nzchar(names(pal))])
    if (!is.character(pal)) cli_abort("Named colors must be character strings, not {.type {pal}}.")
    return (pal)
  }

  if (is.factor(keys))         { keys <- levels(keys)
  } else if (is.character(keys)) { keys <- unique(keys)
  } else                         { keys <- NULL }

  if (is_null(n)) n <- length(keys)
  stopifnot(is_scalar_integerish(n) && !is_na(n) && n >= 0)

  colors <- if (is_null(pal)) {
    get_n_colors(max(n, 1))

  } else if (is_palette(pal)) {
    get_palette(pal)

  } else if (is.character(pal) && !anyNA(pal)) {
    pal

  } else {
    cli_abort("`pal` must be a palette name or a vector of colors, not {.type {pal}}.")
  }

  if (n > 0)            colors <- rep_len(colors, n)
  if (length(keys) > 0) colors <- setNames(colors, keys)

  return (colors)
}
