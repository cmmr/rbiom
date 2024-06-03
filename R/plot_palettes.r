
PALETTES <- list(
  
  # Discrete color palettes from various sources.
  'bright'     = c("4477AA", "EE6677", "228833", "CCBB44", "66CCEE", "AA3377", "BBBBBB"),
  'vibrant'    = c("EE7733", "0077BB", "33BBEE", "EE3377", "CC3311", "009988", "BBBBBB"),
  'okabe'      = c("E69F00", "56B4E9", "009E73", "F0E442", "0072B2", "D55E00", "CC79A7", "999999"),
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
    if (n <= 10)        { palette <- "classic"
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
  
  if (is_null(values))
    values <- local({
      if (n <= 3)  return (15:17)
      if (n <= 15) return (c(0:14))
      return (c(65:90, 97:122))
    })
  
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







# get_palette <- function (pal, n) {
#   if (pal == "color")   return (color_palette(n = n))
#   if (pal == "shape")   return (shape_palette(n = n))
#   if (pal == "pattern") return (pattern_palette(n = n))
# }
# 
# 
# shape_palette <- function (n) {
#   
#   pal <- if (n <= 6)  { c(16, 17, 15, 3, 7, 8)
#   } else if (n <= 15) { c(0:14)
#   } else              { c(65:90, 97:122) }
#   
#   pal[((seq_len(n) - 1) %% (length(pal) + 1)) + 1]
# }
# 
# 
# pattern_palette <- function (n) {
#   pal <- c(
#     "bricks", "hexagons", "horizontalsaw", "hs_fdiagonal", 
#     "fishscales", "verticalsaw", "checkerboard", "octagons", 
#     "right45", "hs_cross", "hs_bdiagonal", "hs_diagcross", 
#     "hs_horizontal", "hs_vertical", "left45", "leftshingle", 
#     "rightshingle", "verticalbricks", "verticalleftshingle", 
#     "verticalrightshingle" )
#   
#   pal[((seq_len(n) - 1) %% (length(pal) + 1)) + 1]
# }



# Colorblind-friendly color sets and gradients.
# Hardcoded here to avoid additional dependencies on khroma and/or other R packages.

# Example calls/returns:
# 
# color_palette(pal = "okabe") or color_palette(n = 8)
# c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
#
# color_palette(pal = "r4", keys = c("Male", "Female"))
# c(Male = "#DF536B", Female = "#61D04F")
# 
# color_palette(n = 3)
# c("#E69F00", "#56B4E9", "#009E73")
# 
# color_palette(pal = "-ha")
# c("#B2F2FD", "#66E7D2", "#6BD38C", "#89BB48", "#9B951B", "#9A6D27", "#964C3E", "#922C54", "#8C0172")
# 
# color_palette(pal = "Orange")
# c("#FFFFFF", "Orange")

color_palette <- function (pal = NULL, n = NULL, keys = NULL) {
  
  if (!is_null(names(pal)))  return (unlist(pal[nzchar(names(pal))]))
  if (!is_null(names(keys))) return (unlist(keys[nzchar(names(keys))]))
  
  
  # Handle case where keys is an entire metadata column.
  if (!is_null(keys)) {
    
    if (is.factor(keys)) {
      keys <- levels(keys)
      n    <- length(keys)
      
    } else if (is.character(keys)) {
      if (any(duplicated(keys))) keys <- sort(unique(keys))
      n <- length(keys)
      
    } else {
      
      # Cycle through simple gradients for default use.
      if (is_null(pal)) pal <- 0
      if (is.numeric(pal))
        pal <- c('reds', 'oranges', 'greens', 'grays', 'purples')[(pal %% 5) + 1]
      keys <- NULL
    }
  }
  
  
  # Pick a discrete palette based on the number of unique colors required.
  if (!is.character(pal) && !is_null(n)) {
    
    if (is_null(pal)) {
      pal <- if (n <= 8)  { "okabe"
      } else if (n <= 12) { "tol"
      } else if (n <= 20) { "tableau20"
      } else              { "alphabet" } # 26
      
    } else if (is.numeric(pal)) {
      pals <- c(
        'okabe', 'carto', 'r4', 'polychrome', 'tol', 'bright', 'light', 'muted', 
        'vibrant', 'tableau', 'classic', 'alphabet', 'tableau20', 'kelly', 'fishy' )
      pals <- sapply(pals, function (p) { length(PALETTES[[p]]) - n })
      pals <- names(head(sort(pals[pals > -1])))
      if (length(pals) == 0) pals <- "alphabet"
      if (pal < 0)           pals <- paste0("-", pals)
      pal <- pals[(as.integer(pal) %% length(pals)) + 1]
      
    } else {
      stop("Invalid color_palette() `pal` argument: ", as.character(pal))
    }
  }
  
  
  # Look up this palette name (or partial name).
  if (length(pal) == 1) {
    
    pal_rev <- startsWith(pal, '-');
    pal     <- sub('-', '', pal)
    
    if (hasName(PALETTES, pal)) {
      pal <- paste0('#', PALETTES[[pal]])
    } else {
      pal_index <- pmatch(trimws(tolower(pal)), names(PALETTES))
      if (!is.na(pal_index)) { pal <- paste0('#', PALETTES[[pal_index]])
      } else                 { pal <- c("#FFFFFF", pal) }
    }
    
    if (pal_rev) pal <- rev(pal)
  }
    
  if (!is_null(n))    pal <- pal[((seq_len(n) - 1) %% (length(pal) + 1)) + 1]
  if (!is_null(keys)) pal <- setNames(pal, keys)
  
  return (pal)
}



