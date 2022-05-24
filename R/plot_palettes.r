
get_palette <- function (pal, n) {
  if (pal == "color")   return (color_palette(n))
  if (pal == "shape")   return (shape_palette(n))
  if (pal == "pattern") return (pattern_palette(n))
}


color_palette <- function (n) {
  
  pal <- if (n <= 2) {
    c('#00B9EB', '#ED5F16')
  } else if (n <= 8) {
    # Color blind friendly palette by Bang Wong. https://doi.org/10.1038/nmeth.1618
    c('#0072B2', '#D55E00', '#009E73', '#E69F00', '#CC79A7', '#F0E442', '#56B4E9', '#000000')
  } else if (n <= 11) {
    c('#66CDAA', '#FF8C69', '#8DB6CD', '#008B8B', '#FF6EB4', '#A2CD5A', '#FF6347', '#FFC125', 
      '#EEC591', '#BEBEBE', '#CD96CD')
  } else if (n <= 20) {
    c('#8B4500', '#CD853F', '#FF1493', '#FFB5C5', '#CDC5BF', '#6C7B8B', '#B22222', '#FF0000', 
      '#FF7F24', '#FFD700', '#00FF00', '#2E8B57', '#00FFFF', '#63B8FF', '#0000FF', '#191970', 
      '#8B008B', '#A020F0', '#DA70D6', '#F08080')
  } else {
    c('#1C86EE', '#E31A1C', '#008B00', '#6A3D9A', '#A52A2A', '#FF7F00', '#FFD700', '#7EC0EE', 
      '#FB9A99', '#90EE90', '#CAB2D6', '#FDBF6F', '#B3B3B3', '#EEE685', '#B03060', '#FF83FA', 
      '#FF1493', '#0000FF', '#36648B', '#00CED1', '#00FF00', '#8B8B00', '#CDCD00', '#8B4500')
  }
  
  pal[((seq_len(n) - 1) %% (length(pal) + 1)) + 1]
}


shape_palette <- function (n) {
  
  pal <- if (n <= 6)  { c(16, 17, 15, 3, 7, 8)
  } else if (n <= 15) { c(0:14)
  } else              { c(65:90, 97:122) }
  
  pal[((seq_len(n) - 1) %% (length(pal) + 1)) + 1]
}


pattern_palette <- function (n) {
  pal <- unique(c(
    'bricks', 'fishscales', 'right45', 'horizonal_saw', 
    'hs_cross', 'crosshatch45', 'crosshatch',
    gridpattern::names_magick ))
  
  pal[((seq_len(n) - 1) %% (length(pal) + 1)) + 1]
}
