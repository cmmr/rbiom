
library(ggplot2)
library(ggpattern)
library(shadowtext)
library(ggstar)
library(magick)

ggplot() + 
  geom_polygon_pattern(
    mapping = aes(
      x = 1700 * sqrt(3)/2 * c(0, 1, 1, 0, -1, -1), 
      y = 1700 * .5 * c(2, 1, -1, -2, -1, 1) ),
    pattern          = 'image',
    pattern_type     = 'expand',
    pattern_filename = 'logo/bg.jpg',
    color            = '#0078a5',
    linewidth        = 6 ) + 
  geom_segment( # horizontal grid lines
    linewidth = 1,
    alpha     = 0.6,
    mapping   = aes(
      x    = rep(265 * -5, 5),
      y    = 250 - 265 * 0:4,
      xend = rep(265 * 5, 5),
      yend = 250 - 265 * 0:4 )) + 
  geom_segment( # vertical grid lines
    linewidth = 1,
    alpha     = 0.6,
    mapping   = aes(
      x    =  265 * -5:5,
      y    =  rep(250, 11),
      xend =  265 * -5:5,
      yend = rep(-810, 11) )) +
  geom_path( # significance line
    linewidth = 2.5,
    color     = 'black',
    mapping   = aes(
      x = 265 * c(2.75, 3.5, 3.5, 2.75), 
      y = 250 - 265 * c(1, 1, 3, 3) ) ) +
  annotate(
    geom  = 'star',
    color = 'black',
    fill  = 'red',
    size  = 9,
    starstroke = 2,
    x     = 265 * 4.25,
    y     = 250 - 265 * 2 ) +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(rect = element_rect(fill = 'transparent')) +
  annotate(
    geom      = 'shadowtext',
    label     = 'rbiom',
    family    = 'Milky Moringa',
    color     = 'white',
    bg.colour = 'black',
    bg.r      = 0.03,
    size      = 28,
    x         = 0,
    y         = 600 ) +
  annotation_raster(
    raster = image_read('logo/original.png'),
    xmin = -1307, xmax = 607,
    ymin = -854, ymax = 314 )



ggsave(
  path     = 'logo',
  filename = 'rbiom.png', 
  device   = 'png',
  width    = 2400, 
  height   = 2400, 
  dpi      = 400,
  units    = 'px',
  bg       = 'transparent' )


image_read('logo/rbiom.png') |>
  image_trim() |>
  image_resize('x200') |>
  image_write('man/figures/logo.png')
