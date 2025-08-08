library(ggplot2)
library(ggpattern)
library(shadowtext)
library(magick)

ggplot() + 
  geom_polygon_pattern(
    mapping = aes(
      x = 1150 * sqrt(3)/2 * c(0, 1, 1, 0, -1, -1), 
      y = 1150 * .5 * c(2, 1, -1, -2, -1, 1) ),
    pattern          = 'image',
    pattern_type     = 'tile',
    pattern_filename = 'logo/bg.jpg',
    color            = '#0078a5',
    linewidth        = NA ) + 
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(rect = element_rect(fill = 'transparent')) +
  annotate(
    geom      = 'shadowtext',
    bg.colour = 'black',
    bg.r      = 0.03,
    label     = 'rbiom',
    family    = 'Boogie Boys',
    size      = 42, 
    y         = 600,
    color     = 'white',
    x         = 0 ) +
  annotation_raster(
    raster = image_read('logo/cartoon.png'),
    xmin = -957, xmax = 957,
    ymin = -804, ymax = 364 )



ggsave(
  path     = 'logo',
  filename = 'rbiom.png', 
  device   = 'png',
  width    = 2300, 
  height   = 2300, 
  dpi      = 400,
  units    = 'px',
  bg       = 'transparent' )


image_read('logo/rbiom.png') |>
  image_trim() |>
  image_resize('x200') |>
  image_write('man/figures/logo.png')
