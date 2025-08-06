
library(ggplot2)
library(ggpattern)
library(magick)
library(magrittr)

ggplot() + 
  geom_polygon_pattern(
    mapping = aes(
      x = 1150 * sqrt(3)/2 * c(0, 1, 1, 0, -1, -1), 
      y = 1150 * .5 * c(2, 1, -1, -2, -1, 1) ),
    pattern          = 'image',
    pattern_type     = 'expand',
    pattern_filename = 'logo/bg.jpg',
    color            = '#000000',
    linewidth        = 5 ) + 
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(rect = element_rect(fill = 'transparent')) +
  annotate(
    geom   = 'text', 
    label  = 'rbiom',
    family = 'Milky Moringa',
    color  = 'white',
    size   = 36,
    x      = 0, 
    y      = 490 ) +
  annotation_raster(
    raster = image_read('logo/shadowed_cells.png'),
    xmin = -957, xmax = 957,
    ymin = -854, ymax = 314 )

ggsave(
  path     = 'logo',
  filename = 'rbiom.png', 
  device   = 'png',
  width    = 2200, 
  height   = 2200, 
  dpi      = 320,
  units    = 'px',
  bg       = 'transparent' )


# pkgdown website sets logo width to 120px
magick::image_read('logo/rbiom.png') %>%
  image_trim() %>%
  image_resize('120x') %>%
  image_write('man/figures/logo.png')


# height = 200px (150pt) for joss paper
magick::image_read('logo/rbiom.png') %>%
  image_trim() %>%
  image_resize('x200') %>%
  image_write('joss/figures/logo.png')
