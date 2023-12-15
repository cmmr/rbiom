# Restart
unlink("man", recursive = TRUE)
unlink("docs", recursive = TRUE)
pkgdown::init_site()
devtools::build_readme()
devtools::document()
# Install
pkgdown::build_site()
