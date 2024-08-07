unlink("man", recursive = TRUE)
unlink("docs", recursive = TRUE)
devtools::document()
# Restart
# Install
genby  <- paste("rbiom", utils::packageVersion("rbiom"))
hmp50  <- as_rbiom(as.list(hmp50),  generated_by = genby)
babies <- as_rbiom(as.list(babies), generated_by = genby)
gems   <- as_rbiom(as.list(gems),   generated_by = genby)
usethis::use_data(hmp50, babies, gems, overwrite = TRUE)
# Restart
# Install
pkgdown::init_site()
cat(" ", file = "README.Rmd", append = TRUE)
devtools::build_readme()
pkgdown::build_site(preview = FALSE, devel = TRUE)
pkgdown::build_reference(topics = "rare_multiplot")
