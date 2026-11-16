
# Named `colors` / `shapes` / `patterns` vectors: names for groups that
# aren't drawn are ignored, and groups without a name get defaults.

scale_values <- function (p, aes = "colour") {
  for (s in p$scales$scales)
    if (aes %in% s$aesthetics) return (s$palette(0))
  NULL
}

expect_replays <- function (p) {
  env <- new.env(parent = globalenv())
  p2  <- suppressWarnings(eval(parse(text = p$code), envir = env))
  expect_equal(
    object   = ggplot2::ggplot_build(p2)$data,
    expected = ggplot2::ggplot_build(p)$data,
    ignore_attr = TRUE )
}


test_that("resolve_named_aes", {

  dflt <- get_n_colors(3)
  lvls <- c("A", "B", "C")

  # Complete, reordered.
  expect_identical(resolve_named_aes(c(C = "c", A = "a", B = "b"), lvls, dflt), c("a", "b", "c"))

  # Partial: defaults in level order, skipping colors already in use.
  expect_identical(resolve_named_aes(c(B = "red"), lvls, dflt), c(dflt[[1]], "red", dflt[[2]]))
  expect_identical(resolve_named_aes(c(B = dflt[[1]]), lvls, dflt), c(dflt[[2]], dflt[[1]], dflt[[3]]))
  expect_identical(resolve_named_aes(c(B = "#0072b2"), lvls, dflt), c(dflt[[2]], "#0072b2", dflt[[3]]))

  # Defaults are recycled when every one of them is in use.
  expect_identical(resolve_named_aes(c(A = "x"), lvls, c("x", "y"), what = "shapes"), c("x", "y", "x"))

  # Extra names: silent when every level is named, or when known.
  expect_silent(x <- resolve_named_aes(c(A = "a", B = "b", C = "c", D = "d"), lvls, dflt))
  expect_identical(x, c("a", "b", "c"))
  expect_silent(resolve_named_aes(c(A = "a", D = "d"), lvls, dflt, known = "D"))
  expect_warning(resolve_named_aes(c(A = "a", D = "d"), lvls, dflt), '"D"')
  expect_warning(resolve_named_aes(c(A = "a", D = "d"), lvls, dflt, by = "Field"), "Field")

  # Abbreviations, but not for names that are known to be something else.
  lvls <- c("Female", "Male", "Female vs Male")
  expect_identical(resolve_named_aes(c(m = "b", f = "r", "female vs" = "p"), lvls, dflt), c("r", "b", "p"))
  expect_identical(
    object   = resolve_named_aes(c(Female = "r", Male = "b"), "Female vs Male", dflt, known = c("Female", "Male")),
    expected = dflt[[1]] )

  # Pairs: as labelled, or in the other order.
  expect_identical(resolve_named_aes(c("Male vs Female" = "p"), lvls, dflt)[[3]], "p")
  expect_silent(resolve_named_aes(c("Male vs Other" = "p"), lvls, dflt, known = c("Male", "Other")))

  # Mixed: unnamed entries go to the unnamed levels, in order, then defaults.
  lvls <- c("A", "B", "C", "D")
  expect_identical(resolve_named_aes(c(C = "c", "x", "y"), lvls, dflt), c("x", "y", "c", dflt[[1]]))
  expect_identical(resolve_named_aes(c(C = 3, 9), lvls, c(1, 2, 3), what = "shapes"), c(9, 1, 3, 2))

  expect_error(resolve_named_aes(c(A = "a", A = "b"), lvls, dflt), "Duplicated")

  expect_identical(shape_names(c(16, 4, 65)), c("circle small", "cross", "A"))
  expect_identical(shape_names("cross"), "cross")
})


test_that("named colors for groups that aren't plotted", {

  sites <- c(Saliva = "red", Stool = "brown", `Buccal mucosa` = "pink", `Mid vagina` = "purple", `Anterior nares` = "green")

  # A project-wide map, after filtering samples.
  sub <- subset(rare50, `Body Site` %in% c("Saliva", "Buccal mucosa"))
  expect_silent(p <- adiv_boxplot(sub, x = "Body Site", colors = sites, test = "none"))
  expect_identical(scale_values(p), c("pink", "red"))
  expect_match(p$code, 'scale_color_manual(values = c("pink", "red"))', fixed = TRUE)
  expect_replays(p)

  # Unnamed vectors, palettes, and TRUE are unchanged.
  expect_identical(scale_values(adiv_boxplot(sub, x = "Body Site", colors = c("red", "blue"), test = "none")), c("red", "blue"))
  expect_identical(scale_values(adiv_boxplot(sub, x = "Body Site", colors = "bright", test = "none")), get_palette("bright")[1:2])
  expect_identical(scale_values(adiv_boxplot(sub, x = "Body Site", test = "none")), get_n_colors(2))
  expect_identical(scale_values(adiv_boxplot(rare50, x = "Body Site", colors = "red", test = "none")), rep("red", 5))
  expect_null(scale_values(adiv_boxplot(sub, x = "Body Site", colors = FALSE, test = "none")))

  # A partial map.
  expect_silent(p <- adiv_boxplot(rare50, x = "Sex", colors = c(Male = "red"), test = "none"))
  expect_identical(scale_values(p), c("#0072B2", "red"))
  expect_identical(scale_values(p, "fill"), c("#0072B2", "red"))

  # A name that matches nothing, with a group left unnamed: likely a typo.
  expect_warning(p <- adiv_boxplot(rare50, x = "Sex", colors = c(Male = "red", Femael = "blue"), test = "none"), "Femael")
  expect_identical(scale_values(p), c("#0072B2", "red"))

  # Other plot types go through the same code.
  expect_silent(p <- adiv_corrplot(sub, x = "Age", stat.by = "Body Site", colors = sites, test = "none"))
  expect_identical(scale_values(p), c("pink", "red"))
  expect_silent(p <- bdiv_ord_plot(sub, stat.by = "Body Site", colors = sites, shapes = c(Stool = 4, Saliva = 16, `Buccal mucosa` = 3), test = "none"))
  expect_identical(scale_values(p), c("pink", "red"))
  expect_identical(scale_values(p, "shape"), c(3, 16))

  # Taxa: names for taxa that aren't among the most abundant are ignored.
  expect_silent(p <- taxa_stacked(rare50, rank = "Phylum", taxa = 3, colors = c(Firmicutes = "red", Fusobacteria = "blue")))
  expect_identical(scale_values(p, "fill")[levels(p$data$.taxa) == "Firmicutes"], "red")
  expect_false("blue" %in% scale_values(p, "fill"))
})


test_that("named colors for between-group distances", {

  sex <- c(Female = "red", Male = "blue")

  # Only the mixed pair is drawn: the group colors don't apply.
  expect_silent(p <- bdiv_boxplot(rare50, x = "!=Sex", stat.by = "Sex", colors = sex, test = "none"))
  expect_identical(levels(p$data$Sex), "Female vs Male")
  expect_identical(scale_values(p), "#0072B2")

  # Within-group levels use the groups' colors; the pair gets a default.
  expect_silent(p <- bdiv_boxplot(rare50, x = "Sex", stat.by = "Sex", colors = sex, test = "none"))
  expect_identical(levels(p$data$Sex), c("Female", "Male", "Female vs Male"))
  expect_identical(scale_values(p), c("red", "blue", "#0072B2"))

  # A color for the pair itself, in either order.
  for (pair in c("Female vs Male", "Male vs Female")) {
    p <- bdiv_boxplot(rare50, x = "Sex", stat.by = "Sex", colors = c(sex, setNames("purple", pair)), test = "none")
    expect_identical(scale_values(p), c("red", "blue", "purple"))
  }

  expect_replays(p)
})


test_that("named shapes and patterns", {

  # A partial map; defaults skip the shapes in use.
  expect_silent(p <- adiv_boxplot(rare50, x = "Sex", layers = "bd", shapes = c(Male = 16), test = "none"))
  expect_identical(scale_values(p, "shape"), c(17, 16))

  # Named by string: the defaults are strings too.
  expect_silent(p <- adiv_boxplot(rare50, x = "Sex", layers = "bd", shapes = c(Male = "cross"), test = "none"))
  expect_identical(scale_values(p, "shape"), c("circle small", "cross"))
  expect_s3_class(ggplot2::ggplotGrob(p), 'gtable')

  # Extra names.
  expect_silent(p <- adiv_boxplot(rare50, x = "Sex", layers = "bd", shapes = c(Male = 4, Female = 3, Other = 2), test = "none"))
  expect_identical(scale_values(p, "shape"), c(3, 4))
  expect_warning(adiv_boxplot(rare50, x = "Sex", layers = "bd", shapes = c(Male = 4, Other = 2), test = "none"), "Other")

  # Unnamed and TRUE are unchanged.
  expect_identical(scale_values(adiv_boxplot(rare50, x = "Sex", layers = "bd", shapes = c(3, 4), test = "none"), "shape"), c(3, 4))
  expect_identical(scale_values(adiv_boxplot(rare50, x = "Sex", layers = "bd", test = "none"), "shape"), get_n_shapes(2))

  skip_if_not_installed('fillpattern')

  expect_silent(p <- adiv_boxplot(rare50, x = "Sex", layers = "b", patterns = c(Male = "brick"), test = "none"))
  expect_match(p$code, 'scale_fill_pattern(patterns = c("chevron", "brick"))', fixed = TRUE)
  expect_s3_class(ggplot2::ggplotGrob(p), 'gtable')

  expect_silent(p <- adiv_boxplot(rare50, x = "Sex", layers = "b", patterns = c(Male = "brick", Female = "grid", Other = "fish"), test = "none"))
  expect_match(p$code, 'scale_fill_pattern(patterns = c("grid", "brick"))', fixed = TRUE)
  expect_warning(adiv_boxplot(rare50, x = "Sex", layers = "b", patterns = c(Male = "brick", Other = "fish"), test = "none"), "Other")

  p <- adiv_boxplot(rare50, x = "Sex", layers = "b", patterns = c("brick", "grid"), test = "none")
  expect_match(p$code, 'scale_fill_pattern(patterns = c("brick", "grid"))', fixed = TRUE)
})


test_that("named colors for heat map tracks", {

  mtx <- matrix(1:20, nrow = 4, dimnames = list(letters[1:4], LETTERS[1:5]))
  sex <- c("Male", "Female")[c(1, 2, 1, 2, 1)]

  # Extra names are ignored.
  expect_silent(p <- plot_heatmap(mtx, tracks = list(Sex = list(values = sex, colors = c(Male = "blue", Female = "red", Other = "gray")))))
  expect_match(p$code, 'values = c("red", "blue")', fixed = TRUE)

  # Values without a name get default colors.
  expect_silent(p <- plot_heatmap(mtx, tracks = list(Sex = list(values = sex, colors = c(Male = "blue")))))
  expect_match(p$code, 'values = c("#0072B2", "blue")', fixed = TRUE)

  expect_warning(plot_heatmap(mtx, tracks = list(Sex = list(values = sex, colors = c(Male = "blue", Zzz = "red")))), "Zzz")

  # Tracks from an rbiom object's metadata; a bare named vector too.
  sites <- c(Saliva = "red", Stool = "brown", `Buccal mucosa` = "pink", `Mid vagina` = "purple", `Anterior nares` = "green")
  sub   <- subset(rare50, `Body Site` %in% c("Saliva", "Buccal mucosa"))

  expect_silent(p <- taxa_heatmap(sub, rank = "Phylum", tracks = list('Body Site' = sites)))
  expect_match(p$code, 'values = c("pink", "red")', fixed = TRUE)

  expect_silent(p <- bdiv_heatmap(sub, tracks = list('Body Site' = list(colors = sites), Sex = list(colors = c(m = "blue")))))
  expect_match(p$code, 'values = c("pink", "red")', fixed = TRUE)
  expect_match(p$code, 'values = c("#0072B2", "blue")', fixed = TRUE)
  expect_replays(p)
})


test_that("named colors for taxa_forest", {

  biom <- subset(hmp50, `Body Site` %in% c("Saliva", "Stool"))
  cols <- function (p) scale_values(p)[match(c("Higher in Saliva", "Higher in Stool", "Not significant"), levels(p$data$.change))]

  # By group, with names for groups that aren't compared.
  expect_silent(p <- taxa_forest(biom, stat.by = "Body Site", rank = "Genus", p.top = Inf, colors = c(Saliva = "red", Stool = "brown", `Buccal mucosa` = "pink")))
  expect_identical(cols(p), c("red", "brown", "#999999"))
  expect_replays(p)

  # By legend label or direction.
  p <- taxa_forest(biom, stat.by = "Body Site", rank = "Genus", p.top = Inf, colors = c(`Higher in Saliva` = "red", `not significant` = "black"))
  expect_identical(cols(p)[c(1, 3)], c("red", "black"))
  expect_true(cols(p)[[2]] %in% c("#D55E00", "#0072B2"))

  p1 <- taxa_forest(biom, stat.by = "Body Site", rank = "Genus", p.top = Inf, colors = c(Decreased = "red", Increased = "blue"))
  p2 <- taxa_forest(biom, stat.by = "Body Site", rank = "Genus", p.top = Inf, colors = c("blue", "red"))
  expect_identical(scale_values(p1), scale_values(p2))

  expect_warning(taxa_forest(biom, stat.by = "Body Site", rank = "Genus", colors = c(up = "red", Stool = "brown")), '"up"')
  expect_error(taxa_forest(biom, stat.by = "Body Site", rank = "Genus", colors = c(Saliva = "red", "blue")), "named")
})
