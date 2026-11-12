#________________________________________________________
# The argument registry.
#
# One description of every user-facing argument of the functions listed
# in FN_SCHEMA. It is the single source of truth for:
#
#   * the option lists and numeric limits that the validators use
#     (validate_var_choices(), validate_var_range(), validate_layers(),
#     validate_ord(), ... read them from CHOICES and ARGS rather than
#     repeating them);
#   * schema_check(), which each covered function runs on its arguments
#     on entry, and which fn_check() runs against a summary of a dataset;
#   * fn_schema(), the exported, JSON-ready description.
#
# A spec is a plain list. `type` is required; every other field is
# optional and only present when it applies.
#
#   type      biom, tree, file, dots, logical, integer, number, string,
#             choice, layers, field, rank, taxa, metric, aesthetic,
#             palette, tracks
#   r_only    Only an R caller may set it (an R object, a file path, or a
#             host setting). fn_check() refuses it.
#   deprecated  The version that deprecated it.
#   class     For fields: "categorical", "numeric", or "any".
#   coerce    For categorical fields: numeric fields are accepted too, and
#             treated as categories.
#   ops       For fields: the accepted prefixes ("==", "!=").
#   extra     For fields: values accepted besides the dataset's fields
#             (e.g. ".taxa").
#   div       For metrics: "alpha" or "beta".
#   aes, by   For aesthetics: what is mapped (color, shape, pattern) and
#             which argument's values are mapped.
#   forms     The forms accepted, when there are several (see fn_schema()).
#   palettes  Palette kinds suited to the argument (see list_palettes()).
#   needs     An optional package the argument requires.
#   choices   Fixed options. Named character vector for layers.
#   choices_when  Choices that depend on other arguments: a list of
#             `list(when = <conditions>, choices = <chr>)`; the first
#             whose conditions hold applies, else `choices`.
#   min, max  Numeric limits.
#   max_n     How many values may be given (absent: any number).
#   per_side  One value, or two: for rows and for columns.
#   null      Present when NULL is accepted: what NULL means.
#   default_note  What an expression default (e.g. `stat.by = x`) means.
#   when      Conditions under which the argument has any effect.
#   note      Anything else a caller needs to know.
#
# Internal-only fields, not exported by fn_schema():
#
#   matcher   "match.arg" when the function checks the value with
#             match.arg(tolower(x), choices) instead of
#             validate_var_choices().
#   check     Name of a custom check in SCHEMA_CHECKS.
#   all       A value that selects every choice (e.g. ".all").
#________________________________________________________


#________________________________________________________
# A spec, dropping NULL fields.
#________________________________________________________
arg <- function (type, ...) {
  x <- list(type = type, ...)
  x[!vapply(x, is.null, logical(1))]
}


#________________________________________________________
# A copy of `spec` with fields changed. NULL removes a field.
#________________________________________________________
like <- function (spec, ...) {
  mods <- list(...)
  for (i in names(mods)) spec[i] <- list(mods[[i]])
  spec[!vapply(spec, is.null, logical(1))]
}


#________________________________________________________
# Conditions for a spec's `when` field. A list of them must all hold
# for the argument to have any effect. Values are canonical: layer
# names (not letters) and metric ids.
#________________________________________________________
when_is      <- function (param, ...) list(param = param, is      = c(...))
when_is_not  <- function (param, ...) list(param = param, is_not  = c(...))
when_has_any <- function (param, ...) list(param = param, has_any = c(...))
when_null    <- function (param)      list(param = param, is_null  = TRUE)
when_set     <- function (param)      list(param = param, not_null = TRUE)


#________________________________________________________
# Range of a numeric spec, for validate_var_range().
#________________________________________________________
arg_range <- function (spec) {
  c(spec[["min"]] %||% -Inf, spec[["max"]] %||% Inf)
}




#________________________________________________________
# Option lists. The validators and the specs below both use these.
#________________________________________________________
CHOICES <- list(

  ci             = c('ci', 'range', 'sd', 'se', 'mad'),
  p.adj          = stats::p.adjust.methods,
  xlab.angle     = c('auto', '0', '30', '90'),
  fit            = c('lm', 'log', 'gam'),
  alt            = c("!=", ">", "<"),

  adiv.transform = c("none", "rank", "log", "log1p", "sqrt"),
  transform      = c("none", "rank", "log", "log1p", "sqrt", "percent"),
  ties           = c("average", "first", "last", "random", "max", "min"),
  norm           = c("none", "percent", "chord", "binary", "clr", "rclr"),

  unc            = c("asis", "singly", "grouped", "drop"),
  ord            = c("PCoA", "tSNE", "NMDS", "UMAP"),
  dist           = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
  clust          = c("average", "ward.D", "ward.D2", "mcquitty", "single", "median", "complete", "centroid"),

  box.test       = c('auto', 'none', 'wilcox', 'emmeans', 'kruskal'),
  box.bdiv.test  = c('auto', 'none', 'adonis2', 'mrpp', 'betadisper'),
  corr.test      = c('none', 'emmeans', 'emtrends'),
  corr.bdiv.test = c('none', 'mantel'),
  ord.test       = c("adonis2", "mrpp", "betadisper", "none"),
  stats.test     = c("wilcox", "kruskal", "emmeans", "emtrends"),
  bdiv.test      = c("adonis2", "mrpp", "betadisper"),
  bdiv.regr.test = c("mantel", "adonis2"),
  distmat.test   = c("adonis2", "mrpp", "betadisper", "mantel", "none"),
  da.method      = "linda",

  box.y.transform  = c('none', 'log1p', 'sqrt'),
  rare.y.transform = c('none', 'log10'),

  rescale        = c("none", "rows", "cols"),
  legend         = c("right", "bottom"),
  heatmap.label  = c("rows", "cols", "both", "bottom", "right", "none"),
  heatmap.trees  = c("rows", "cols", "both", "left", "top", "none"),

  biom.format    = c("tab", "json", "hdf5"),

  box.layers     = c(
    'x' = "box",   'v' = "violin", 'e' = "errorbar",
    'b' = "bar",   's' = "strip",  'l' = "linerange",
    'd' = "dot",   'c' = "crossbar", 'p' = "pointrange" ),
  corr.layers    = c(
    't' = "trend", 'c' = "confidence", 'p' = "point",
    'n' = "name",  'r' = "residual" ),
  corr.bdiv.layers = c(                # no sample names on sample pairs
    't' = "trend", 'c' = "confidence", 'p' = "point",
    'r' = "residual" ),
  ord.layers     = c(
    'p' = "point",   'n' = "name",  's' = "spider",
    'd' = "density", 't' = "taxon", 'm' = "mean",
    'a' = "arrow",   'e' = "ellipse" )
)




#________________________________________________________
# Specs shared by several functions.
#________________________________________________________
ARGS <- list(

  # Objects only an R caller can supply, and host settings.
  biom      = arg('biom', r_only = TRUE),
  tree      = arg('tree', r_only = TRUE, null = "Use the tree in `biom`."),
  file      = arg('file', r_only = TRUE),
  overwrite = arg('logical', r_only = TRUE, max_n = 1),
  clone     = arg('logical', r_only = TRUE, max_n = 1),
  cpus      = arg('integer', r_only = TRUE, min = 1, max = 256, max_n = 1, null = "All available cores.", default_note = "All available cores."),
  dots      = arg('dots', r_only = TRUE),

  # Diversity metrics.
  adiv       = arg('metric', div = 'alpha'),
  bdiv       = arg('metric', div = 'beta'),
  bdiv1      = arg('metric', div = 'beta', max_n = 1),
  weighted   = arg('logical', deprecated = "2.3.0", null = "Not set: the metric name decides.", note = "Use a metric name such as \"weighted_unifrac\" instead."),
  normalized = arg('logical', deprecated = "2.3.0", null = "Not set: the metric name decides.", note = "Use a metric name such as \"normalized_unifrac\" instead."),
  alpha      = arg('number', min = 0, max = 1, max_n = 1, when = list(when_has_any('bdiv', 'generalized_unifrac'))),

  # Taxa.
  rank      = arg('rank'),
  rank1     = arg('rank', max_n = 1),
  taxa      = arg('taxa', forms = c('count', 'fraction', 'names')),
  unc       = arg('choice', choices = CHOICES$unc, max_n = 1),
  lineage   = arg('logical', max_n = 1),
  other     = arg('logical', forms = c('logical', 'string'), max_n = 1, null = "No \"Other\" taxon.", check = 'check_other', note = "A string names the \"Other\" taxon."),

  # Transformations and randomness.
  adiv.transform = arg('choice', choices = CHOICES$adiv.transform, max_n = 1),
  transform      = arg('choice', choices = CHOICES$transform, max_n = 1),
  ties           = arg('choice', choices = CHOICES$ties, max_n = 1, when = list(when_is('transform', 'rank'))),
  seed           = arg('integer', min = 0, max_n = 1),
  permutations   = arg('integer', min = 1, max_n = 1),

  # Statistics.
  p.adj     = arg('choice', choices = CHOICES$p.adj, max_n = 1),
  level     = arg('number', min = 0.5, max = 1, max_n = 1),
  p.label   = arg('number', min = 0, max = 1, note = "Several values give the breaks for asterisks."),
  p.top     = arg('number', min = 0, max_n = 1),
  fit       = arg('choice', choices = CHOICES$fit, max_n = 1),
  at        = arg('number', null = "The median of `x`."),
  alt       = arg('choice', choices = CHOICES$alt, max_n = 1),
  mu        = arg('number', max_n = 1),

  # Plot options.
  caption    = arg('logical', max_n = 1),
  check      = arg('logical', max_n = 1),
  xlab.angle = arg('choice', choices = CHOICES$xlab.angle, max_n = 1, check = 'check_xlab_angle', note = "The numbers 0, 30, and 90 are accepted too."),
  colors     = arg('aesthetic', aes = 'color', by = 'stat.by', forms = c('logical', 'palette', 'values', 'named_values'), palettes = 'discrete', null = "No colors."),
  shapes     = arg('aesthetic', aes = 'shape', by = 'stat.by', forms = c('logical', 'values', 'named_values'), null = "No shapes."),
  patterns   = arg('aesthetic', aes = 'pattern', by = 'stat.by', forms = c('logical', 'values', 'named_values'), needs = 'fillpattern', null = "No patterns."),
  facet.by   = arg('field', class = 'categorical', null = "No facets."),
  strata     = arg('field', class = 'categorical', max_n = 1, null = "No blocking."),
  within     = arg('field', class = 'categorical', null = "No restriction."),
  between    = arg('field', class = 'categorical', null = "No restriction."),
  rline      = arg('logical', forms = c('logical', 'integer'), max_n = 1, check = 'check_rline', note = "TRUE: an automatically selected depth. A number: that depth. FALSE: no line.")
)




#________________________________________________________
# Box plots.
#________________________________________________________
BOX <- list(
  x        = arg('field', class = 'categorical', coerce = TRUE, max_n = 1, null = "All samples in one group."),
  layers   = arg('layers', choices = CHOICES$box.layers),
  stat.by  = arg('field', class = 'categorical', coerce = TRUE, max_n = 1, null = "No statistics, and no grouping by color, shape, or pattern.", default_note = "Same as `x`."),
  facet.by = like(ARGS$facet.by, coerce = TRUE),
  shapes   = like(ARGS$shapes,   when = list(when_has_any('layers', 'dot', 'strip', 'pointrange'))),
  patterns = like(ARGS$patterns, when = list(when_has_any('layers', 'box', 'bar', 'violin'))),
  test     = arg('choice', choices = CHOICES$box.test, max_n = 1),
  flip     = arg('logical', max_n = 1),
  stripe   = arg('logical', max_n = 1, null = "Same as `flip`."),
  ci       = arg('choice', choices = CHOICES$ci, max_n = 1, when = list(when_has_any('layers', 'crossbar', 'errorbar', 'linerange', 'pointrange'))),
  outliers = arg('logical', max_n = 1, null = "Hidden when a dot or strip layer is drawn, shown otherwise.", when = list(when_has_any('layers', 'box'))),
  strata   = like(ARGS$strata, coerce = TRUE)
)

BOX$bdiv <- list(
  x        = like(BOX$x,        ops = c('==', '!=')),
  stat.by  = like(BOX$stat.by,  ops = c('==', '!=')),
  facet.by = like(BOX$facet.by, ops = c('==', '!=')),
  test     = arg('choice', choices = CHOICES$box.bdiv.test, max_n = 1),
  permutations = like(ARGS$permutations, when = list(
    when_is('test', 'auto', 'adonis2', 'mrpp', 'betadisper'),
    when_set('stat.by') ))
)




#________________________________________________________
# Scatter plots with trendlines.
#________________________________________________________
CORR <- list(
  x        = arg('field', class = 'numeric', max_n = 1),
  layers   = arg('layers', choices = CHOICES$corr.layers),
  stat.by  = arg('field', class = 'categorical', coerce = TRUE, max_n = 1, null = "No groups."),
  facet.by = like(ARGS$facet.by, coerce = TRUE),
  shapes   = like(ARGS$shapes, when = list(when_has_any('layers', 'point'))),
  test     = arg('choice', choices = CHOICES$corr.test, max_n = 1),
  strata   = like(ARGS$strata, coerce = TRUE)
)




#________________________________________________________
# Heat maps.
#________________________________________________________
HEAT <- list(
  tracks       = arg('tracks', class = 'any', forms = c('fields', 'spec'), null = "No tracks.", check = 'check_tracks'),
  grid         = arg('palette', forms = c('palette', 'values', 'spec'), palettes = c('sequential', 'diverging')),
  label        = arg('logical', per_side = TRUE, forms = c('logical', 'choice'), choices = CHOICES$heatmap.label, check = 'check_sides'),
  label_size   = arg('number', per_side = TRUE, null = "Sized to fit the grid.", check = 'check_sides'),
  rescale      = arg('choice', choices = CHOICES$rescale, max_n = 1),
  clust        = arg('choice', choices = CHOICES$clust, per_side = TRUE, forms = c('choice', 'logical'), check = 'check_clust', note = "FALSE (or NA) leaves the rows or columns in their original order."),
  trees        = arg('logical', per_side = TRUE, forms = c('logical', 'choice'), choices = CHOICES$heatmap.trees, check = 'check_sides'),
  dist         = arg('choice', choices = CHOICES$dist, per_side = TRUE, check = 'check_dist'),
  asp          = arg('number', max_n = 1),
  tree_height  = arg('number', per_side = TRUE, check = 'check_sides'),
  track_height = arg('number', per_side = TRUE, check = 'check_sides'),
  legend       = arg('choice', choices = CHOICES$legend, max_n = 1)
)




#________________________________________________________
# Statistics tables.
#________________________________________________________
STATS <- list(
  regr       = arg('field', class = 'numeric', max_n = 1, null = "No regression.", when = list(when_is('test', 'emmeans', 'emtrends'))),
  stat.by    = arg('field', class = 'categorical', max_n = 1, null = "No groups."),
  split.by   = arg('field', class = 'categorical', null = "No split."),
  covariates = arg('field', class = 'any', null = "No covariates."),
  test       = arg('choice', choices = CHOICES$stats.test, max_n = 1),
  fit        = like(ARGS$fit, when = list(when_is('test', 'emmeans', 'emtrends'))),
  at         = like(ARGS$at, null = "The median of `regr`.", when = list(when_set('regr')))
)

STATS$rules <- list(
  list(when = list(when_is('test', 'wilcox', 'kruskal')), forbids  = c('strata', 'covariates')),
  list(when = list(when_is('test', 'kruskal')),           requires = 'stat.by'),
  list(when = list(when_is('test', 'emtrends')),          requires = 'regr') )




#________________________________________________________
# Return values.
#________________________________________________________
RETURNS <- list(
  plot_stats = list(type = "plot", has = c("data", "code", "stats")),
  plot       = list(type = "plot", has = c("data", "code")),
  table      = list(type = "table", has = c("code")),
  factor     = list(type = "factor"),
  matrix     = list(type = "matrix", has = c("code")),
  dist       = list(type = "dist", has = c("code")),
  biom       = list(type = "biom", has = c("code")),
  file       = list(type = "file")
)




#________________________________________________________
# One entry per covered function, with its arguments in the order of
# formals(). To cover another function, add an entry here and call
# schema_check() on entry to the function.
#
# `rules` are constraints between arguments, exported as they are:
#   one_of    Exactly one of these arguments must be set (non-NULL).
#   any_of    At least one of these arguments must be set.
#   when + requires / forbids
#             While the `when` conditions hold, these arguments must be
#             set (requires) or must be NULL (forbids).
#
# Internal fields of an entry:
#   check     The check in SCHEMA_CHECKS that enforces the rules, with
#             the function's own messages.
#   drop_ops  `==` / `!=` prefixes are removed from fields, rather than
#             moved to `within` / `between`.
#________________________________________________________
FN_SCHEMA <- list(


  #________________________________________________________
  # Alpha diversity.
  #________________________________________________________

  adiv_boxplot = list(
    returns = RETURNS$plot_stats,
    params  = list(
      biom       = ARGS$biom,
      x          = BOX$x,
      adiv       = ARGS$adiv,
      layers     = BOX$layers,
      stat.by    = BOX$stat.by,
      facet.by   = BOX$facet.by,
      colors     = ARGS$colors,
      shapes     = BOX$shapes,
      patterns   = BOX$patterns,
      test       = BOX$test,
      flip       = BOX$flip,
      stripe     = BOX$stripe,
      ci         = BOX$ci,
      level      = ARGS$level,
      p.adj      = ARGS$p.adj,
      outliers   = BOX$outliers,
      xlab.angle = ARGS$xlab.angle,
      p.label    = ARGS$p.label,
      transform  = ARGS$adiv.transform,
      strata     = BOX$strata,
      caption    = ARGS$caption,
      '...'      = ARGS$dots )),

  adiv_corrplot = list(
    returns = RETURNS$plot_stats,
    params  = list(
      biom       = ARGS$biom,
      x          = CORR$x,
      adiv       = ARGS$adiv,
      layers     = CORR$layers,
      stat.by    = CORR$stat.by,
      facet.by   = CORR$facet.by,
      colors     = ARGS$colors,
      shapes     = CORR$shapes,
      test       = CORR$test,
      fit        = ARGS$fit,
      at         = ARGS$at,
      level      = ARGS$level,
      p.adj      = ARGS$p.adj,
      transform  = ARGS$adiv.transform,
      alt        = ARGS$alt,
      mu         = ARGS$mu,
      strata     = CORR$strata,
      caption    = ARGS$caption,
      check      = ARGS$check,
      '...'      = ARGS$dots )),

  adiv_stats = list(
    returns = RETURNS$table,
    rules   = STATS$rules,
    check   = 'check_stats_test',
    params  = list(
      biom       = ARGS$biom,
      regr       = STATS$regr,
      stat.by    = STATS$stat.by,
      adiv       = ARGS$adiv,
      split.by   = STATS$split.by,
      strata     = ARGS$strata,
      covariates = STATS$covariates,
      transform  = ARGS$adiv.transform,
      test       = STATS$test,
      fit        = STATS$fit,
      at         = STATS$at,
      level      = ARGS$level,
      alt        = ARGS$alt,
      mu         = ARGS$mu,
      p.adj      = ARGS$p.adj )),

  adiv_matrix = list(
    returns = RETURNS$matrix,
    params  = list(
      biom       = ARGS$biom,
      adiv       = ARGS$adiv,
      tree       = ARGS$tree,
      transform  = ARGS$adiv.transform,
      ties       = ARGS$ties,
      seed       = ARGS$seed,
      cpus       = ARGS$cpus )),


  #________________________________________________________
  # Beta diversity.
  #________________________________________________________

  bdiv_boxplot = list(
    returns = RETURNS$plot_stats,
    params  = list(
      biom         = ARGS$biom,
      x            = BOX$bdiv$x,
      bdiv         = ARGS$bdiv,
      layers       = BOX$layers,
      weighted     = ARGS$weighted,
      tree         = ARGS$tree,
      within       = ARGS$within,
      between      = ARGS$between,
      stat.by      = BOX$bdiv$stat.by,
      facet.by     = BOX$bdiv$facet.by,
      colors       = ARGS$colors,
      shapes       = BOX$shapes,
      patterns     = BOX$patterns,
      test         = BOX$bdiv$test,
      flip         = BOX$flip,
      stripe       = BOX$stripe,
      ci           = BOX$ci,
      level        = ARGS$level,
      p.adj        = ARGS$p.adj,
      outliers     = BOX$outliers,
      xlab.angle   = ARGS$xlab.angle,
      p.label      = ARGS$p.label,
      transform    = ARGS$transform,
      strata       = ARGS$strata,
      seed         = ARGS$seed,
      permutations = BOX$bdiv$permutations,
      caption      = ARGS$caption,
      alpha        = ARGS$alpha,
      cpus         = ARGS$cpus,
      '...'        = ARGS$dots )),

  bdiv_corrplot = list(
    returns = RETURNS$plot_stats,
    params  = list(
      biom         = ARGS$biom,
      x            = CORR$x,
      bdiv         = ARGS$bdiv,
      layers       = arg('layers', choices = CHOICES$corr.bdiv.layers),
      weighted     = ARGS$weighted,
      tree         = ARGS$tree,
      within       = ARGS$within,
      between      = ARGS$between,
      stat.by      = like(CORR$stat.by,  ops = c('==', '!=')),
      facet.by     = like(CORR$facet.by, ops = c('==', '!=')),
      colors       = ARGS$colors,
      shapes       = CORR$shapes,
      test         = arg('choice', choices = CHOICES$corr.bdiv.test, max_n = 1),
      fit          = ARGS$fit,
      at           = ARGS$at,
      level        = ARGS$level,
      p.adj        = ARGS$p.adj,
      transform    = ARGS$transform,
      ties         = ARGS$ties,
      seed         = ARGS$seed,
      permutations = like(ARGS$permutations, when = list(when_is('test', 'mantel'))),
      alt          = ARGS$alt,
      mu           = ARGS$mu,
      strata       = ARGS$strata,
      caption      = ARGS$caption,
      check        = ARGS$check,
      alpha        = ARGS$alpha,
      cpus         = ARGS$cpus,
      '...'        = ARGS$dots )),

  bdiv_heatmap = list(
    returns = RETURNS$plot,
    params  = list(
      biom         = ARGS$biom,
      bdiv         = ARGS$bdiv,
      tree         = ARGS$tree,
      tracks       = HEAT$tracks,
      grid         = HEAT$grid,
      label        = HEAT$label,
      label_size   = HEAT$label_size,
      rescale      = HEAT$rescale,
      clust        = HEAT$clust,
      trees        = HEAT$trees,
      asp          = HEAT$asp,
      tree_height  = HEAT$tree_height,
      track_height = HEAT$track_height,
      legend       = HEAT$legend,
      title        = arg('logical', forms = c('true', 'string'), max_n = 1, null = "No title.", check = 'check_title', note = "TRUE: the metric's name. FALSE is not accepted; use NULL."),
      xlab.angle   = ARGS$xlab.angle,
      alpha        = ARGS$alpha,
      cpus         = ARGS$cpus,
      '...'        = ARGS$dots )),

  bdiv_ord_plot = list(
    returns = RETURNS$plot_stats,
    params  = list(
      biom         = ARGS$biom,
      bdiv         = ARGS$bdiv,
      ord          = arg('choice', choices = CHOICES$ord, all = '.all', note = "\"tSNE\" needs the Rtsne package and \"UMAP\" the uwot package."),
      layers       = arg('layers', choices = CHOICES$ord.layers),
      stat.by      = arg('field', class = 'categorical', max_n = 1, null = "No groups."),
      facet.by     = ARGS$facet.by,
      colors       = ARGS$colors,
      shapes       = like(ARGS$shapes, when = list(when_has_any('layers', 'point'))),
      tree         = ARGS$tree,
      test         = arg('choice', choices = CHOICES$ord.test, max_n = 1),
      seed         = ARGS$seed,
      permutations = like(ARGS$permutations, when = list(when_is_not('test', 'none'), when_set('stat.by'))),
      rank         = like(ARGS$rank, null = "No taxa biplot."),
      taxa         = like(ARGS$taxa, when = list(when_set('rank'))),
      p.top        = ARGS$p.top,
      p.adj        = ARGS$p.adj,
      unc          = ARGS$unc,
      strata       = ARGS$strata,
      caption      = ARGS$caption,
      alpha        = ARGS$alpha,
      cpus         = ARGS$cpus,
      '...'        = ARGS$dots )),

  bdiv_clusters = list(
    returns = RETURNS$factor,
    params  = list(
      biom       = ARGS$biom,
      bdiv       = ARGS$bdiv1,
      weighted   = ARGS$weighted,
      normalized = ARGS$normalized,
      tree       = ARGS$tree,
      k          = arg('integer', min = 1, max_n = 1),
      seed       = ARGS$seed,
      alpha      = ARGS$alpha,
      cpus       = ARGS$cpus,
      '...'      = ARGS$dots )),

  bdiv_distmat = list(
    returns = RETURNS$dist,
    params  = list(
      biom        = ARGS$biom,
      bdiv        = ARGS$bdiv1,
      weighted    = ARGS$weighted,
      normalized  = ARGS$normalized,
      tree        = ARGS$tree,
      within      = ARGS$within,
      between     = ARGS$between,
      norm        = arg('choice', choices = CHOICES$norm, max_n = 1, matcher = 'match.arg'),
      pseudocount = arg('number', max_n = 1, null = "None (a warning is given when one is needed)."),
      power       = arg('number', max_n = 1, when = list(when_has_any('bdiv', 'minkowski'))),
      alpha       = ARGS$alpha,
      transform   = ARGS$transform,
      ties        = ARGS$ties,
      seed        = ARGS$seed,
      cpus        = ARGS$cpus )),


  bdiv_stats = list(
    returns  = RETURNS$table,
    rules    = list(list(any_of = c('stat.by', 'regr'))),
    check    = 'check_regr_or_stat_by',
    drop_ops = TRUE,
    params   = list(
      biom         = ARGS$biom,
      regr         = arg('field', class = 'numeric', max_n = 1, ops = c('==', '!='), null = "No regression.", note = "A `==` or `!=` prefix is accepted and ignored."),
      stat.by      = arg('field', class = 'categorical', max_n = 1, ops = c('==', '!='), null = "No groups.", note = "A `==` or `!=` prefix is accepted and ignored. With `regr`, the groups split the samples instead."),
      bdiv         = ARGS$bdiv,
      weighted     = ARGS$weighted,
      tree         = ARGS$tree,
      within       = like(ARGS$within,  deprecated = "4.0.0", null = "Not set.", note = "Has no effect: the tests run on samples, not sample pairs. Use `split.by` to test within each group of a field."),
      between      = like(ARGS$between, deprecated = "4.0.0", null = "Not set.", note = "Has no effect: the tests run on samples, not sample pairs."),
      split.by     = like(STATS$split.by, ops = c('==', '!='), note = "A `==` or `!=` prefix is accepted and ignored."),
      strata       = ARGS$strata,
      test         = arg(
        type         = 'choice',
        choices      = c(CHOICES$bdiv.test, setdiff(CHOICES$bdiv.regr.test, CHOICES$bdiv.test)),
        choices_when = list(
          list(when = list(when_null('regr')), choices = CHOICES$bdiv.test),
          list(when = list(when_set('regr')),  choices = CHOICES$bdiv.regr.test) ),
        max_n        = 1,
        default_note = "\"adonis2\" without `regr`, \"mantel\" with it." ),
      seed         = ARGS$seed,
      permutations = ARGS$permutations,
      p.adj        = ARGS$p.adj,
      alpha        = ARGS$alpha,
      cpus         = ARGS$cpus )),

  distmat_stats = list(
    returns = RETURNS$table,
    params  = list(
      dm           = arg('dist', r_only = TRUE),
      groups       = arg('vector', r_only = TRUE),
      test         = arg('choice', choices = CHOICES$distmat.test, max_n = 1, note = "Categorical `groups`: \"adonis2\", \"mrpp\", or \"betadisper\". Numeric `groups`: \"adonis2\" or \"mantel\". \"none\" skips the test."),
      seed         = like(ARGS$seed,         when = list(when_is_not('test', 'none'))),
      permutations = like(ARGS$permutations, when = list(when_is_not('test', 'none'))),
      strata       = arg('vector', r_only = TRUE, null = "No strata.") )),


  #________________________________________________________
  # Rarefaction.
  #________________________________________________________

  rare_corrplot = list(
    returns = RETURNS$plot_stats,
    params  = list(
      biom       = ARGS$biom,
      adiv       = ARGS$adiv,
      layers     = CORR$layers,
      rline      = like(ARGS$rline, null = "No line."),
      stat.by    = arg('field', class = 'categorical', max_n = 1, null = "No groups."),
      facet.by   = ARGS$facet.by,
      colors     = ARGS$colors,
      shapes     = CORR$shapes,
      test       = CORR$test,
      fit        = ARGS$fit,
      at         = ARGS$at,
      level      = ARGS$level,
      p.adj      = ARGS$p.adj,
      transform  = ARGS$adiv.transform,
      alt        = ARGS$alt,
      mu         = ARGS$mu,
      caption    = ARGS$caption,
      check      = ARGS$check,
      cpus       = ARGS$cpus,
      '...'      = ARGS$dots )),

  rare_multiplot = list(
    returns = RETURNS$plot_stats,
    params  = list(
      biom       = ARGS$biom,
      adiv       = ARGS$adiv,
      layers     = CORR$layers,
      rline      = like(ARGS$rline, null = "No line."),
      stat.by    = arg('field', class = 'categorical', max_n = 1, null = "No groups."),
      facet.by   = ARGS$facet.by,
      colors     = ARGS$colors,
      shapes     = CORR$shapes,
      test       = CORR$test,
      fit        = ARGS$fit,
      at         = ARGS$at,
      level      = ARGS$level,
      p.adj      = ARGS$p.adj,
      transform  = ARGS$adiv.transform,
      alt        = ARGS$alt,
      mu         = ARGS$mu,
      caption    = ARGS$caption,
      check      = ARGS$check,
      cpus       = ARGS$cpus,
      '...'      = ARGS$dots )),

  rare_stacked = list(
    returns = RETURNS$plot,
    params  = list(
      biom        = ARGS$biom,
      rline       = ARGS$rline,
      counts      = arg('logical', max_n = 1),
      labels      = arg('logical', max_n = 1),
      y.transform = arg('choice', choices = CHOICES$rare.y.transform, max_n = 1),
      '...'       = ARGS$dots )),

  rarefy = list(
    returns = RETURNS$biom,
    params  = list(
      biom    = ARGS$biom,
      depth   = arg('integer', min = 1, max_n = 1, null = "Chosen automatically, as by suggest_rarefy_depth()."),
      seed    = ARGS$seed,
      inflate = arg('logical', max_n = 1),
      clone   = ARGS$clone,
      cpus    = ARGS$cpus )),


  #________________________________________________________
  # Taxa abundance.
  #________________________________________________________

  taxa_boxplot = list(
    returns = RETURNS$plot_stats,
    params  = list(
      biom        = ARGS$biom,
      x           = like(BOX$x, null = "Taxa along the x-axis."),
      rank        = ARGS$rank,
      layers      = BOX$layers,
      taxa        = ARGS$taxa,
      unc         = ARGS$unc,
      other       = ARGS$other,
      p.top       = like(ARGS$p.top, min = NULL),
      stat.by     = BOX$stat.by,
      facet.by    = BOX$facet.by,
      colors      = ARGS$colors,
      shapes      = BOX$shapes,
      patterns    = BOX$patterns,
      test        = BOX$test,
      flip        = BOX$flip,
      stripe      = BOX$stripe,
      ci          = BOX$ci,
      level       = ARGS$level,
      p.adj       = ARGS$p.adj,
      outliers    = BOX$outliers,
      xlab.angle  = ARGS$xlab.angle,
      p.label     = ARGS$p.label,
      transform   = ARGS$transform,
      y.transform = arg('choice', choices = CHOICES$box.y.transform, max_n = 1),
      strata      = BOX$strata,
      caption     = ARGS$caption,
      '...'       = ARGS$dots )),

  taxa_corrplot = list(
    returns = RETURNS$plot_stats,
    params  = list(
      biom       = ARGS$biom,
      x          = CORR$x,
      rank       = ARGS$rank,
      layers     = CORR$layers,
      taxa       = ARGS$taxa,
      lineage    = ARGS$lineage,
      unc        = ARGS$unc,
      other      = ARGS$other,
      p.top      = like(ARGS$p.top, min = NULL),
      stat.by    = CORR$stat.by,
      facet.by   = CORR$facet.by,
      colors     = ARGS$colors,
      shapes     = CORR$shapes,
      test       = CORR$test,
      fit        = ARGS$fit,
      at         = ARGS$at,
      level      = ARGS$level,
      p.adj      = ARGS$p.adj,
      transform  = ARGS$transform,
      ties       = ARGS$ties,
      seed       = ARGS$seed,
      alt        = ARGS$alt,
      mu         = ARGS$mu,
      strata     = CORR$strata,
      caption    = ARGS$caption,
      check      = ARGS$check,
      '...'      = ARGS$dots )),

  taxa_heatmap = list(
    returns = RETURNS$plot,
    params  = list(
      biom         = ARGS$biom,
      rank         = ARGS$rank,
      taxa         = ARGS$taxa,
      tracks       = HEAT$tracks,
      grid         = HEAT$grid,
      other        = ARGS$other,
      unc          = ARGS$unc,
      lineage      = ARGS$lineage,
      label        = HEAT$label,
      label_size   = HEAT$label_size,
      rescale      = HEAT$rescale,
      trees        = HEAT$trees,
      clust        = HEAT$clust,
      dist         = HEAT$dist,
      asp          = HEAT$asp,
      tree_height  = HEAT$tree_height,
      track_height = HEAT$track_height,
      legend       = HEAT$legend,
      title        = arg('logical', forms = c('logical', 'string'), max_n = 1, null = "No title.", check = 'check_title', note = "TRUE: the rank's name. FALSE: no title."),
      xlab.angle   = ARGS$xlab.angle,
      '...'        = ARGS$dots )),

  taxa_stacked = list(
    returns = RETURNS$plot,
    params  = list(
      biom       = ARGS$biom,
      rank       = ARGS$rank,
      taxa       = ARGS$taxa,
      colors     = like(ARGS$colors,   by = 'taxa', forms = c('logical', 'values', 'named_values'), palettes = NULL),
      patterns   = like(ARGS$patterns, by = 'taxa'),
      label.by   = arg('field', class = 'any', max_n = 1, null = "No sample labels."),
      order.by   = arg('field', class = 'any', null = "Samples are arranged by similarity (see `dist` and `clust`).", note = "Samples that tie on `order.by` are arranged by similarity."),
      facet.by   = ARGS$facet.by,
      dist       = arg('choice', choices = CHOICES$dist,  max_n = 1),
      clust      = arg('choice', choices = CHOICES$clust, max_n = 1),
      other      = ARGS$other,
      unc        = ARGS$unc,
      lineage    = ARGS$lineage,
      xlab.angle = ARGS$xlab.angle,
      '...'      = ARGS$dots )),

  taxa_forest = list(
    returns = RETURNS$plot_stats,
    rules   = list(list(one_of = c('stat.by', 'regr'))),
    check   = 'check_stat_by_regr',
    params  = list(
      biom       = ARGS$biom,
      stat.by    = arg('field', class = 'categorical', max_n = 1, null = "Not used; give `regr` instead."),
      regr       = arg('field', class = 'numeric', max_n = 1, null = "Not used; give `stat.by` instead."),
      rank       = ARGS$rank,
      taxa       = like(ARGS$taxa, null = "All taxa."),
      method     = arg('choice', choices = CHOICES$da.method, max_n = 1),
      covariates = arg('field', class = 'any', null = "No covariates."),
      split.by   = arg('field', class = 'categorical', null = "No split."),
      strata     = ARGS$strata,
      p.top      = ARGS$p.top,
      p.sig      = arg('number', min = 0, max = 1, max_n = 1),
      colors     = arg('aesthetic', aes = 'color', forms = c('logical', 'values'), max_n = 3, note = "Up to three colors: taxa more abundant, less abundant, and not significant."),
      caption    = ARGS$caption,
      prev.cut   = arg('number', min = 0, max = 1, max_n = 1),
      level      = ARGS$level,
      p.adj      = ARGS$p.adj,
      lineage    = ARGS$lineage,
      unc        = ARGS$unc,
      '...'      = ARGS$dots )),

  taxa_da = list(
    returns = RETURNS$table,
    rules   = list(list(one_of = c('stat.by', 'regr'))),
    check   = 'check_stat_by_regr',
    params  = list(
      biom       = ARGS$biom,
      stat.by    = arg('field', class = 'categorical', max_n = 1, null = "Not used; give `regr` instead."),
      regr       = arg('field', class = 'numeric', max_n = 1, null = "Not used; give `stat.by` instead."),
      rank       = ARGS$rank,
      taxa       = like(ARGS$taxa, null = "All taxa."),
      method     = arg('choice', choices = CHOICES$da.method, max_n = 1),
      covariates = arg('field', class = 'any', null = "No covariates."),
      split.by   = arg('field', class = 'categorical', null = "No split."),
      strata     = ARGS$strata,
      prev.cut   = arg('number', min = 0, max = 1, max_n = 1),
      level      = ARGS$level,
      p.adj      = ARGS$p.adj,
      lineage    = ARGS$lineage,
      unc        = ARGS$unc )),

  taxa_clusters = list(
    returns = RETURNS$factor,
    params  = list(
      biom   = ARGS$biom,
      rank   = ARGS$rank1,
      k      = arg('integer', min = 1, max_n = 1),
      seed   = ARGS$seed,
      '...'  = ARGS$dots )),

  taxa_stats = list(
    returns = RETURNS$table,
    rules   = STATS$rules,
    check   = 'check_stats_test',
    params  = list(
      biom       = ARGS$biom,
      regr       = STATS$regr,
      stat.by    = like(STATS$stat.by,  extra = ".taxa", note = "\".taxa\" compares the taxa with one another."),
      rank       = ARGS$rank,
      taxa       = ARGS$taxa,
      lineage    = ARGS$lineage,
      unc        = ARGS$unc,
      other      = ARGS$other,
      split.by   = like(STATS$split.by, extra = ".taxa"),
      strata     = ARGS$strata,
      covariates = STATS$covariates,
      transform  = ARGS$transform,
      test       = STATS$test,
      fit        = STATS$fit,
      at         = STATS$at,
      level      = ARGS$level,
      alt        = ARGS$alt,
      mu         = ARGS$mu,
      p.adj      = ARGS$p.adj )),

  taxa_matrix = list(
    returns = RETURNS$matrix,
    params  = list(
      biom      = ARGS$biom,
      rank      = ARGS$rank1,
      taxa      = like(ARGS$taxa, null = "All taxa."),
      lineage   = ARGS$lineage,
      sparse    = arg('logical', max_n = 1),
      unc       = ARGS$unc,
      other     = ARGS$other,
      transform = ARGS$transform,
      ties      = ARGS$ties,
      seed      = ARGS$seed )),


  #________________________________________________________
  # Export.
  #________________________________________________________

  write_biom = list(
    returns = like(RETURNS$file, note = "The path written to; or, when `file` is NULL, the content as a string."),
    params  = list(
      biom      = ARGS$biom,
      file      = ARGS$file,
      format    = arg('choice', choices = CHOICES$biom.format, max_n = 1, matcher = 'match.arg'),
      overwrite = ARGS$overwrite )),

  write_xlsx = list(
    returns = RETURNS$file,
    params  = list(
      biom      = ARGS$biom,
      file      = ARGS$file,
      depth     = arg('integer', min = 0, max_n = 1, null = "Chosen automatically, as by suggest_rarefy_depth().", note = "0 skips rarefaction."),
      seed      = ARGS$seed,
      unc       = like(ARGS$unc, matcher = 'match.arg'),
      overwrite = ARGS$overwrite ))
)




#________________________________________________________
# Check a function's arguments against its schema entry.
#
# `env` holds the arguments, already evaluated (a params environment
# from eval_envir(), or the function's own frame). The values are
# copied into a scratch environment and checked there, so `env` is
# never modified: the function goes on to convert them (abbreviations,
# rank numbers, ...) exactly as before. Returns the scratch environment,
# which holds the converted values.
#
# The checks are the validators the functions themselves use, so the
# error messages are the same.
#________________________________________________________
schema_check <- function (fn, env = parent.frame()) {

  def     <- FN_SCHEMA[[fn]]
  specs   <- def$params
  scratch <- rlang::new_environment(parent = rlang::ns_env('rbiom'))

  # Copy the values to check, plus the objects the checks consult.
  keep <- names(specs)[vapply(specs, function (s) !isTRUE(s$r_only), logical(1))]
  keep <- intersect(c('biom', 'tree', keep), names(specs))

  for (i in keep)
    if (exists(i, envir = env, inherits = FALSE))
      assign(i, get(i, envir = env, inherits = FALSE), envir = scratch)


  withCallingHandlers(
    expr = {

      # `==` / `!=` prefixes move the field to `within` / `between`, or
      # are dropped.
      ops <- names(specs)[vapply(specs, function (s) !is.null(s$ops), logical(1))]
      ops <- intersect(ops, ls(scratch))
      if (isTRUE(def$drop_ops)) {
        for (i in ops)
          if (is.character(scratch[[i]])) scratch[[i]] <- sub("^(==|!=)", "", scratch[[i]])
      } else if (length(ops) > 0) {
        validate_var_cmp(ops, env = scratch)
      }

      for (i in intersect(names(specs), ls(scratch, all.names = TRUE)))
        schema_check_arg(i, specs[[i]], scratch)

      # Constraints between arguments (`rules`).
      if (!is.null(def$check))
        SCHEMA_CHECKS[[def$check]](NULL, def, scratch)
    },

    # The function itself warns about deprecated arguments.
    lifecycle_warning_deprecated = function (w) invokeRestart("muffleWarning") )

  return (invisible(scratch))
}


#________________________________________________________
# Check one argument in `env` against its spec.
#________________________________________________________
schema_check_arg <- function (name, spec, env) {

  if (isTRUE(spec$r_only) || !is.null(spec$deprecated)) return (invisible())

  x       <- get(name, envir = env, inherits = FALSE)
  null_ok <- !is.null(spec$null)

  if (is.null(x) && null_ok) return (invisible())

  if (!is.null(spec$check))
    return (SCHEMA_CHECKS[[spec$check]](name, spec, env))

  switch(
    EXPR = spec$type,

    logical = validate_bool(
      var      = name,
      env      = env,
      null_ok  = null_ok,
      multiple = !identical(spec$max_n, 1) ),

    integer = ,
    number  = validate_var_range(
      var     = name,
      range   = arg_range(spec),
      env     = env,
      n       = if (identical(spec$max_n, 1)) 1 else NA,
      int     = identical(spec$type, 'integer'),
      null_ok = null_ok ),

    string = {
      if (!is_scalar_character(x) || is.na(x))
        cli_abort("`{name}` must be a character string, not {.type {x}}.")
    },

    choice = {

      choices <- spec[['choices']]
      for (cw in spec[['choices_when']])
        if (all(vapply(cw[['when']], when_holds, logical(1), env = env))) {
          choices <- cw[['choices']]
          break
        }

      if (identical(spec$matcher, 'match.arg')) {
        assign(name, match.arg(tolower(x), choices), envir = env)
      } else {
        validate_var_choices(
          var        = name,
          choices    = choices,
          env        = env,
          max        = spec$max_n %||% Inf,
          null_ok    = null_ok,
          all_option = spec$all )
      }
    },

    layers = validate_layers(name, spec$choices, env = env),

    field = {

      # Values such as ".taxa" aren't dataset fields; check the rest.
      extra <- if (is.character(x)) intersect(x, spec$extra)
      if (length(extra) > 0) {
        if (length(x) > (spec$max_n %||% Inf))
          cli_abort("`{name}` cannot be length {length(x)}; maximum is {spec$max_n}.")
        if (length(x) == length(extra)) return (invisible())
        assign(name, setdiff(x, extra), envir = env)
      }

      validate_biom_field(
        var      = name,
        env      = env,
        null_ok  = null_ok,
        max      = spec$max_n %||% Inf,
        col_type = if (isTRUE(spec$coerce)) NULL else switch(
          EXPR        = spec$class,
          categorical = 'cat',
          numeric     = 'num',
          NULL ) )

      if (length(extra) > 0)
        assign(name, c(get(name, envir = env), extra), envir = env)
    },

    rank = validate_rank(
      var     = name,
      env     = env,
      null_ok = null_ok,
      max     = spec$max_n %||% Inf ),

    taxa = validate_taxa(name, env = env, null_ok = null_ok),

    metric = {
      multiple <- !identical(spec$max_n, 1)
      if (identical(spec$div, 'alpha')) { validate_adiv(name, env = env, multiple = multiple) }
      else                              { validate_bdiv(name, env = env, multiple = multiple) }
      validate_phylo_tree(name, spec$div, env)
    },

    # Checked when the plot is built: values against the levels of `by`,
    # palette names, colors, and the parts of a list specification.
    aesthetic = NULL,
    palette   = NULL,

    cli_abort("Internal error: unknown argument type {.val {spec$type}}.")
  )

  return (invisible())
}


#________________________________________________________
# Does a `when` condition hold for the arguments in `env`?
#________________________________________________________
when_holds <- function (cond, env) {
  # [[ ]] rather than $, which would match `is` to `is_null`.
  x <- get0(cond[['param']], envir = env, inherits = FALSE)
  if (!is.null(cond[['is']]))       return (length(x) == 1 && isTRUE(x %in% cond[['is']]))
  if (!is.null(cond[['is_not']]))   return (!(length(x) == 1 && isTRUE(x %in% cond[['is_not']])))
  if (!is.null(cond[['has_any']]))  return (any(x %in% cond[['has_any']]))
  if (!is.null(cond[['is_null']]))  return (is.null(x))
  if (!is.null(cond[['not_null']])) return (!is.null(x))
  cli_abort("Internal error: unknown `when` condition.")
}


#________________________________________________________
# Checks for arguments that accept several forms. Each takes the
# argument's name, its spec, and the environment holding it.
#________________________________________________________
SCHEMA_CHECKS <- list(

  check_other = function (name, spec, env) {
    validate_string(name, env = env, null_ok = TRUE, if_true = 'Other', if_false = NULL)
  },

  check_rline = function (name, spec, env) {
    validate_rline(name, env = env, null_ok = !is.null(spec$null))
  },

  check_title = function (name, spec, env) {
    validate_title(name, env = env, false_ok = 'logical' %in% spec$forms)
  },

  check_xlab_angle = function (name, spec, env) {
    assign(name, as.character(get(name, envir = env)), envir = env)
    validate_var_choices(name, spec$choices, env = env)
  },

  check_sides = function (name, spec, env) {
    validate_sides(name, env = env, choices = spec$choices, numeric = spec$type == 'number')
  },

  check_clust = function (name, spec, env) {
    x <- get(name, envir = env)
    if (!any(is_na(x), is_false(x), inherits(x, 'hclust')))
      validate_clust(name, env = env, max = 2, na_ok = TRUE)
  },

  check_dist = function (name, spec, env) {
    if (!inherits(get(name, envir = env), 'dist'))
      validate_dist(name, env = env, max = 2)
  },

  check_stat_by_regr = function (name, spec, env) {
    validate_stat_by_regr(env = env)
  },

  check_regr_or_stat_by = function (name, spec, env) {
    validate_regr_or_stat_by(env = env)
  },

  check_stats_test = function (name, spec, env) {
    validate_stats_test(env = env)
  },

  check_tracks = function (name, spec, env) {
    x <- get(name, envir = env)
    if (is.character(x)) { validate_biom_field(name, env = env, max = Inf) }
    else if (!is.list(x)) { cli_abort("Invalid `{name}` specification.") }
  }
)
