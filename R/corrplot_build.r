
#________________________________________________________
# Build a scatterplot/trendline figure from a params environment.
#
# Inputs read from `params` (see stats_corrplot() for their meaning):
#   df, x, y, layers, stat.by, facet.by, colors, shapes, test, fit, at,
#   level, p.adj, p.top, alt, mu, strata, caption, check, .dots, and
#   optionally rline (rare_corrplot), .labs (default labels), and for beta
#   diversity corrplots biom, bdiv, tree, alpha, cpus, seed, permutations.
#
# The trendline (`fit`) and residual tables are computed with dplyr
# pipelines whose code text is recorded alongside, so that the generated
# code reproduces them from `data`.
#________________________________________________________

corrplot_build <- function (params) {

  stopifnot(is_bare_environment(params))

  dv <- params$.data_var %||% "data"
  params$.dv <- dv
  params$.data_code  %<>% if.null(character())
  params$.dots_used  %<>% if.null(character())
  default('p.top', Inf, env = params)
  default('strata', NULL, env = params)


  with(params, {


    #________________________________________________________
    # Validate and pre-process user's arguments.
    #________________________________________________________
    if (!inherits(df, 'data.frame'))
      cli_abort("`df` must be a data.frame, not {.type {df}}.")

    # Enforce unnamed vectors.
    for (.i in colnames(df)) df[[.i]] %<>% unname()
    remove(".i")

    .df0 <- df

    if (any(c(stat.by, facet.by, strata) %in% c(x, y)))
      cli_abort("`stat.by`, `facet.by`, and `strata` must differ from `x` and `y`.")

    validate_df_field('x',        col_type = "num")
    validate_df_field('y',        col_type = "num")
    validate_df_field('stat.by',  col_type = "cat", null_ok = TRUE)
    validate_df_field('facet.by', col_type = "cat", null_ok = TRUE, max = Inf)

    # `strata` is a sample field: it's in `df` except for beta diversity
    # corrplots, where `df` holds sample pairs and the tests use `biom`.
    if (eq(y, '.distance') && exists('biom', inherits = FALSE)) {
      validate_biom_field('strata', col_type = "cat", null_ok = TRUE)
    } else {
      validate_df_field('strata', col_type = "cat", null_ok = TRUE)
    }

    validate_var_choices('test',  choices = c('none', 'emmeans', 'emtrends', 'mantel'))
    validate_var_choices('fit',   choices = c('lm', 'log', 'gam'))
    validate_var_choices('p.adj', choices = p.adjust.methods)
    validate_var_choices('alt',   choices = c("!=", ">", "<"))

    validate_var_range('at',    null_ok = TRUE)
    validate_var_range('level', n = 1, range = c(0.5, 1))
    validate_var_range('mu',    n = 1)

    validate_bool('caption')
    validate_bool('check')


    #________________________________________________________
    # Consistent naming for plot_* functions.
    #________________________________________________________
    .ggdata <- as_tibble(df)
    .xcol   <- x
    .ycol   <- y
    .xmode  <- "numeric"


    #________________________________________________________
    # Drop rows without finite x and y values. Record the cleaning (NA
    # removal and factor conversion by the validators above, plus
    # this) in the code.
    #________________________________________________________
    .ggdata <- .ggdata[is.finite(.ggdata[[.xcol]]) & is.finite(.ggdata[[.ycol]]),]

    .data_code <- c(.data_code, plot_data_prep_code(
      before   = .df0,
      after    = .ggdata,
      cat_cols = unique(c(stat.by, facet.by, if (!eq(.ycol, '.distance')) strata)),
      num_cols = c(.xcol, .ycol),
      dv       = .dv ))

    remove("df", "x", "y", ".df0")


    #________________________________________________________
    # Set lower/upper bounds for confidence interval range.
    #________________________________________________________
    if (is.null(.dots$ci.min)) .dots$ci.min <- -Inf
    if (is.null(.dots$ci.max)) .dots$ci.max <-  Inf
    .dots_used <- c(.dots_used, 'ci.min', 'ci.max')

  })



  #________________________________________________________
  # Convert user's `layers` spec to layer names.
  #________________________________________________________
  init_layers(
    params  = params,
    choices = c(
      't' = "trend", 'c' = "confidence", 'p' = "point",
      'n' = "name", 'r' = "residual" ) )


  #________________________________________________________
  # Default labels: the response label from the table function, plus
  # any set by the calling function. The user's `labs.*` win.
  #________________________________________________________
  if (eq(params$.ycol, attr(params$.ggdata, 'response', exact = TRUE)))
    set_layer(params, 'labs', y = attr(params$.ggdata, 'resp_label', exact = TRUE))

  if (length(params$.labs) > 0)
    set_layer(params, 'labs', params$.labs)



  #________________________________________________________
  # All layers will have a common `x` and `y` field name.
  #________________________________________________________
  set_layer(
    params = params,
    layer  = 'ggplot',
    'mapping|x' = params$.xcol,
    'mapping|y' = params$.ycol )



  #________________________________________________________
  # Display sample names at x,y coordinates.
  #________________________________________________________
  if (has_layer(params, 'name')) {

    if (!hasName(params$.ggdata, '.sample'))
      cli_abort("The 'name' layer needs a `.sample` column in `df`.")

    set_layer(params, 'name', 'mapping|label' = ".sample")
  }



  #________________________________________________________
  # Rarefaction curve's vertical line.
  #________________________________________________________
  if (!is_null(params$rline))
    set_layer(
      params     = params,
      layer      = 'vline',
      xintercept = params$rline,
      color      = "red",
      linetype   = "dashed" )



  #________________________________________________________
  # Theming specific to a numeric x-axis.
  #________________________________________________________
  .element_blank <- P('ggplot2::element_blank')
  set_layer(params, 'theme', panel.grid.minor.x = .element_blank())




  #________________________________________________________
  # Build the plot.
  #________________________________________________________
  # The statistics run first so that `p.top` has dropped taxa before
  # the trendlines are fitted.
  p <- params %>%
    corrplot_stats() %>%
    corrplot_models() %>%
    plot_facets() %>%
    plot_build()


  return (p)
}




#________________________________________________________
# Fit one model per stat.by / facet.by group and derive from it:
#
#   attr(.ggdata, 'fit')      - the trendline (and confidence band) at
#                               100 x positions ('trend', 'confidence').
#   attr(.ggdata, 'residual') - fitted values and residuals at the
#                               observed x positions ('residual', check).
#
# Each is computed with a dplyr pipeline, and the same pipeline is
# recorded as code in params$.derived_code. Trendlines describe the
# data, so they are fitted without the `strata` term.
#________________________________________________________
corrplot_models <- function (params) {

  want_fit <- any(has_layer(params, c('confidence', 'trend')))
  want_res <- has_layer(params, 'residual') || isTRUE(params$check)

  if (!want_fit && !want_res)
    return (invisible(params))


  ggdata   <- params$.ggdata
  xcol     <- params$.xcol
  ycol     <- params$.ycol
  fit      <- params$fit
  level    <- params$level
  ci.min   <- params$.dots$ci.min
  ci.max   <- params$.dots$ci.max
  dv       <- params$.data_var %||% "data"
  by       <- unique(c(params$stat.by, params$facet.by))
  band     <- has_layer(params, 'confidence')
  n_points <- if (fit == 'lm' && !band) 2 else 100


  #________________________________________________________
  # Clamping of the fitted values to [ci.min, ci.max].
  #________________________________________________________
  clamp <- function (v) pmin(pmax(v, ci.min), ci.max)

  clamp_code <- function (expr) {
    if (is.finite(ci.min)) expr <- sprintf("pmax(%s, %s)", expr, deparse_vec(ci.min))
    if (is.finite(ci.max)) expr <- sprintf("pmin(%s, %s)", expr, deparse_vec(ci.max))
    expr
  }


  #________________________________________________________
  # Shared prologue: the model for one group.
  #________________________________________________________
  model_fun <- function (d) {
    m <- data.frame(.regr = d[[xcol]], .resp = d[[ycol]])
    list(m = m, model = suppressWarnings(stats_fit_model(m, fit = fit)))
  }

  model_code <- c(
    sprintf("m <- data.frame(.regr = d[[%s]], .resp = d[[%s]])", double_quote(xcol), double_quote(ycol)),
    strsplit(stats_model_code("m", fit), "\n", fixed = TRUE)[[1]] )


  #________________________________________________________
  # Run a per-group pipeline, collecting the groups whose model failed.
  #________________________________________________________
  failed <- character()

  pipeline <- function (fun) {
    res <- ggdata
    if (length(by) > 0) res <- dplyr::group_by(res, !!!rlang::syms(by))
    res <- res %>%
      dplyr::group_modify(function (d, keys) {
        tryCatch(
          error = function (e) {
            failed <<- c(failed, paste(as.character(unlist(keys)), collapse = " / "))
            tibble() },
          expr  = fun(d) )
      }) %>%
      dplyr::ungroup()
    as_tibble(res)
  }

  pipeline_code <- function (var, body) {
    paste(collapse = "\n", c(
      sprintf("%s <- %s %%>%%", var, dv),
      if (length(by) > 0) sprintf("  dplyr::group_by(%s) %%>%%", paste(coan(by), collapse = ", ")),
      "  dplyr::group_modify(function (d, keys) {",
      "    tryCatch(error = function (e) tibble::tibble(), expr = {",
      paste0("      ", c(model_code, body)),
      "    })",
      "  }) %>%",
      "  dplyr::ungroup()" ))
  }


  #________________________________________________________
  # Trendline: means (and optionally CI) at `n_points` x positions.
  #________________________________________________________
  if (want_fit) {

    fit_tbl <- pipeline(function (d) {

      mm <- model_fun(d)
      at <- list(.regr = seq(from = min(mm$m$.regr), to = max(mm$m$.regr), length.out = n_points))

      if (band) {
        em <- summary(emmeans::emmeans(mm$model, specs = ".regr", at = at, infer = TRUE, level = level))
        tibble(
          !!xcol := em$.regr,
          !!ycol := clamp(em$emmean),
          .ymin   = pmax(em$lower.CL, ci.min),
          .ymax   = pmin(em$upper.CL, ci.max) )
      } else {
        em <- summary(emmeans::emmeans(mm$model, specs = ".regr", at = at))
        tibble(
          !!xcol := em$.regr,
          !!ycol := clamp(em$emmean) )
      }
    })

    attr(ggdata, 'fit') <- fit_tbl

    params$.derived_code[['fit']] <- paste(collapse = "\n", c(
      "# Trendline: estimated marginal means from the model fitted to each group.",
      pipeline_code("fit", c(
        sprintf("at <- list(.regr = seq(from = min(m$.regr), to = max(m$.regr), length.out = %d))", n_points),
        if (band) c(
          sprintf("em <- summary(emmeans::emmeans(model, specs = \".regr\", at = at, infer = TRUE, level = %s))", deparse_vec(level)),
          "tibble::tibble(",
          sprintf("  %s = em$.regr,", coan(xcol)),
          sprintf("  %s = %s,", coan(ycol), clamp_code("em$emmean")),
          sprintf("  .ymin = %s,", if (is.finite(ci.min)) sprintf("pmax(em$lower.CL, %s)", deparse_vec(ci.min)) else "em$lower.CL"),
          sprintf("  .ymax = %s )", if (is.finite(ci.max)) sprintf("pmin(em$upper.CL, %s)", deparse_vec(ci.max)) else "em$upper.CL") )
        else c(
          "em <- summary(emmeans::emmeans(model, specs = \".regr\", at = at))",
          sprintf("tibble::tibble(%s = em$.regr, %s = %s)", coan(xcol), coan(ycol), clamp_code("em$emmean")) ) )) ))

    if (band)
      set_layer(
        params = params,
        layer  = 'confidence',
        'mapping|ymin' = ".ymin",
        'mapping|ymax' = ".ymax",
        'alpha'        = 0.25,
        'linewidth'    = 0.2 )

    if (is_empty_df(fit_tbl))
      del_layer(params, c('confidence', 'trend'))
  }


  #________________________________________________________
  # Residuals: fitted values at the observed x positions.
  #________________________________________________________
  if (want_res) {

    res_tbl <- pipeline(function (d) {

      mm <- model_fun(d)
      m  <- mm$m
      pr <- stats::predict(mm$model, newdata = m, se.fit = TRUE)

      tibble(
        !!xcol    := m$.regr,
        !!ycol    := m$.resp,
        .fitted    = as.vector(pr$fit),
        .residual  = m$.resp - as.vector(pr$fit),
        .std.res   = (m$.resp - as.vector(pr$fit)) / (sqrt(nrow(m)) * as.vector(pr$se.fit)),
        .yend      = clamp(as.vector(pr$fit)),
        .se        = as.vector(pr$se.fit),
        .sd        = sqrt(nrow(m)) * as.vector(pr$se.fit) )
    })

    attr(ggdata, 'residual') <- res_tbl

    params$.derived_code[['residual']] <- paste(collapse = "\n", c(
      "# Residuals: distance from each point to the trendline fitted to its group.",
      pipeline_code("residuals", c(
        "pr <- stats::predict(model, newdata = m, se.fit = TRUE)",
        "tibble::tibble(",
        sprintf("  %s = m$.regr,", coan(xcol)),
        sprintf("  %s = m$.resp,", coan(ycol)),
        "  .fitted   = as.vector(pr$fit),",
        "  .residual = m$.resp - as.vector(pr$fit),",
        "  .std.res  = (m$.resp - as.vector(pr$fit)) / (sqrt(nrow(m)) * as.vector(pr$se.fit)),",
        sprintf("  .yend     = %s,", clamp_code("as.vector(pr$fit)")),
        "  .se       = as.vector(pr$se.fit),",
        "  .sd       = sqrt(nrow(m)) * as.vector(pr$se.fit) )" )) ))

    if (has_layer(params, 'residual'))
      set_layer(params, 'residual', 'mapping|yend' = ".yend")

    if (is_empty_df(res_tbl))
      del_layer(params, 'residual')
  }


  if (length(failed <- unique(failed)) > 0)
    cli_warn(c(
      '!' = "The {fit} model could not be fitted for {qty(length(failed))}{?this/these} group{?s}: {.val {failed}}.",
      'i' = "No trendline is drawn for {qty(length(failed))}{?it/them}. Try `fit = \"lm\"`, or more samples per group." ))


  params$.ggdata <- ggdata

  return (invisible(params))
}
