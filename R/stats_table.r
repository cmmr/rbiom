
#' Run statistics on a data.frame.
#'
#' A simple interface to lower-level statistics functions:
#' [stats::wilcox.test()], [stats::kruskal.test()], [emmeans::emmeans()],
#' and [emmeans::emtrends()] on models fitted with [stats::lm()] or
#' [mgcv::gam()]. The default test (`"emmeans"`) is parametric; use
#' `"wilcox"` or `"kruskal"` for rank-based tests.
#'
#' @inherit documentation_default
#'
#' @family stats_tables
#'
#' @param strata   Dataset field with a blocking variable, such as a subject
#'        identifier for repeated measures. Must be categorical. When given,
#'        regression models (`test = "emmeans"` or `"emtrends"`) include a
#'        random intercept per stratum, fitted with [mgcv::gam()] as
#'        `s(strata, bs = "re")` (a linear mixed model for `fit = "lm"` and
#'        `"log"`). Not supported by the rank tests. Default: `NULL`
#'
#' @param covariates   Dataset field(s) to adjust for, such as age or
#'        antibiotic use. Each is added to the regression model as a linear
#'        term (numeric fields) or a factor (categorical fields), and the
#'        reported means, slopes, and differences are adjusted for them:
#'        [emmeans::emmeans()] holds numeric covariates at their mean and
#'        averages over the levels of categorical ones with equal weights.
#'        Only for `test = "emmeans"` or `"emtrends"`; the rank tests can't
#'        adjust for covariates. Rows with an `NA` in a covariate are
#'        dropped. A covariate that has the same value in every observation
#'        of a comparison (e.g. within a pair of groups) gives an `NA` row
#'        and a warning naming it. Default: `NULL`
#'
#' @param at   Position(s) along `regr` where group means or local slopes
#'        are evaluated. Only used by `test = "emmeans"` with `stat.by`, and
#'        by `test = "emtrends"` with `fit = "gam"` (the other combinations
#'        test a single, global slope). Default: `NULL`, which evaluates at
#'        one prespecified position: the median of `regr` (computed within
#'        each `split.by` stratum). A vector evaluates at each position; the
#'        returned table then has one row per position, and the p-values are
#'        adjusted across all of them.
#'
#' @param p.adj   Method to use for multiple comparisons adjustment of
#'        p-values. Run `p.adjust.methods` for a list of available
#'        options. The adjustment spans the entire returned table: every
#'        combination of `split.by` strata, `stat.by` pairs, and `at`
#'        positions. Default: `"fdr"`
#'
#' @section Tests:
#'
#' Which quantity is tested depends on `test` and on whether `regr` and
#' `stat.by` are given. Rows with an `NA` in any of the `regr`, `resp`,
#' `stat.by`, `split.by`, `strata`, or `covariates` fields are dropped
#' (with a message).
#'
#' | `regr` | `stat.by` | `test = "emmeans"` | `test = "emtrends"` |
#' | ------ | --------- | ------------------ | ------------------- |
#' | no     | no        | mean of `resp` vs `mu` | (not allowed) |
#' | no     | yes       | pairwise differences in group means (`.mean.diff`) | (not allowed) |
#' | yes    | no        | `fit = "lm"`/`"log"`: the slope's t-test (`.slope`). `fit = "gam"`: the smooth term's approximate significance (`.edf`, `.f`) | the slope at `at` (`.slope`; global for `"lm"`/`"log"`, local for `"gam"`) |
#' | yes    | yes       | pairwise differences in group means at `at` (`.mean.diff`) | pairwise differences in slopes (`.slope.diff`; global for `"lm"`/`"log"`, local at `at` for `"gam"`) |
#'
#' For `fit = "log"`, the slope is with respect to `log(regr)`.
#' `mu` and `alt` apply to the mean or slope being tested; they are ignored
#' for pairwise comparisons (always `!= 0`) and for the smooth-term test.
#' Pairwise comparisons fit a separate model to each pair of groups.
#' Without `regr`, `fit` is irrelevant: the model is `resp ~ stat.by`.
#' `covariates` are added to every model (in the formula, as `.cov1`,
#' `.cov2`, ...; the table header names them), so that each of the
#' quantities above is adjusted for them. For example,
#' `stat.by = "Treatment", covariates = "Age"` compares the groups' means
#' at the average age.
#'
#' `test = "wilcox"` runs [stats::wilcox.test()] (one-sample against `mu`,
#' or two-sample for each pair of `stat.by` groups) with `exact = NULL`, so
#' R chooses an exact or normal-approximation p-value based on the sample
#' sizes and ties. `test = "kruskal"` runs [stats::kruskal.test()] across
#' all `stat.by` groups.
#'
#' @return A tibble data.frame with the fields below. The R code used to
#' generate the statistics is in `$code`.
#'
#' | **Field**     | **Description**                                      |
#' | ------------- | ---------------------------------------------------- |
#' | .mean         | Estimated marginal mean. See [emmeans::emmeans()].   |
#' | .mean.diff    | Difference in means (first group minus second).      |
#' | .slope        | Trendline slope. See [emmeans::emtrends()].          |
#' | .slope.diff   | Difference in slopes (first group minus second).     |
#' | .loc.shift    | Hodges-Lehmann estimate of the location shift between two groups (wilcox). |
#' | .pseudomedian | Hodges-Lehmann estimate of the location (one-sample wilcox). |
#' | .h1           | Alternative hypothesis.                              |
#' | .p.val        | The probability, under the null hypothesis, of a test statistic at least as extreme as the one observed. |
#' | .adj.p        | `.p.val` after adjusting for multiple comparisons across the whole table. |
#' | .effect.size  | Cohen's d effect size. See [emmeans::eff_size()].    |
#' | .lower        | Confidence interval lower bound.                     |
#' | .upper        | Confidence interval upper bound.                     |
#' | .se           | Standard error.                                      |
#' | .n            | Number of observations used.                         |
#' | .df           | Degrees of freedom.                                  |
#' | .stat         | Wilcoxon or Kruskal-Wallis rank sum statistic.       |
#' | .t.ratio      | t statistic: estimate divided by `.se`.              |
#' | .edf, .ref.df | Effective and reference degrees of freedom of the smooth term (gam). |
#' | .f            | F statistic for the smooth term (gam).               |
#' | .r.sqr        | R-squared (lm); adjusted R-squared (gam).            |
#' | .adj.r        | Adjusted R-squared (lm).                             |
#' | .dev.expl     | Proportion of deviance explained (gam).              |
#' | .aic          | Akaike Information Criterion.                        |
#' | .bic          | Bayesian Information Criterion.                      |
#' | .loglik       | Log-likelihood of the fitted model.                  |
#'
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(hmp50)
#'
#'     df <- taxa_table(biom, rank = "Family")
#'     stats_table(df, stat.by = "Body Site")[,1:6]
#'
#'     df <- adiv_table(biom)
#'     stats_table(df, stat.by = "Sex", split.by = "Body Site")[,1:7]
#'
#'     stats_table(df, regr = "Age", fit = "lm")
#'
#'     # Repeated measures: a random intercept per subject.
#'     df <- adiv_table(babies)
#'     stats_table(df, regr = "Age", stat.by = "Delivery", strata = "Subject", fit = "lm")
#'
#'     # Adjusting for covariates.
#'     df <- adiv_table(biom)
#'     stats_table(df, stat.by = "Sex", covariates = c("Age", "Body Site"))

stats_table <- function (
    df, regr = NULL, resp = attr(df, 'response'),
    stat.by = NULL, split.by = NULL, strata = NULL, covariates = NULL,
    test = "emmeans", fit = "gam", at = NULL,
    level = 0.95, alt = "!=", mu = 0, p.adj = "fdr" ) {


  if (!inherits(df, 'data.frame'))
    cli_abort("`df` must be a data.frame, not {.type {df}}.")

  params <- eval_envir(environment())
  cmd    <- current_cmd('stats_table', params)


  #________________________________________________________
  # See if this result is already in the cache. Provenance attributes
  # on `df` are excluded from the key so that identical data from
  # different pipelines share an entry.
  #________________________________________________________
  cache_file <- local({
    key <- as.list(params)
    key$df <- strip_prov(df)
    get_cache_file('stats_table', key)
  })

  stats <- get_cache_value(cache_file)


  if (is.null(stats)) {

    with(params, {

      #________________________________________________________
      # Validate the arguments. `df` is left untouched: the tests do
      # their own NA filtering and factor conversion, exactly as the
      # recorded code does.
      #________________________________________________________
      validate_df_field('regr',     col_type = "num", null_ok = TRUE,            force = FALSE, drop_na = FALSE)
      validate_df_field('resp',     col_type = "num",                            force = FALSE, drop_na = FALSE)
      validate_df_field('stat.by',  col_type = "cat", null_ok = TRUE,            force = FALSE, drop_na = FALSE)
      validate_df_field('split.by', col_type = "cat", null_ok = TRUE, max = Inf, force = FALSE, drop_na = FALSE)
      validate_df_field('strata',   col_type = "cat", null_ok = TRUE,            force = FALSE, drop_na = FALSE)
      validate_df_field('covariates',                  null_ok = TRUE, max = Inf, force = FALSE, drop_na = FALSE)

      validate_var_choices('test',  c("wilcox", "kruskal", "emmeans", "emtrends"))
      validate_var_choices('alt',   c("!=", ">", "<"))
      validate_var_choices('fit',   c("lm", "log", "gam"))
      validate_var_choices('p.adj', p.adjust.methods)

      validate_var_range('at',    null_ok = TRUE)
      validate_var_range('level', n = 1, range = c(0.5, 1))
      validate_var_range('mu',    n = 1)


      #________________________________________________________
      # Warn about invalid combinations of parameters.
      #________________________________________________________
      if (as.logical(anyDuplicated(c(regr, resp, stat.by, split.by, strata, covariates))))
        cli_abort("`regr`, `resp`, `stat.by`, `split.by`, `strata`, and `covariates` must all be unique.")

      if (is.null(regr) && !is.null(at))
        cli_warn('`at` is ignored when `regr` = NULL.')

      if (test %in% c('wilcox', 'kruskal') && !is.null(covariates))
        cli_abort(c(
          'x' = "`covariates` is not supported by rank tests (`test = \"{test}\"`).",
          'i' = "Use `test = \"emmeans\"` to adjust for covariates in a regression model." ))

      if (test %in% c('wilcox', 'kruskal')) {
        if (!eq(fit, "gam")) cli_warn('`fit` is ignored when `test` = "{test}".')
        if (!is.null(regr))  cli_warn('`regr` is ignored when `test` = "{test}".')
      }

      if (!is.null(stat.by)) {
        if (!eq(alt, "!=")) cli_warn('`alt` must be "!=" when using `stat.by`.')
        if (!eq(mu, 0))     cli_warn('`mu` must be 0 when using `stat.by`.')
      }

    })


    #________________________________________________________
    # Dispatch to relevant stats_* function.
    #________________________________________________________
    stats <- switch(
      EXPR = params$test,
      'wilcox'   = do.call(stats_wilcox,   fun_params(stats_wilcox,   params)),
      'kruskal'  = do.call(stats_kruskal,  fun_params(stats_kruskal,  params)),
      'emmeans'  = do.call(stats_emmeans,  fun_params(stats_emmeans,  params)),
      'emtrends' = do.call(stats_emtrends, fun_params(stats_emtrends, params)) )

    set_cache_value(cache_file, stats)
  }


  return (prov_attach(stats, cmd, data = df))
}





#' Test alpha diversity for associations with metadata.
#'
#' A convenience wrapper for [adiv_table()] + [stats_table()].
#'
#' @inheritParams stats_table
#' @inherit documentation_default
#' @inherit stats_table return
#' @inheritSection stats_table Tests
#'
#' @family alpha_diversity
#' @family stats_tables
#'
#' @param transform   Transformation to apply to the diversity values.
#'        Options are: `c("none", "rank", "log", "log1p", "sqrt")`. `"rank"`
#'        is useful for correcting for non-normally distributions before
#'        applying regression statistics. Default: `"none"`
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(hmp50)
#'
#'     adiv_stats(biom, stat.by = "Sex")[,1:6]
#'
#'     adiv_stats(biom, stat.by = "Sex", split.by = "Body Site")[,1:6]
#'
#'     adiv_stats(biom, stat.by = "Body Site", test = "kruskal")
#'
#'     adiv_stats(biom, regr = "Age", fit = "lm")

adiv_stats <- function (
    biom, regr = NULL, stat.by = NULL, adiv = "Shannon",
    split.by = NULL, strata = NULL, covariates = NULL, transform = "none",
    test = "emmeans", fit = "gam", at = NULL,
    level = 0.95, alt = "!=", mu = 0, p.adj = "fdr" ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  cmd    <- current_cmd('adiv_stats', params)


  #________________________________________________________
  # Compute alpha diversity values
  #________________________________________________________
  df <- adiv_table(
    biom      = biom,
    adiv      = adiv,
    md        = unique(c(regr, stat.by, split.by, strata, covariates)),
    transform = transform )

  if (nlevels(df$.adiv) > 1)
    split.by %<>% c('.adiv')


  #________________________________________________________
  # Run statistics
  #________________________________________________________
  stats <- stats_table(
    df         = df,
    regr       = regr,
    stat.by    = stat.by,
    split.by   = split.by,
    strata     = strata,
    covariates = covariates,
    test       = test,
    fit        = fit,
    at         = at,
    level      = level,
    alt        = alt,
    mu         = mu,
    p.adj      = p.adj )


  return (prov_attach(stats, cmd, data = df))
}




#' Test beta diversity for associations with metadata.
#'
#' Pairwise distances are not independent observations, so beta diversity
#' is tested with permutation tests in which samples, not distances, are
#' permuted: PERMANOVA and friends for a categorical `stat.by`, or a Mantel
#' test for a numeric `regr`. See [distmat_stats()] for the tests.
#'
#' @details
#' With `strata` (e.g. repeated samples per subject), the permutation
#' scheme depends on the variable being tested. A variable that varies
#' within strata (e.g. time point) is permuted within each stratum. A
#' variable that is constant within every stratum (e.g. delivery mode, with
#' subject as the stratum) is permuted at the stratum level: whole strata
#' exchange their values, so the effective sample size is the number of
#' strata. When the number of distinct arrangements is no larger than
#' `permutations`, all of them are enumerated and the p-value is exact. The
#' scheme is chosen separately for each row of the table (a pair of groups
#' can be a between-strata comparison even when the omnibus test is not)
#' and is reported in the table header. See [distmat_stats()] for details.
#'
#' @inheritParams distmat_stats
#' @inherit documentation_default
#'
#' @family beta_diversity
#' @family stats_tables
#'
#' @param stat.by   Dataset field with the groups to compare. Must be
#'        categorical. Any `==`/`!=` prefix is ignored. Default: `NULL`
#'
#' @param regr   Dataset field with a numeric variable to correlate the
#'        distances with (Mantel test). When both `regr` and `stat.by` are
#'        given, the Mantel test is run within each `stat.by` group.
#'        Default: `NULL`
#'
#' @param test   The permutation test. With `stat.by`: `"adonis2"`
#'        (PERMANOVA; the default), `"mrpp"`, or `"betadisper"` (homogeneity
#'        of dispersion). With `regr`: `"mantel"` (the default) or
#'        `"adonis2"` (PERMANOVA with a numeric covariate). See
#'        [distmat_stats()].
#'
#' @param strata   Dataset field with a blocking variable, such as a
#'        subject identifier for repeated measures. Must be categorical.
#'        When `stat.by` (or `regr`) varies within strata, permutations are
#'        restricted to within each stratum; when it is constant within
#'        every stratum (a between-subject variable), whole strata are
#'        permuted instead. See the Details of [distmat_stats()].
#'        `test = "mrpp"` and `"betadisper"` are not available for a
#'        between-strata variable. Default: `NULL`
#'
#' @param within,between   Accepted for compatibility with [bdiv_table()]
#'        and [bdiv_boxplot()], where they select which sample *pairs* are
#'        displayed. The tests here are computed on samples, so these
#'        arguments have no effect.
#'
#' @param p.adj   Method to use for multiple comparisons adjustment of
#'        p-values. Run `p.adjust.methods` for a list of available
#'        options. The adjustment spans the entire returned table.
#'        Default: `"fdr"`
#'
#' @return A tibble data.frame with one row per test. With a categorical
#'         `stat.by` and more than two groups, each `split.by` stratum
#'         (and `bdiv` metric) gets an omnibus test across all groups plus
#'         one pairwise test per pair of groups, computed on that pair's
#'         samples. The columns are:
#'        \describe{
#'          \item{`split.by` fields, `.bdiv` - }{ The stratum. `.bdiv` is
#'                only present when several `bdiv` metrics were requested. }
#'          \item{`stat.by` field - }{ The groups compared, joined by
#'                `" vs "`. The omnibus row lists all groups. }
#'          \item{`.groups` - }{ The number of groups in the comparison:
#'                2 for pairwise rows, more for the omnibus row. }
#'          \item{`.n`, `.stat`, `.r.sqr`, `.z`, `.p.val` - }{
#'                See [distmat_stats()]. }
#'          \item{`.adj.p` - }{ `.p.val` adjusted for multiple comparisons
#'                across the whole table. }
#'        }
#'         With a numeric `regr`, one Mantel test per stratum, with `.n`,
#'         `.stat` (Spearman's r), `.p.val`, and `.adj.p`.
#'         The R code used to compute the table is in `$code`.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(hmp50)
#'
#'     bdiv_stats(biom, stat.by = "Sex", bdiv = c("bray", "w_unifrac"))
#'
#'     biom <- subset(biom, `Body Site` %in% c('Saliva', 'Stool', 'Buccal mucosa'))
#'     bdiv_stats(biom, stat.by = "Body Site", split.by = "Sex")
#'
#'     bdiv_stats(biom, regr = "Age")
#'
#'     # Repeated measures: delivery mode is constant within each subject,
#'     # so whole subjects are permuted.
#'     b <- rarefy(subset(babies, `Age (days)` <= 30))
#'     bdiv_stats(b, stat.by = "Delivery mode", strata = "Subject ID")

bdiv_stats <- function (
    biom, regr = NULL, stat.by = NULL, bdiv = "bray",
    weighted = NULL, tree = NULL, within = NULL, between = NULL,
    split.by = NULL, strata = NULL,
    test = if (is.null(regr)) "adonis2" else "mantel",
    seed = 0, permutations = 999, p.adj = "fdr",
    alpha = 0.5, cpus = n_cpus() ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  cmd    <- current_cmd('bdiv_stats', cmd_params(params))


  #________________________________________________________
  # Validate and restructure user's arguments. '==' / '!=' prefixes
  # only select sample pairs for display, so they're dropped here.
  #________________________________________________________
  strip_cmp <- function (x) if (is.null(x)) NULL else sub("^(==|!=)", "", x)
  regr     <- strip_cmp(regr)
  stat.by  <- strip_cmp(stat.by)
  split.by <- strip_cmp(split.by)

  validate_bdiv(multiple = TRUE)
  validate_tree(null_ok = TRUE)

  validate_biom_field('regr',     null_ok = TRUE, col_type = "num")
  validate_biom_field('stat.by',  null_ok = TRUE, col_type = "cat")
  validate_biom_field('split.by', null_ok = TRUE, col_type = "cat", max = Inf)
  validate_biom_field('strata',   null_ok = TRUE, col_type = "cat")

  validate_var_choices('p.adj', p.adjust.methods)
  validate_var_range('permutations', int = TRUE, range = c(1, Inf), n = 1)
  validate_seed()
  validate_alpha()
  validate_cpus()

  if (is.null(regr) && is.null(stat.by))
    cli_abort("Either `stat.by` (groups to compare) or `regr` (a numeric covariate) is required.")

  if (is.null(regr)) {
    validate_var_choices('test', c("adonis2", "mrpp", "betadisper"))
  } else {
    validate_var_choices('test', c("mantel", "adonis2"))
    split.by <- unique(c(split.by, stat.by))
    stat.by  <- NULL
  }


  #________________________________________________________
  # Drop samples with NA values; check group sizes.
  #________________________________________________________
  fields <- unique(c(regr, stat.by, split.by, strata))
  sub    <- stats_validate(biom$metadata, regr, resp = NULL, stat.by, split.by, strata)

  na_line <- NULL
  if (nrow(sub) < biom$n_samples) {
    biom    <- na.omit(biom, fields = fields)
    na_line <- sprintf("biom <- na.omit(biom, fields = %s)", as.args(list(fields)))
  }


  #________________________________________________________
  # MRPP and betadisper can't permute whole strata.
  #________________________________________________________
  if (!is.null(strata) && test %in% c("mrpp", "betadisper") &&
      eq(strata_scheme(sub[[stat.by]], sub[[strata]]), "between"))
    cli_abort(c(
      'x' = "{.field {stat.by}} is constant within each {.field {strata}}, so observations aren't exchangeable within strata.",
      'i' = "Use `test = \"adonis2\"`, which permutes whole strata." ))


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file <- get_cache_file('bdiv_stats', params)
  stats      <- get_cache_value(cache_file)

  if (is.null(stats)) {

    #________________________________________________________
    # One stratum, one metric. The code text further down mirrors
    # this function line for line.
    #________________________________________________________
    # The permutation scheme of each successful test, for the header.
    schemes <- NULL
    note    <- function (res, groups, strata) {
      force(res)
      schemes <<- rbind(schemes, perm_info(perm_plan(groups, strata, permutations)))
      res
    }

    run <- function (biom, bdiv) {

      dm <- bdiv_distmat(biom, bdiv = bdiv, tree = tree, alpha = alpha, cpus = cpus)

      if (is.null(regr)) {

        groups  <- droplevels(as.factor(pull(biom, stat.by)))
        strata  <- if (!is.null(strata)) pull(biom, strata)
        lvls    <- levels(groups)
        sets    <- c(list(lvls), if (length(lvls) > 2) combn(lvls, 2, simplify = FALSE))

        dplyr::bind_rows(lapply(sets, function (lvls) {
          keep <- names(groups)[groups %in% lvls]
          res  <- tryCatch(
            error = function (e) { warning(conditionMessage(e), call. = FALSE); tibble::tibble(.p.val = NA_real_) },
            expr  = note(
              res = distmat_stats(
                dm           = stats::as.dist(as.matrix(dm)[keep, keep]),
                groups       = groups[keep],
                strata       = if (!is.null(strata)) strata[keep],
                test         = test,
                seed         = seed,
                permutations = permutations ),
              groups = groups[keep],
              strata = if (!is.null(strata)) strata[keep] ))
          tibble::tibble(!!stat.by := paste(lvls, collapse = ' vs '), .groups = length(lvls), res)
        }))

      } else {

        note(
          res = distmat_stats(
            dm           = dm,
            groups       = pull(biom, regr),
            strata       = if (!is.null(strata)) pull(biom, strata),
            test         = test,
            seed         = seed,
            permutations = permutations ),
          groups = pull(biom, regr),
          strata = if (!is.null(strata)) pull(biom, strata) )
      }
    }

    # With several strata or metrics, one that fails (e.g. a stratum with
    # a single group) becomes an NA row and is reported once at the end.
    multi <- !is.null(split.by) || length(bdiv) > 1

    guarded <- function (biom, bdiv) {
      tryCatch(
        error = function (e) { warning(conditionMessage(e), call. = FALSE); tibble::tibble(.p.val = NA_real_) },
        expr  = run(biom, bdiv) )
    }

    msgs  <- character()
    stats <- withCallingHandlers(
      warning = function (w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart('muffleWarning') },
      expr    = {

        stats <- if (multi) {
          bdply(
            biom   = biom,
            vars   = split.by,
            iters  = if (length(bdiv) > 1) list(bdiv = bdiv) else list(),
            prefix = TRUE,
            FUN    = if (length(bdiv) > 1) guarded else function (biom) guarded(biom, bdiv) )
        } else {
          run(biom, bdiv)
        }

        .p.val <- NULL # for CRAN check only
        stats <- dplyr::mutate(stats, .adj.p = p.adjust(.p.val, p.adj), .after = .p.val)
        stats <- dplyr::arrange(stats, .p.val)
      })

    stats_report(stats, unique(msgs), key_cols = c(split.by, '.bdiv', stat.by))


    #________________________________________________________
    # The code text.
    #________________________________________________________
    dm_args <- c(
      sprintf("bdiv = %s", if (length(bdiv) > 1) "bdiv" else double_quote(bdiv)),
      if (!is.null(tree))  "tree = tree",
      if (!eq(alpha, 0.5)) sprintf("alpha = %s", deparse_vec(alpha)) )

    ds_args <- c(
      sprintf("test         = %s,", double_quote(test)),
      sprintf("seed         = %s,", deparse_vec(seed)),
      sprintf("permutations = %s )", deparse_vec(permutations)) )

    if (is.null(regr)) {

      inner <- c(
        sprintf("dm     <- bdiv_distmat(biom, %s)", paste(dm_args, collapse = ", ")),
        sprintf("groups <- droplevels(as.factor(pull(biom, %s)))", double_quote(stat.by)),
        if (!is.null(strata)) sprintf("strata <- pull(biom, %s)", double_quote(strata)),
        "lvls   <- levels(groups)",
        "sets   <- c(list(lvls), if (length(lvls) > 2) combn(lvls, 2, simplify = FALSE))",
        "dplyr::bind_rows(lapply(sets, function (lvls) {",
        "  keep <- names(groups)[groups %in% lvls]",
        "  res  <- tryCatch(",
        "    error = function (e) { warning(conditionMessage(e), call. = FALSE); tibble::tibble(.p.val = NA_real_) },",
        "    expr  = distmat_stats(",
        "      dm           = stats::as.dist(as.matrix(dm)[keep, keep]),",
        "      groups       = groups[keep],",
        if (!is.null(strata)) "      strata       = strata[keep],",
        paste0("      ", sub(" )$", " ))", ds_args)),
        sprintf("  tibble::tibble(%s = paste(lvls, collapse = ' vs '), .groups = length(lvls), res)", coan(stat.by)),
        "}))" )

    } else {

      inner <- c(
        sprintf("dm     <- bdiv_distmat(biom, %s)", paste(dm_args, collapse = ", ")),
        "distmat_stats(",
        "  dm           = dm,",
        sprintf("  groups       = pull(biom, %s),", double_quote(regr)),
        if (!is.null(strata)) sprintf("  strata       = pull(biom, %s),", double_quote(strata)),
        paste0("  ", ds_args) )
    }


    if (multi)
      inner <- c(
        "tryCatch(",
        "  error = function (e) { warning(conditionMessage(e), call. = FALSE); tibble::tibble(.p.val = NA_real_) },",
        "  expr  = {",
        paste0("    ", inner),
        "  })" )

    body <- bdply_code(inner, split.by, bdiv)


    code <- paste(collapse = "\n", c(
      if (!is.null(tree)) sprintf("# tree <- %s   # not reproducible from code", obj_repr(tree)),
      na_line,
      if (!is.null(tree) || !is.null(na_line)) "",
      body,
      "",
      sprintf("stats <- dplyr::mutate(stats, .adj.p = p.adjust(.p.val, %s), .after = .p.val)", double_quote(p.adj)),
      "stats <- dplyr::arrange(stats, .p.val)" ))


    test_label <- switch(
      EXPR       = test,
      adonis2    = "PERMANOVA",
      mrpp       = "MRPP (vegan::mrpp)",
      betadisper = "Dispersion (vegan::betadisper)",
      mantel     = "Mantel (Spearman)" )

    scheme_text <- if (is.null(schemes)) {
      sprintf("%s permutations", permutations)
    } else {
      perm_scheme_summary(schemes, coan(strata %||% "strata"), permutations)
    }

    stats <- as_rbiom_tbl(stats)
    attr(stats, 'tbl_sum') <- c('Test' = sprintf(
      "%s ~ %s; %s.", test_label, coan(if (is.null(regr)) stat.by else regr), scheme_text ))
    attr(stats, 'code') <- code

    set_cache_value(cache_file, stats)
  }


  return (prov_attach(stats, cmd, biom = params$biom))
}





#' Test taxa abundances for associations with metadata.
#'
#' A convenience wrapper for [taxa_table()] + [stats_table()].
#'
#' @inheritParams stats_table
#' @inherit documentation_default
#' @inherit stats_table return
#' @inheritSection stats_table Tests
#'
#' @family taxa_abundance
#' @family stats_tables
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- rarefy(hmp50)
#'
#'     taxa_stats(biom, stat.by = "Body Site", rank = "Family")[,1:6]

taxa_stats <- function (
    biom, regr = NULL, stat.by = NULL, rank = -1, taxa = 6,
    lineage = FALSE, unc = "singly", other = FALSE,
    split.by = NULL, strata = NULL, covariates = NULL, transform = "none",
    test = "emmeans", fit = "gam", at = NULL,
    level = 0.95, alt = "!=", mu = 0, p.adj = "fdr" ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  cmd    <- current_cmd('taxa_stats', params)


  #________________________________________________________
  # Compute taxa abundance values
  #________________________________________________________
  df <- taxa_table(
    biom      = biom,
    rank      = rank,
    taxa      = taxa,
    lineage   = lineage,
    md        = setdiff(unique(c(regr, stat.by, split.by, strata, covariates)), ".taxa"),
    unc       = unc,
    other     = other,
    transform = transform )

  if (nlevels(df$.rank) > 1) split.by %<>% c('.rank')
  if (!eq(stat.by, ".taxa")) split.by %<>% c('.taxa')


  #________________________________________________________
  # Run statistics
  #________________________________________________________
  stats <- stats_table(
    df         = df,
    regr       = regr,
    stat.by    = stat.by,
    split.by   = split.by,
    strata     = strata,
    covariates = covariates,
    test       = test,
    fit        = fit,
    at         = at,
    level      = level,
    alt        = alt,
    mu         = mu,
    p.adj      = p.adj )


  return (prov_attach(stats, cmd, data = df))
}
