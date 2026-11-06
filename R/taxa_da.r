#' Differential abundance of taxa.
#'
#' Finds the taxa whose abundance differs between groups (`stat.by`) or
#' changes with a numeric variable (`regr`). Unlike [taxa_stats()], which
#' tests each taxon's abundance on its own, these methods account for the
#' compositional nature of sequencing data: a sample's counts are a
#' sample from the community, so an increase in one taxon lowers the
#' observed proportions of all the others. The method is LinDA (Zhou et
#' al. 2022), built into rbiom.
#'
#' LinDA expects unrarefied counts: it models differences in sequencing
#' depth itself, and rarefying discards data.
#'
#' @inherit documentation_default
#'
#' @family taxa_abundance
#' @family stats_tables
#'
#' @param stat.by   Dataset field with the groups to compare. Must be
#'        categorical. Each pair of groups is compared, in the order of the
#'        field's levels: with levels `A`, `B`, and `C`, the comparisons
#'        are `"B - A"`, `"C - A"`, and `"C - B"`, and a positive log fold
#'        change means the taxon is more abundant in the first-named group.
#'        Give either `stat.by` or `regr`, not both. Default: `NULL`
#'
#' @param regr   Dataset field with a numeric variable. The log fold change
#'        is then per unit increase in `regr`. Give either `stat.by` or
#'        `regr`, not both. Default: `NULL`
#'
#' @param rank   What rank(s) of taxa to test. E.g. `"Phylum"`, `"Genus"`,
#'        `".otu"`, etc. An integer vector can also be given, where `1` is
#'        the highest rank, `-1` is the lowest rank, and `0` is the OTU
#'        "rank". Each rank is analyzed separately. Default: `-1`
#'
#' @param taxa   Which taxa to report. An integer value reports the top n
#'        most abundant taxa, a value 0 <= n < 1 reports the taxa with that
#'        mean relative abundance or greater, and a character vector reports
#'        the named taxa. Every taxon that passes `prev.cut` is still part
#'        of the analysis, since LinDA estimates the compositional bias
#'        from all of them, and p-values are adjusted for all of them.
#'        Default: `NULL` (all tested taxa).
#'
#' @param method   The differential abundance method. Currently only
#'        `"linda"`. See Details. Default: `"linda"`
#'
#' @param covariates   Dataset field(s) to adjust for, such as age or
#'        antibiotic use. They are included in every model, but their own
#'        effects are not reported. Default: `NULL`
#'
#' @param split.by   Dataset field(s) that the samples should be split by.
#'        Each subset is analyzed separately. Must be categorical.
#'        Default: `NULL`
#'
#' @param strata   Dataset field with a blocking variable, such as a subject
#'        identifier for repeated measures. Must be categorical. Adds a
#'        random intercept per stratum to every model, fitted as linear
#'        mixed models with `lmerTest::lmer()` (requires the lmerTest
#'        package). Default: `NULL`
#'
#' @param prev.cut   Prevalence filter: taxa found in fewer than this
#'        proportion of the samples are not tested. Default: `0.1`
#'
#' @param level   The confidence level for the confidence intervals
#'        (`.lower` and `.upper`). Default: `0.95`
#'
#' @param p.adj   Method to use for multiple comparisons adjustment of
#'        p-values. Run `p.adjust.methods` for a list of available options.
#'        As is standard for differential abundance, p-values are adjusted
#'        across taxa separately for each comparison (and each rank and
#'        `split.by` subset). Default: `"fdr"`
#'
#' @details
#' **Counts.** The counts are summed to `rank` with [taxa_matrix()] (using
#' `unc` and `lineage`) and every sample with a value for all of the fields
#' used is included. Samples with an `NA` in `stat.by`, `regr`,
#' `covariates`, `split.by`, or `strata` are dropped (with a message).
#'
#' **LinDA** (Zhou et al. 2022) fits one linear model per taxon to the
#' centered log2-ratios of the counts. Compositionality shifts every
#' taxon's coefficient by the same unknown amount; assuming that most taxa
#' don't change, that shift is the mode of the coefficients across taxa,
#' and subtracting it gives the log2 fold changes. rbiom's implementation
#' follows the published algorithm and the defaults of
#' `MicrobiomeStat::linda()`:
#'
#' 1. Taxa present in fewer than `prev.cut` of the samples are removed,
#'    and then any samples without counts.
#' 2. Each taxon's relative abundances are winsorized at its 97th
#'    percentile, and converted back to (rounded) counts.
#' 3. Zeros: when the sequencing depth (log of the sample sums) is
#'    associated with the model's variables (any p <= 0.1), each zero is
#'    imputed as the sample's depth divided by the largest depth among
#'    the samples where the taxon is absent. Otherwise a pseudo-count of
#'    0.5 is added to every count.
#' 4. Centered log2-ratio transformation (each sample's log2 counts minus
#'    their mean).
#' 5. A linear model per taxon, `clr ~ stat.by + covariates` (or
#'    `regr`), fitted with [stats::lm()], or with `lmerTest::lmer()` and a
#'    random intercept per stratum when `strata` is given. Numeric
#'    variables are standardized for fitting; `regr`'s fold changes are
#'    reported per unit.
#' 6. Each comparison is a coefficient (`regr`, or a group versus the
#'    first group) or a difference of two coefficients (two other
#'    groups). Its bias is its mode across taxa, estimated by the
#'    mean-shift algorithm with a Gaussian kernel (the same estimator as
#'    `modeest::mlv(method = "meanshift")`), and is subtracted.
#' 7. The standard error comes from the coefficients' covariance, and the
#'    p-value from a t-test with the model's residual degrees of freedom
#'    (Satterthwaite's for mixed models).
#'
#' For a numeric `regr` and for comparisons with the first group, this
#' gives the results of `MicrobiomeStat::linda()`. For a comparison
#' between two other groups, `MicrobiomeStat::linda.wald.test()` would
#' subtract the difference of the two coefficients' modes; rbiom takes the
#' mode of the difference instead, which assumes that most taxa don't
#' differ between those two groups (rather than between each of them and
#' the first group), and doesn't depend on which group comes first.
#'
#' One deliberate difference from `MicrobiomeStat::linda()` (version 1.4):
#' when it reports that sequencing depth is associated with the model, it
#' still adds the pseudo-count, because it compares
#' `"Imputation"` with `"imputation"`; rbiom imputes, as the paper
#' describes.
#'
#' @return A tibble data.frame with one row per taxon and comparison, and
#'         the columns:
#'        \describe{
#'          \item{`.rank` - }{ The rank. Only present when several ranks
#'                were requested. }
#'          \item{`split.by` fields - }{ The subset of samples. }
#'          \item{`stat.by` field - }{ The comparison, e.g. `"Stool - Saliva"`
#'                (positive fold changes: more abundant in Stool). Absent
#'                when `regr` is given. }
#'          \item{`.taxa` - }{ The taxon. }
#'          \item{`.lfc` - }{ The bias-corrected log2 fold change: the
#'                difference in log2 abundance between the two groups, or
#'                the change per unit of `regr`. }
#'          \item{`.se` - }{ The standard error of `.lfc`. }
#'          \item{`.lower`, `.upper` - }{ The confidence interval of `.lfc`
#'                at `level`. }
#'          \item{`.t.ratio` - }{ The test statistic, `.lfc / .se`. }
#'          \item{`.df` - }{ Degrees of freedom of the t-test. }
#'          \item{`.p.val` - }{ The p-value. }
#'          \item{`.adj.p` - }{ `.p.val` adjusted for multiple comparisons
#'                across taxa (see `p.adj`). }
#'          \item{`.n` - }{ The number of samples analyzed. }
#'        }
#'         Rows are sorted by `.p.val`. The R code that reproduces the
#'         table is in `$code`.
#'
#' @references
#' Zhou H, He K, Chen J, Zhang X (2022). LinDA: linear models for
#' differential abundance analysis of microbiome compositional data.
#' *Genome Biology* 23:95. \doi{10.1186/s13059-022-02655-5}
#'
#' @seealso [taxa_forest()] to plot the results.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     biom <- subset(hmp50, `Body Site` %in% c("Saliva", "Stool"))
#'
#'     da <- taxa_da(biom, stat.by = "Body Site", rank = "Genus")
#'     da
#'
#'     # Adjusting for a covariate.
#'     taxa_da(biom, regr = "Age", covariates = "Sex", rank = "Phylum")
#'
#'     # See the R code that reproduces the table.
#'     da$code

taxa_da <- function (
    biom, stat.by = NULL, regr = NULL, rank = -1, taxa = NULL,
    method = "linda", covariates = NULL, split.by = NULL, strata = NULL,
    prev.cut = 0.1, level = 0.95, p.adj = "fdr",
    lineage = FALSE, unc = "singly" ) {

  biom   <- as_rbiom(biom)
  params <- eval_envir(environment())
  cmd    <- current_cmd('taxa_da', params)


  #________________________________________________________
  # Validate the arguments.
  #________________________________________________________
  validate_rank(max = Inf)
  validate_taxa(null_ok = TRUE)
  validate_unc()
  validate_bool('lineage')

  validate_biom_field('stat.by',    null_ok = TRUE, col_type = "cat")
  validate_biom_field('regr',       null_ok = TRUE, col_type = "num")
  validate_biom_field('covariates', null_ok = TRUE, max = Inf)
  validate_biom_field('split.by',   null_ok = TRUE, col_type = "cat", max = Inf)
  validate_biom_field('strata',     null_ok = TRUE, col_type = "cat")

  validate_var_choices('method', "linda")
  validate_var_choices('p.adj',  p.adjust.methods)
  validate_var_range('prev.cut', range = c(0, 1),   n = 1)
  validate_var_range('level',    range = c(0.5, 1), n = 1)

  if (is.null(stat.by) == is.null(regr))
    cli_abort(c(
      'x' = "Give either `stat.by` (groups to compare) or `regr` (a numeric variable), not {if (is.null(stat.by)) 'neither' else 'both'}.",
      'i' = "To adjust for other variables, use `covariates`." ))

  if (as.logical(anyDuplicated(c(stat.by, regr, covariates, split.by, strata))))
    cli_abort("`stat.by`, `regr`, `covariates`, `split.by`, and `strata` must all be different fields.")

  if (isTRUE(biom$relative))
    cli_abort(c(
      'x' = "`biom` holds relative abundances, but differential abundance methods need counts.",
      'i' = "Use the unrarefied counts." ))

  if (!is.null(strata))
    require_package('lmerTest', reason = 'for `taxa_da()` with `strata`')

  if (!is.null(biom$depth))
    cli::cli_inform(
      .frequency    = "once",
      .frequency_id = "rbiom_taxa_da_rarefied",
      message       = c('i' = "`biom` is rarefied. LinDA models sequencing depth itself, so unrarefied counts are better." ))


  #________________________________________________________
  # Drop samples with NA values in the fields used.
  #________________________________________________________
  fields <- c(stat.by, regr, covariates, split.by, strata)
  keep   <- stats::complete.cases(biom$metadata[fields])

  if (!any(keep))
    cli_abort("No samples remain after removing those with NA values in {.field {fields}}.")

  if (!all(keep)) {
    cli::cli_inform(c('i' = "Dropped {sum(!keep)} of {length(keep)} sample{?s} with NA values in {.field {fields}}."))
    biom <- na.omit(biom, fields = fields)
  }


  #________________________________________________________
  # The full table (every tested taxon) is cached without `taxa`,
  # which only selects the rows to report.
  #________________________________________________________
  cache_file <- local({
    key <- as.list(params)
    key$taxa <- NULL
    get_cache_file('taxa_da', key)
  })

  stats <- get_cache_value(cache_file)

  if (is.null(stats)) {

    stats <- da_run(
      biom       = biom,
      stat.by    = stat.by,
      regr       = regr,
      rank       = rank,
      method     = method,
      covariates = covariates,
      split.by   = split.by,
      strata     = strata,
      prev.cut   = prev.cut,
      level      = level,
      p.adj      = p.adj,
      lineage    = lineage,
      unc        = unc )

    set_cache_value(cache_file, stats)
  }


  #________________________________________________________
  # Report only the requested taxa.
  #________________________________________________________
  tbl_sum <- attr(stats, 'tbl_sum', exact = TRUE)
  method  <- attr(stats, 'method',  exact = TRUE)

  if (!is.null(taxa)) {

    keep <- unlist(lapply(rank, function (r) {
      mtx <- taxa_matrix(biom, rank = r, lineage = lineage, unc = unc, sparse = TRUE)
      paste(r, top_taxa(mtx, taxa), sep = "\r")
    }))

    r     <- if (hasName(stats, '.rank')) as.character(stats[['.rank']]) else rep(rank[[1]], nrow(stats))
    stats <- stats[paste(r, stats[['.taxa']], sep = "\r") %in% keep, , drop = FALSE]
  }


  #________________________________________________________
  # The code: a description of the method and the call that
  # reproduces this table.
  #________________________________________________________
  stats <- as_rbiom_tbl(stats)
  attr(stats, 'tbl_sum') <- tbl_sum
  attr(stats, 'code')    <- paste(collapse = "\n", c(method, paste("stats <-", cmd)))

  return (prov_attach(stats, cmd, biom = params$biom))
}




#________________________________________________________
# Run the analysis for every rank and split.by subset, and combine the
# results into one table.
#________________________________________________________
da_run <- function (
    biom, stat.by, regr, rank, method, covariates, split.by, strata,
    prev.cut, level, p.adj, lineage, unc ) {


  #________________________________________________________
  # Standardized, syntactic column names for the model formulas.
  #________________________________________________________
  voi  <- if (is.null(stat.by)) ".regr" else ".stat.by"
  covs <- if (length(covariates) > 0) paste0(".cov", seq_along(covariates))

  engine <- switch(method, linda = da_linda)


  #________________________________________________________
  # One rank of one subset of samples.
  #________________________________________________________
  run <- function (biom, rank) {

    counts <- taxa_matrix(biom, rank = rank, lineage = lineage, unc = unc, sparse = FALSE)
    counts <- strip_prov(unclass(counts))
    md     <- biom$metadata[match(colnames(counts), biom$metadata[['.sample']]), , drop = FALSE]

    df <- data.frame(row.names = colnames(counts), check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(stat.by)) df[['.stat.by']] <- droplevels(as.factor(md[[stat.by]]))
    if (!is.null(regr))    df[['.regr']]    <- as.numeric(md[[regr]])
    for (i in seq_along(covariates)) {
      x <- md[[covariates[[i]]]]
      df[[covs[[i]]]] <- if (is.numeric(x)) as.numeric(x) else droplevels(as.factor(x))
    }
    if (!is.null(strata)) df[['.strata']] <- droplevels(as.factor(md[[strata]]))

    if (!is.null(stat.by) && nlevels(df[['.stat.by']]) < 2)
      stop("`stat.by` has fewer than two groups in this subset.")

    res <- engine(
      counts   = counts,
      md       = df,
      voi      = voi,
      covs     = covs,
      strata   = !is.null(strata),
      prev.cut = prev.cut,
      level    = level,
      p.adj    = p.adj )

    if (!is.null(stat.by))
      names(res)[names(res) == '.cmp'] <- stat.by

    return (res)
  }


  #________________________________________________________
  # Every subset and rank; one that fails is reported, not fatal,
  # when there are several.
  #________________________________________________________
  multi <- !is.null(split.by) || length(rank) > 1

  ranks <- rank

  guarded <- function (biom, rank) {

    if (!multi) return (run(biom, rank))

    label <- c(
      if (length(ranks) > 1) rank,
      vapply(split.by, function (v) as.character(pull(biom, v)[[1]]), character(1)) )

    tryCatch(
      error = function (e) {
        warning(sprintf("%s: %s", paste(label, collapse = " / "), conditionMessage(e)), call. = FALSE)
        NULL },
      expr  = run(biom, rank) )
  }

  msgs  <- character()
  stats <- withCallingHandlers(
    warning = function (w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart('muffleWarning') },
    expr    = {
      if (multi) {
        bdply(
          biom   = biom,
          vars   = split.by,
          iters  = list(rank = rank),
          prefix = TRUE,
          FUN    = function (biom, rank) {
            res <- guarded(biom, rank)
            if (is.null(res)) da_empty(stat.by) else res
          })
      } else {
        run(biom, rank)
      }
    })

  if (length(msgs) > 0) {
    msgs <- unique(msgs)
    msgs <- gsub("}", "}}", gsub("{", "{{", msgs, fixed = TRUE), fixed = TRUE)
    cli_warn(c(
      '!' = "Messages from the differential abundance analysis:",
      stats::setNames(msgs, rep('*', length(msgs))) ))
  }


  #________________________________________________________
  # Tidy up: `.rank` only for several ranks; comparisons as factors in
  # the order they were made; sorted by p-value.
  #________________________________________________________
  stats <- as_tibble(stats)

  if (hasName(stats, '.rank')) {
    if (length(rank) > 1) { stats[['.rank']] <- factor(stats[['.rank']], levels = rank) }
    else                  { stats[['.rank']] <- NULL                                     }
  }

  if (!is.null(stat.by))
    stats[[stat.by]] <- factor(stats[[stat.by]], levels = unique(stats[[stat.by]]))

  stats <- as_rbiom_tbl(stats[order(stats[['.p.val']]), , drop = FALSE])


  #________________________________________________________
  # Table header and the description of the method for `$code`.
  #________________________________________________________
  model <- da_formula_str(c(stat.by, regr), covariates, strata)

  attr(stats, 'tbl_sum') <- c(
    'Method' = switch(method, linda = "LinDA (Zhou et al. 2022)"),
    'Model'  = sprintf("%s, per taxon", model),
    'p.adj'  = sprintf("%s, across taxa within each comparison", p_adj_label(p.adj) %||% "none") )

  attr(stats, 'method') <- da_method_comment(method, stat.by, regr, covariates, strata)

  return (stats)
}



#________________________________________________________
# A zero-row result, for a subset that couldn't be analyzed.
#________________________________________________________
da_empty <- function (stat.by) {
  res <- tibble::tibble(
    .taxa = character(), .lfc = numeric(), .se = numeric(),
    .lower = numeric(), .upper = numeric(), .t.ratio = numeric(),
    .df = numeric(), .p.val = numeric(), .adj.p = numeric(), .n = integer() )
  if (!is.null(stat.by)) res <- tibble::add_column(res, !!stat.by := character(), .before = 1)
  res
}



#________________________________________________________
# The model, for the table header, e.g.
# "log abundance ~ `Body Site` + Age + (1 | Subject)".
#________________________________________________________
da_formula_str <- function (voi, covariates, strata) {
  rhs <- paste(coan(c(voi, covariates)), collapse = " + ")
  if (!is.null(strata)) rhs <- sprintf("%s + (1 | %s)", rhs, coan(strata))
  paste("log abundance ~", rhs)
}



#________________________________________________________
# Comment lines for `$code` that describe what taxa_da() computes.
# Field names appear only as escaped, backticked symbols (coan()).
#________________________________________________________
da_method_comment <- function (method, stat.by, regr, covariates, strata) {

  switch(
    EXPR  = method,
    linda = c(
      "# Differential abundance with LinDA (Zhou et al. 2022), as implemented in",
      "# rbiom; see ?taxa_da for the algorithm. For each rank, the counts from",
      "# taxa_matrix() are modeled per taxon as centered log2-ratios with",
      sprintf("#   %s", da_formula_str(c(stat.by, regr), covariates, strata)),
      "# and each coefficient is corrected by its mode across taxa." ))
}




#________________________________________________________
# The comparisons to report: for a categorical variable of interest,
# every pair of levels (later level minus earlier level); for a numeric
# one, its coefficient. `coefs` are the model matrix's column names.
# Returns a list of list(label, L) where L is a contrast vector over
# `coefs`.
#________________________________________________________
da_contrasts <- function (md, voi, coefs) {

  unit <- function (nm) stats::setNames(as.numeric(coefs == nm), coefs)

  if (voi == ".regr")
    return (list(list(label = NA_character_, L = unit(".regr"))))

  lvls <- levels(md[[voi]])
  coef <- function (lvl) if (lvl == lvls[[1]]) unit("") else unit(paste0(voi, lvl))

  lapply(utils::combn(seq_along(lvls), 2, simplify = FALSE), function (ij) {
    a <- lvls[[ij[[1]]]]
    b <- lvls[[ij[[2]]]]
    list(label = paste(b, "-", a), L = coef(b) - coef(a))
  })
}



#________________________________________________________
# The mode of `x` by the mean-shift algorithm with a Gaussian kernel,
# started from the mean of the shortest half of the data (the shorth)
# with a bandwidth from stats::bw.nrd0(). Equivalent to
# modeest::mlv(x, method = "meanshift", kernel = "gaussian").
#________________________________________________________
da_mode <- function (x, iter = 1000, tolerance = sqrt(.Machine$double.eps)) {

  x <- x[is.finite(x)]
  n <- length(x)

  if (n == 0) return (NA_real_)
  if (n == 1) return (x)


  # Starting point: the shorth.
  y     <- sort(x)
  k     <- ceiling(n / 2) - 1
  diffs <- y[(k + 1):n] - y[1:(n - k)]
  i     <- which(diffs == min(diffs))
  if (length(i) > 1) i <- mean(i)
  par   <- mean(y[i:(i + k)])


  bw <- stats::bw.nrd0(x)

  for (j in seq_len(iter)) {

    w <- stats::dnorm((x - par) / bw)
    if (sum(w) == 0) break

    M  <- sum(x * w) / sum(w)
    th <- abs(M / par - 1)

    if (!is.finite(th)) { if (M == par) break else { par <- M; next } }
    if (th < tolerance) break

    par <- M
  }

  return (par)
}



#________________________________________________________
# Finish an engine's results: confidence intervals, p-value adjustment
# within each comparison, and the column order.
#________________________________________________________
da_finish <- function (res, voi, level, p.adj) {

  q <- stats::qt(1 - (1 - level) / 2, df = res$.df)

  res$.lower <- res$.lfc - q * res$.se
  res$.upper <- res$.lfc + q * res$.se

  if (nrow(res) > 0) {
    grp <- ifelse(is.na(res$.cmp), "", res$.cmp)
    res$.adj.p <- stats::ave(res$.p.val, grp, FUN = function (p) stats::p.adjust(p, p.adj))
  }

  cols <- c('.cmp', '.taxa', '.lfc', '.se', '.lower', '.upper', '.t.ratio', '.df', '.p.val', '.adj.p', '.n')
  if (voi == ".regr") cols <- setdiff(cols, '.cmp')

  res <- tibble::as_tibble(res[, intersect(cols, names(res)), drop = FALSE])
  res$.n <- as.integer(res$.n)

  return (res)
}



#________________________________________________________
# Keep the samples that have counts, and drop the factor levels that
# no longer have samples.
#________________________________________________________
da_samples <- function (counts, md) {
  counts <- counts[, colSums(counts) > 0, drop = FALSE]
  md     <- md[colnames(counts), , drop = FALSE]
  md[]   <- lapply(md, function (x) if (is.factor(x)) droplevels(x) else x)
  if (hasName(md, '.stat.by') && nlevels(md[['.stat.by']]) < 2)
    stop("`stat.by` has fewer than two groups among the samples with counts.")
  list(counts = counts, md = md)
}




#________________________________________________________
# LinDA (Zhou et al. 2022). See the Details of ?taxa_da.
#
#   counts - taxa x samples matrix of counts.
#   md     - data.frame (one row per sample, in the order of the counts'
#            columns) with the columns `voi`, `covs`, and `.strata`.
#________________________________________________________
da_linda <- function (counts, md, voi, covs, strata, prev.cut, level, p.adj) {


  #________________________________________________________
  # Prevalence filter, then drop samples without counts.
  #________________________________________________________
  s      <- da_samples(counts, md)
  counts <- s$counts[rowMeans(s$counts > 0) >= prev.cut, , drop = FALSE]

  if (nrow(counts) < 2)
    stop("Fewer than two taxa are present in at least ", prev.cut * 100, "% of the samples.")


  #________________________________________________________
  # Winsorize each taxon's relative abundances at its 97th percentile.
  #________________________________________________________
  N      <- colSums(counts)
  P      <- t(t(counts) / N)
  P      <- pmin(P, apply(P, 1, stats::quantile, 0.97))
  counts <- round(t(t(P) * N))

  s      <- da_samples(counts, s$md)
  counts <- s$counts
  md     <- s$md
  n      <- ncol(counts)


  #________________________________________________________
  # Standardize numeric variables; build the model formula.
  #________________________________________________________
  fixed  <- c(voi, covs)
  scales <- numeric()

  for (col in fixed)
    if (is.numeric(md[[col]])) {
      scales[[col]] <- stats::sd(md[[col]])
      if (!is.finite(scales[[col]]) || scales[[col]] == 0)
        stop("`", if (col == ".regr") "regr" else "covariates", "` has no variation.")
      md[[col]] <- (md[[col]] - mean(md[[col]])) / scales[[col]]
    }

  rhs <- Reduce(function (a, b) call('+', a, b), lapply(fixed, as.name))
  if (isTRUE(strata)) rhs <- call('+', rhs, quote((1 | .strata)))

  frm <- function (lhs) structure(call('~', as.name(lhs), rhs), class = 'formula', .Environment = baseenv())


  #________________________________________________________
  # Zeros: impute when sequencing depth is associated with the model's
  # variables, otherwise add a pseudo-count of 0.5.
  #________________________________________________________
  if (any(counts == 0)) {

    N  <- colSums(counts)
    dd <- data.frame(md, .log.depth = log(N), check.names = FALSE)

    depth_fit <- if (isTRUE(strata)) {
      suppressMessages(suppressWarnings(lmerTest::lmer(frm('.log.depth'), data = dd)))
    } else {
      stats::lm(frm('.log.depth'), data = dd)
    }

    depth_p <- stats::coef(summary(depth_fit))[-1, "Pr(>|t|)"]

    if (any(depth_p <= 0.1, na.rm = TRUE)) {
      Nz <- matrix(N, nrow = nrow(counts), ncol = ncol(counts), byrow = TRUE)
      Nz[counts > 0] <- 0
      mx <- apply(Nz, 1, max)
      mx[mx == 0] <- 1
      counts <- counts + Nz / mx
    } else {
      counts <- counts + 0.5
    }
  }


  #________________________________________________________
  # Centered log2-ratios: samples x taxa.
  #________________________________________________________
  clr <- log2(counts)
  clr <- t(clr) - colMeans(clr)


  #________________________________________________________
  # Per-taxon models: coefficients (coefs x taxa), each taxon's
  # covariance matrix, and degrees of freedom.
  #________________________________________________________
  if (!isTRUE(strata)) {

    X   <- stats::model.matrix(stats::delete.response(stats::terms(frm('.y'))), data = md)
    fit <- stats::lm.fit(X, clr)

    if (fit$rank < ncol(X))
      stop("The model can't be estimated: `stat.by`, `regr`, and `covariates` are confounded.")

    p      <- ncol(X)
    beta   <- fit$coefficients
    XtXinv <- chol2inv(fit$qr$qr[seq_len(p), seq_len(p), drop = FALSE])
    sigma2 <- colSums(fit$residuals ^ 2) / (n - p)
    coefs  <- colnames(X)

    vcov_i <- function (i) sigma2[[i]] * XtXinv
    df_i   <- function (i, L) n - p

  } else {

    fits <- lapply(seq_len(ncol(clr)), function (i) {
      d <- data.frame(md, .y = clr[, i], check.names = FALSE)
      suppressMessages(suppressWarnings(lmerTest::lmer(frm('.y'), data = d)))
    })

    beta  <- sapply(fits, lme4::fixef)
    coefs <- rownames(beta)

    vcov_i <- function (i) as.matrix(stats::vcov(fits[[i]]))
    df_i   <- function (i, L) lmerTest::contest1D(fits[[i]], L, ddf = "Satterthwaite")[['df']]
  }

  colnames(beta) <- rownames(counts)


  #________________________________________________________
  # Each comparison is a contrast of the coefficients (for the first
  # group, a single coefficient), corrected by its mode across taxa.
  # Taking the mode of the contrast itself, rather than combining the
  # coefficients' modes, makes a comparison between two groups
  # independent of which group is the reference level.
  #________________________________________________________
  res <- lapply(da_contrasts(md, voi, coefs), function (cmp) {

    L   <- cmp$L
    est <- as.vector(crossprod(L, beta))
    lfc <- est - da_mode(sqrt(n) * est) / sqrt(n)
    se  <- vapply(seq_len(ncol(beta)), function (i) sqrt(as.numeric(crossprod(L, vcov_i(i) %*% L))), numeric(1))
    df  <- vapply(seq_len(ncol(beta)), function (i) as.numeric(df_i(i, L)), numeric(1))
    t   <- lfc / se

    # A numeric `regr` was standardized; report per unit.
    s <- if (voi == ".regr") scales[[voi]] else 1

    data.frame(
      .cmp     = cmp$label,
      .taxa    = colnames(beta),
      .lfc     = lfc / s,
      .se      = se / s,
      .t.ratio = t,
      .df      = df,
      .p.val   = 2 * stats::pt(-abs(t), df = df),
      .n       = n,
      stringsAsFactors = FALSE )
  })

  res <- do.call(rbind, res)

  return (da_finish(res, voi, level, p.adj))
}
