
#' Run permutation statistics on a distance matrix.
#'
#' Tests whether the distances in `dm` are structured by `groups`, using
#' permutation tests. Samples, not distances, are the unit of permutation, so
#' the n(n-1)/2 pairwise distances are never treated as independent
#' observations.
#'
#' @inherit documentation_default
#'
#' @family beta_diversity
#' @family stats_tables
#'
#' @param groups  A vector of grouping values, one per sample. Either named
#'        with the sample names in `attr(dm, 'Labels')`, or unnamed and in
#'        the same order. Categorical values define groups; numeric values
#'        are used as a covariate (`test = "adonis2"`) or correlated with
#'        the distances (`test = "mantel"`). Samples with an `NA` value are
#'        dropped (with a message). Required.
#'
#' @param test  The permutation test to run. Options are:
#'        \describe{
#'          \item{`"adonis2"` - }{ PERMANOVA; [vegan::adonis2()]. Accepts
#'                categorical or numeric `groups`. }
#'          \item{`"mrpp"` - }{ Multiple response permutation procedure;
#'                [vegan::mrpp()]. Categorical `groups` only. }
#'          \item{`"betadisper"` - }{ Homogeneity of multivariate dispersion,
#'                i.e. whether groups differ in their spread rather than
#'                their location; [vegan::betadisper()] followed by
#'                [vegan::permutest()]. Categorical `groups` only. A
#'                significant PERMANOVA together with a significant
#'                dispersion test means that group differences may be due
#'                to differences in spread. }
#'          \item{`"mantel"` - }{ Mantel test of the correlation (Spearman)
#'                between `dm` and the pairwise absolute differences in a
#'                numeric `groups`; [vegan::mantel()]. Numeric `groups` only. }
#'          \item{`"none"` - }{ Don't run any statistics; returns `NULL`. }
#'        }
#'        Abbreviations are allowed. Default: `"adonis2"`
#'
#' @param strata  A blocking variable for repeated measures, e.g. a subject
#'        identifier. Given like `groups` (one value per sample, named or in
#'        `attr(dm, 'Labels')` order). How permutations are restricted
#'        depends on whether `groups` varies within strata; see Details.
#'        Default: `NULL` (samples are exchangeable).
#'
#' @details
#' **Permutation schemes with `strata`.** Samples from the same stratum
#' (e.g. repeated samples of one subject) are not independent, so they
#' can't be shuffled freely. Which restriction is valid depends on `groups`:
#'
#' * **Within-strata variable** (`groups` varies within at least one
#'   stratum, e.g. time point or a within-subject treatment): labels are
#'   permuted only within each stratum, using
#'   `permute::how(nperm = permutations, blocks = strata)` with the
#'   \pkg{vegan} function for `test`.
#'
#' * **Between-strata variable** (`groups` is constant within every
#'   stratum, e.g. delivery mode with subject as the stratum): permuting
#'   within strata would never change the labels (every permutation would
#'   reproduce the observed statistic, giving p = 1). Instead, the
#'   stratum-level values are permuted among whole strata (strata may have
#'   different sizes), and the statistic is recomputed on all samples. The
#'   p-value is `(1 + #{permuted >= observed}) / (1 + permutations)`. When
#'   there are no more distinct arrangements of the values among strata
#'   than `permutations`, all of them are enumerated instead and the
#'   p-value is exact: the proportion of arrangements (including the
#'   observed one) with a statistic at least as large. The statistics are
#'   computed by rbiom with the same formulas as \pkg{vegan}: the
#'   PERMANOVA pseudo-F and R-squared from the Gower-centred matrix of
#'   squared distances, and Spearman's r between the distances and
#'   pairwise differences in `groups`. The effective sample size is the
#'   number of strata, so a warning is given when the null distribution is
#'   coarse (e.g. a group present in only one stratum). Only
#'   `test = "adonis2"` and `test = "mantel"` support this scheme;
#'   `"mrpp"` and `"betadisper"` give an error.
#'
#' The scheme used, and the number of permutations or enumerated
#' arrangements, is reported in the table header.
#'
#' @return A tibble data.frame with one row and the columns:
#'        \describe{
#'          \item{`.n` - }{ The number of samples in the distance matrix. }
#'          \item{`.stat` - }{ The observed test statistic: the pseudo-F
#'                ratio (adonis2, betadisper), the weighted mean within-group
#'                distance delta (mrpp), or Spearman's r (mantel). }
#'          \item{`.r.sqr` - }{ adonis2 only: the proportion of the total
#'                sum of squares explained by `groups`. }
#'          \item{`.z` - }{ mrpp only: the observed statistic minus the mean
#'                of the permuted statistics, divided by their standard
#'                deviation. }
#'          \item{`.p.val` - }{ The permutation p-value: the probability,
#'                under the null hypothesis, of a statistic at least as
#'                extreme as the observed one, estimated from the
#'                permutations (plus the observed data). Its resolution is
#'                limited by `permutations` (or by the number of distinct
#'                arrangements of whole strata; see Details). }
#'        }
#'        The R code used to calculate these statistics is in `$code`.
#'        The `groups` (and `strata`) values are embedded in that code so
#'        that it can be re-run as is.
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     hmp10 <- hmp50[1:10]
#'
#'     dm <- bdiv_distmat(hmp10, 'w_unifrac')
#'
#'     distmat_stats(dm, groups = pull(hmp10, 'Body Site'))
#'
#'     distmat_stats(dm, groups = pull(hmp10, 'Age'), test = "mantel")
#'
#'     # See the R code used to calculate these statistics:
#'     stats <- distmat_stats(dm, groups = pull(hmp10, 'Age'), test = "mantel")
#'     stats$code
#'
#'     # Repeated measures: delivery mode is constant within each subject,
#'     # so whole subjects are permuted (all 220 arrangements of the 12
#'     # infants are enumerated, giving an exact p-value).
#'     b  <- rarefy(subset(babies, `Age (days)` <= 10))
#'     dm <- bdiv_distmat(b, 'bray')
#'     distmat_stats(dm, pull(b, 'Delivery mode'), strata = pull(b, 'Subject ID'))
#'

distmat_stats <- function (dm, groups, test = "adonis2", seed = 0, permutations = 999, strata = NULL) {

  params <- eval_envir(environment())

  cmd <- local({
    p <- as.list(params)
    p$dm     <- aa(NA, display = "dm")
    p$groups <- aa(NA, display = "groups")
    if (!is.null(p$strata)) p$strata <- aa(NA, display = "strata")
    current_cmd('distmat_stats', p)
  })
  schema_check('distmat_stats', params)


  #________________________________________________________
  # Sanity checking.
  #________________________________________________________
  if (!inherits(dm, 'dist'))
    cli_abort("`dm` must be a dist object, not {.type {dm}}.")

  validate_var_choices('test', CHOICES$distmat.test)
  validate_seed()
  validate_var_range('permutations', int = TRUE, range = arg_range(ARGS$permutations), n = 1)

  if (test == "none") return (NULL)

  labels <- attr(dm, 'Labels')
  size   <- attr(dm, 'Size')

  align <- function (x, what) {
    if (is.null(names(x))) {
      if (length(x) != size)
        cli_abort("`{what}` must be named with sample names, or have one value per sample ({size}), not {length(x)}.")
      names(x) <- labels
    }
    if (length(missing <- setdiff(labels, names(x))) > 0)
      cli_abort("`{what}` is missing value{?s} for sample{?s} {.val {missing}}.")
    x[labels]
  }

  groups <- align(groups, "groups")
  strata <- if (!is.null(strata)) align(strata, "strata")

  if (is.numeric(groups) && test %in% c("mrpp", "betadisper"))
    cli_abort(c(
      'x' = "`test = \"{test}\"` requires categorical `groups`, not numeric.",
      'i' = "Use `test = \"adonis2\"` or `test = \"mantel\"` for a numeric covariate." ))

  if (!is.numeric(groups) && test == "mantel")
    cli_abort(c(
      'x' = "`test = \"mantel\"` requires numeric `groups`.",
      'i' = "Use `test = \"adonis2\"`, `\"mrpp\"`, or `\"betadisper\"` for categorical groups." ))


  #________________________________________________________
  # Drop samples with NA group / stratum values.
  #________________________________________________________
  keep <- !is.na(groups)
  if (!is.null(strata)) keep <- keep & !is.na(strata)

  if (!any(keep))
    cli_abort("All samples have NA values in `groups` or `strata`.")

  if (!all(keep))
    cli::cli_inform(c('i' = "Dropped {sum(!keep)} of {size} sample{?s} with NA values in `groups` or `strata`."))

  if (!is.numeric(groups) && length(unique(groups[keep])) < 2)
    cli_abort("`groups` must have at least two distinct values.")


  #________________________________________________________
  # Which permutation scheme? See Details.
  #________________________________________________________
  plan   <- perm_plan(groups[keep], if (!is.null(strata)) strata[keep], permutations)
  scheme <- plan$scheme

  if (scheme == "between" && test %in% c("mrpp", "betadisper"))
    cli_abort(c(
      'x' = "`groups` is constant within each stratum, so observations aren't exchangeable within strata.",
      'i' = "Use `test = \"adonis2\"`, which permutes whole strata." ))


  #________________________________________________________
  # See if this result is already in the cache. Provenance attributes
  # on `dm` are excluded from the key so that identical distances from
  # different pipelines share an entry - except whether `dm` has a `cmd`,
  # which changes the cached `code`.
  #________________________________________________________
  cache_file <- local({
    key <- as.list(params, all.names = TRUE)
    key$dm     <- strip_prov(params$dm)
    key$dm_cmd <- is.null(attr(params$dm, 'cmd', exact = TRUE))
    get_cache_file('distmat_stats', key)
  })
  stats      <- get_cache_value(cache_file)

  if (is.null(stats)) {

    #________________________________________________________
    # Align everything to the (possibly reduced) distance matrix.
    #________________________________________________________
    dm_kept <- dm
    if (!all(keep)) {
      ids     <- names(groups)[keep]
      dm_kept <- stats::as.dist(as.matrix(dm)[ids, ids])
    }

    grouping <- groups[attr(dm_kept, 'Labels')]
    blocks   <- if (!is.null(strata)) factor(strata[attr(dm_kept, 'Labels')])
    n        <- attr(dm_kept, 'Size')


    #________________________________________________________
    # Run the test.
    #________________________________________________________
    msgs <- character()

    stats <- withCallingHandlers(
      warning = function (w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart('muffleWarning') },
      expr    = if (scheme == "between") {

        with_seed(seed, between_strata_test(dm_kept, grouping, blocks, test, permutations))

      } else {

        perms <- if (is.null(blocks)) permutations else permute::how(nperm = permutations, blocks = blocks)

        with_seed(seed, switch(
          EXPR = test,

          adonis2 = {
            ptest <- vegan::adonis2(dm_kept ~ grouping, permutations = perms)
            tibble::tibble(
              .n     = n,
              .stat  = ptest$F[[1]],
              .r.sqr = ptest$R2[[1]],
              .p.val = ptest$`Pr(>F)`[[1]] )
          },

          mrpp = {
            ptest  <- vegan::mrpp(dm_kept, grouping, permutations = perms)
            pstats <- summary(vegan::permustats(ptest))
            tibble::tibble(
              .n     = n,
              .stat  = ptest$delta,
              .z     = unname(pstats$z),
              .p.val = ptest$Pvalue )
          },

          betadisper = {
            bdisp <- vegan::betadisper(dm_kept, grouping)
            ptest <- vegan::permutest(bdisp, permutations = perms)
            tibble::tibble(
              .n     = n,
              .stat  = ptest$tab[['F']][[1]],
              .p.val = ptest$tab[['Pr(>F)']][[1]] )
          },

          mantel = {
            ptest <- vegan::mantel(dm_kept, stats::dist(grouping), method = 'spearman', permutations = perms)
            tibble::tibble(
              .n     = n,
              .stat  = ptest$statistic,
              .p.val = ptest$signif )
          } ))
      })

    stats_report(stats, unique(msgs))


    #________________________________________________________
    # The code text. `groups`/`strata` are embedded as literals so that
    # it is self-contained; the dist object is `data` (from `$code`).
    #________________________________________________________
    code_head <- c(
      if (is.null(attr(dm, 'cmd'))) sprintf("# data <- %s   # not reproducible from code", obj_repr(dm)),
      "dm     <- data",
      sprintf("groups <- %s", deparse_vec(groups)),
      if (!is.null(strata)) sprintf("strata <- %s", deparse_vec(strata)) )

    if (scheme == "between") {

      code <- paste(collapse = "\n", c(
        code_head,
        "",
        "# `groups` is constant within each stratum, so rbiom permutes whole strata (see ?distmat_stats).",
        sprintf(
          "stats <- distmat_stats(dm, groups, test = %s, seed = %s, permutations = %s, strata = strata)",
          double_quote(test), deparse_vec(seed), deparse_vec(permutations) )))

    } else {

      perm_code <- if (is.null(strata)) {
        deparse_vec(permutations)
      } else {
        sprintf("permute::how(nperm = %s, blocks = blocks)", deparse_vec(permutations))
      }

      test_code <- switch(
        EXPR = test,
        adonis2 = c(
          sprintf("ptest <- vegan::adonis2(dm ~ grouping, permutations = %s)", perm_code),
          "stats <- tibble::tibble(",
          "  .n     = attr(dm, 'Size'),",
          "  .stat  = ptest$F[[1]],",
          "  .r.sqr = ptest$R2[[1]],",
          "  .p.val = ptest$`Pr(>F)`[[1]] )" ),
        mrpp = c(
          sprintf("ptest  <- vegan::mrpp(dm, grouping, permutations = %s)", perm_code),
          "pstats <- summary(vegan::permustats(ptest))",
          "stats  <- tibble::tibble(",
          "  .n     = attr(dm, 'Size'),",
          "  .stat  = ptest$delta,",
          "  .z     = unname(pstats$z),",
          "  .p.val = ptest$Pvalue )" ),
        betadisper = c(
          "bdisp <- vegan::betadisper(dm, grouping)",
          sprintf("ptest <- vegan::permutest(bdisp, permutations = %s)", perm_code),
          "stats <- tibble::tibble(",
          "  .n     = attr(dm, 'Size'),",
          "  .stat  = ptest$tab[['F']][[1]],",
          "  .p.val = ptest$tab[['Pr(>F)']][[1]] )" ),
        mantel = c(
          sprintf("ptest <- vegan::mantel(dm, stats::dist(grouping), method = 'spearman', permutations = %s)", perm_code),
          "stats <- tibble::tibble(",
          "  .n     = attr(dm, 'Size'),",
          "  .stat  = ptest$statistic,",
          "  .p.val = ptest$signif )" ))

      code <- paste(collapse = "\n", c(
        code_head,
        if (!all(keep)) c(
          if (is.null(strata)) "ids    <- names(groups)[!is.na(groups)]" else "ids    <- names(groups)[!is.na(groups) & !is.na(strata)]",
          "dm     <- stats::as.dist(as.matrix(dm)[ids, ids])" ),
        "",
        "grouping <- groups[attr(dm, 'Labels')]",
        if (!is.null(strata)) "blocks   <- factor(strata[attr(dm, 'Labels')])",
        sprintf("set.seed(%s)", deparse_vec(seed)),
        "",
        test_code ))
    }


    stats <- as_rbiom_tbl(stats)
    attr(stats, 'tbl_sum') <- c('Test' = sprintf(
      "%s ~ groups; %s.",
      perm_test_label(test, scheme),
      perm_scheme_text(scheme, plan$n, plan$exact) ))
    attr(stats, 'code') <- code

    set_cache_value(cache_file, stats)
  }


  return (prov_attach(stats, cmd, data = params$dm))
}




#________________________________________________________
# "between" when `groups` is constant within every stratum (ignoring NA),
# otherwise "within".
#________________________________________________________
strata_scheme <- function (groups, strata) {

  codes  <- match(groups, unique(groups))
  n_vals <- tapply(codes, as.character(strata), function (g) length(unique(g[!is.na(g)])))

  if (all(n_vals <= 1)) "between" else "within"
}


#________________________________________________________
# How distmat_stats() permutes for these `groups` and `strata` (NULL for
# none): list(scheme = "free"|"within"|"between", n = the number of
# permutations, or of arrangements when exact, exact = TRUE when all
# distinct arrangements of whole strata are enumerated).
#________________________________________________________
perm_plan <- function (groups, strata, permutations) {

  if (is.null(strata))
    return (list(scheme = "free", n = permutations, exact = FALSE))

  ok     <- !is.na(groups) & !is.na(strata)
  groups <- groups[ok]
  strata <- as.character(strata[ok])

  if (strata_scheme(groups, strata) == "within")
    return (list(scheme = "within", n = permutations, exact = FALSE))

  arr <- stratum_arrangements(groups[match(unique(strata), strata)], permutations)
  list(scheme = "between", n = if (arr$exact) arr$n else permutations, exact = arr$exact)
}


# The number of distinct arrangements of stratum-level values `vals`
# among the strata, and whether it is small enough to enumerate.
stratum_arrangements <- function (vals, permutations) {
  counts <- tabulate(match(vals, unique(vals)))
  log_nd <- lfactorial(length(vals)) - sum(lfactorial(counts))
  exact  <- log_nd <= log(permutations) + 1e-9
  list(
    counts = counts,
    log_n  = log_nd,
    exact  = exact,
    n      = if (exact) as.integer(round(exp(log_nd))) else NA_integer_ )
}



#________________________________________________________
# Labels for table headers.
#________________________________________________________
perm_test_label <- function (test, scheme = "free") {

  if (scheme == "between")
    return (switch(test, adonis2 = "PERMANOVA (adonis2 pseudo-F)", mantel = "Mantel (Spearman)", test))

  switch(
    EXPR       = test,
    adonis2    = "PERMANOVA (vegan::adonis2)",
    mrpp       = "MRPP (vegan::mrpp)",
    betadisper = "Dispersion (vegan::betadisper)",
    mantel     = "Mantel (vegan::mantel)",
    test )
}


# `strata_name` is a quoted field name, or NULL for a vector of strata.
perm_scheme_text <- function (scheme, n, exact, strata_name = NULL) {
  switch(
    EXPR    = scheme,
    free    = sprintf("%s permutations", n),
    within  = sprintf("%s permutations within %s", n, strata_name %||% "strata"),
    between = sprintf(
      "whole %s permuted (between-strata variable; %s)",
      if (is.null(strata_name)) "strata" else paste(strata_name, "strata"),
      if (isTRUE(exact)) sprintf("all %s arrangements, exact", n) else sprintf("%s permutations", n) ))
}


# Header text summarizing the schemes of several distmat_stats() results.
# `info` is a data.frame with columns scheme, n, exact.
perm_scheme_summary <- function (info, strata_name, permutations) {

  parts <- character()

  if (any(info$scheme == "free"))
    parts <- c(parts, sprintf("%s permutations", permutations))

  if (any(info$scheme == "within"))
    parts <- c(parts, sprintf("%s permutations within %s", permutations, strata_name))

  # A single between-strata test: report its exact count.
  if (nrow(info) == 1 && info$scheme == "between")
    return (perm_scheme_text("between", info$n, info$exact, strata_name))

  b <- info[info$scheme == "between", , drop = FALSE]
  if (nrow(b) > 0) {
    detail <- if (all(b$exact) && length(unique(b$n)) == 1) {
      sprintf("all %s arrangements, exact", b$n[[1]])
    } else if (any(b$exact)) {
      sprintf("%s permutations, or all arrangements when fewer", permutations)
    } else {
      sprintf("%s permutations", permutations)
    }
    parts <- c(parts, sprintf("whole %s strata permuted (between-strata variable; %s)", strata_name, detail))
  }

  paste(parts, collapse = "; ")
}


# A perm_plan() as a one-row data.frame, for perm_scheme_summary().
perm_info <- function (plan) {
  data.frame(scheme = plan$scheme, n = plan$n, exact = isTRUE(plan$exact))
}




#________________________________________________________
# Permutation test for a variable that is constant within each stratum:
# the stratum-level values are permuted among whole strata.
#
# Every arrangement only changes which strata share a value, so the
# statistic depends on the data only through stratum x stratum block sums,
# which are computed once:
#
#   adonis2 - pseudo-F = [tr(HG) / (p - 1)] / [(tr(G) - tr(HG)) / (n - p)]
#             and R2 = tr(HG) / tr(G), where G = -1/2 J D^2 J is the
#             Gower-centred matrix of squared distances and H the hat
#             matrix of model.matrix(~ groups). For a factor,
#             tr(HG) = sum_g S_g / n_g (S_g: sum of G over samples in group
#             g); for a numeric covariate x, tr(HG) = sum(G) / n +
#             xc' G xc / xc' xc with xc = x - mean(x). Identical to vegan's.
#   mantel  - Spearman's r between the distances and |x_i - x_j|: ranks of
#             the distances are fixed; ranks of |x_i - x_j| only depend on
#             the pair of strata, so they are computed with weights.
#
# Returns a one-row tibble.
#________________________________________________________
between_strata_test <- function (dm, grouping, blocks, test, permutations) {

  EPS    <- sqrt(.Machine$double.eps)
  blocks <- factor(blocks)
  idx    <- as.integer(blocks)
  n      <- attr(dm, 'Size')
  n_s    <- nlevels(blocks)

  if (!is.numeric(grouping)) grouping <- droplevels(as.factor(grouping))

  if (n_s < 2)
    cli_abort("At least two strata are needed to permute whole strata.")

  # One value per stratum (constant within strata), and stratum sizes.
  vals  <- grouping[match(seq_len(n_s), idx)]
  sizes <- tabulate(idx, n_s)

  # Block sums of an n x n matrix over pairs of strata.
  block_sums <- function (m) unname(rowsum(t(rowsum(m, idx, reorder = TRUE)), idx, reorder = TRUE))


  if (test == "adonis2") {

    D2  <- as.matrix(dm)^2
    rm  <- rowMeans(D2)
    G   <- -0.5 * (D2 - rm - rep(rm, each = n) + mean(D2))
    trG <- sum(diag(G))
    B   <- block_sums(G)
    totG <- sum(B)

    if (is.factor(vals)) {

      L <- nlevels(vals)

      stat <- function (v) {
        K    <- outer(as.integer(v), seq_len(L), `==`) * 1   # strata x groups
        S    <- colSums(K * (B %*% K))
        n_g  <- colSums(K * sizes)
        trHG <- sum(S / n_g)
        c((trHG / (L - 1)) / ((trG - trHG) / (n - L)), trHG / trG)
      }

    } else {

      stat <- function (v) {
        xc   <- v - sum(sizes * v) / n
        ss   <- sum(sizes * xc^2)
        trHG <- totG / n + drop(crossprod(xc, B %*% xc)) / ss
        c((trHG / 1) / ((trG - trHG) / (n - 2)), trHG / trG)
      }
    }

  } else if (test == "mantel") {

    N   <- n * (n - 1) / 2
    rx  <- rank(as.vector(dm)) - (N + 1) / 2
    Sx  <- sum(rx^2)

    Rm  <- matrix(0, n, n)
    Rm[lower.tri(Rm)] <- rx
    Bs  <- block_sums(Rm)
    Rs  <- Bs + t(Bs)
    diag(Rs) <- diag(Bs)

    Cm  <- outer(sizes, sizes)
    diag(Cm) <- sizes * (sizes - 1) / 2

    LT  <- lower.tri(Cm, diag = TRUE)
    Rs  <- Rs[LT]
    Cw  <- Cm[LT]

    # Ranks (average ties) of values y, each standing for w observations.
    wrank <- function (y, w) {
      u   <- sort(unique(y))
      k   <- match(y, u)
      wu  <- tabulate_w(k, w, length(u))
      cum <- cumsum(wu) - wu
      (cum + (wu + 1) / 2)[k]
    }

    stat <- function (v) {
      y  <- sqrt(outer(v, v, `-`)^2)[LT]        # as stats::dist() computes it
      ry <- wrank(y, Cw) - (N + 1) / 2
      c(sum(ry * Rs) / sqrt(Sx * sum(Cw * ry^2)))
    }

  } else {
    cli_abort("`test = \"{test}\"` can't permute whole strata; use `test = \"adonis2\"`.")
  }


  obs <- stat(vals)


  #________________________________________________________
  # Distinct arrangements of the stratum-level values: enumerate them all
  # when there are no more than `permutations`.
  #________________________________________________________
  codes  <- match(vals, unique(vals))
  arr    <- stratum_arrangements(vals, permutations)
  counts <- arr$counts
  exact  <- arr$exact

  if (is.factor(vals) && any(counts < 2))
    cli_warn(c(
      '!' = "Some `groups` values occur in only one stratum ({.val {as.character(unique(vals)[counts < 2])}}).",
      'i' = "With whole strata permuted, the null distribution is very coarse." ))

  if (arr$log_n < log(20)) {
    n_dist <- round(exp(arr$log_n))
    cli_warn(c(
      '!' = "Only {n_dist} distinct arrangement{?s} of `groups` among {n_s} strata.",
      'i' = "The smallest attainable p-value is {format(1 / n_dist, digits = 2)}." ))
  }

  if (exact) {

    arr   <- multiset_arrangements(counts)
    uvals <- vals[match(seq_along(counts), codes)]
    perm  <- vapply(seq_len(nrow(arr)), FUN.VALUE = 0, function (i) stat(uvals[arr[i, ]])[[1]])
    p_val <- mean(perm >= obs[[1]] - EPS)

  } else {

    perm  <- vapply(seq_len(permutations), FUN.VALUE = 0, function (i) stat(vals[sample.int(n_s)])[[1]])
    p_val <- (1 + sum(perm >= obs[[1]] - EPS)) / (1 + permutations)
  }


  if (test == "adonis2") {
    tibble::tibble(.n = n, .stat = obs[[1]], .r.sqr = obs[[2]], .p.val = p_val)
  } else {
    tibble::tibble(.n = n, .stat = obs[[1]], .p.val = p_val)
  }
}



# Weighted counts per category: sum of w for each value of k in 1..n.
tabulate_w <- function (k, w, n) {
  out <- numeric(n)
  s   <- rowsum(w, k, reorder = TRUE)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}



#________________________________________________________
# All distinct arrangements of a multiset with the given `counts` of each
# value (value k appears counts[k] times). One row per arrangement; the
# entries are value indices.
#________________________________________________________
multiset_arrangements <- function (counts) {

  n   <- sum(counts)
  res <- matrix(0L, nrow = 1, ncol = n)

  for (k in seq_along(counts)[-length(counts)]) {
    res <- do.call(rbind, lapply(seq_len(nrow(res)), function (i) {
      free <- which(res[i, ] == 0L)
      ch   <- utils::combn(length(free), counts[[k]])
      out  <- matrix(res[i, ], nrow = ncol(ch), ncol = n, byrow = TRUE)
      out[cbind(rep(seq_len(ncol(ch)), each = counts[[k]]), free[as.vector(ch)])] <- k
      out
    }))
  }

  res[res == 0L] <- length(counts)
  return (res)
}
