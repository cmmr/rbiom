#' Calculate PCoA and other ordinations, including taxa biplots and statistics.
#'
#' The biplot parameters (\code{taxa}, \code{unc}, \code{p.top}, and
#' \code{p.adj}) only only have an effect when \code{rank} is not `NULL`.
#'
#' Distances are computed from `biom` as given (not relativized first), so
#' they are identical to [bdiv_distmat()] on the same object. Relative
#' abundances are used only for the taxa biplot.
#'
#' @inherit documentation_rank.NULL
#' @inherit documentation_default
#' @inheritParams distmat_ord_table
#'
#' @family beta_diversity
#' @family ordination
#'
#' @param stat.by   The categorical or numeric metadata field over which
#'        statistics should be calculated. Default: `NULL`
#'
#' @param test    Permutational test for the sample groups (`stat.by`).
#'        Options are `"adonis2"` (PERMANOVA; the default), `"mrpp"`,
#'        `"betadisper"`, or `"none"`. See [distmat_stats()]. The taxa
#'        biplot always uses [vegan::envfit()].
#'
#' @param strata   Metadata field with a blocking variable (e.g. a subject
#'        identifier). For the sample statistics, a `stat.by` that varies
#'        within strata is permuted within each stratum, and one that is
#'        constant within every stratum is permuted at the stratum level
#'        (whole strata exchange their values); see the Details of
#'        [distmat_stats()]. For the taxa biplot, taxa abundances are
#'        permuted within each stratum. Default: `NULL`
#'
#' @param seed   Random seed for the permutation tests and the stochastic
#'        ordinations (NMDS, tSNE, UMAP). Default: `0`
#'
#' @param ...  Additional arguments to pass on to [stats::cmdscale()],
#'        [vegan::metaMDS()], [Rtsne::Rtsne()], or [uwot::umap()].
#'
#' @return A data.frame with columns \code{.sample},
#'         \code{.bdiv}, \code{.ord}, \code{.x}, \code{.y}, and (optionally)
#'         \code{.z}. Any columns given by \code{md}, \code{split.by}, and
#'         \code{stat.by} are included as well.
#'         If \code{stat.by} is given, then \code{$stats} and
#'         \code{$stats$code} are set (see [distmat_stats()]).
#'         If \code{rank} is given, then \code{$taxa_coords},
#'         \code{$taxa_stats}, and \code{$taxa_stats$code} are set:
#'         `$taxa_stats` has one row per taxon and ordination with
#'         `.r.sqr` (the squared correlation of the taxon's relative
#'         abundance with the ordination axes, from [vegan::envfit()]),
#'         `.p.val` (its permutation p-value), and `.adj.p` (adjusted
#'         across all tested taxa, before any `p.top` filtering);
#'         `$taxa_coords` holds the arrow for each taxon (`.x0`/`.y0` at
#'         the centroid of the samples, `.x`/`.y` at its tip, scaled to
#'         fit the sample cloud) and `.size`, its mean relative abundance.
#'         The `'axes'` attribute gives the axis labels (see
#'         [distmat_ord_table()]).
#'
#' @export
#' @examples
#'     library(rbiom)
#'
#'     ord <- bdiv_ord_table(hmp50, "bray", "pcoa", stat.by="Body Site", rank="g")
#'     head(ord)
#'
#'     ord$stats
#'
#'     ord$taxa_stats
#'
#'

bdiv_ord_table <- function (
    biom, bdiv = "bray", ord = "PCoA", weighted = NULL, md = NULL, k = 2,
    stat.by = NULL, split.by = NULL, strata = NULL, tree = NULL,
    test = "adonis2", seed = 0, permutations = 999, rank = NULL, taxa = 6,
    p.top = Inf, p.adj = 'fdr', unc = "singly",
    alpha = 0.5, cpus = n_cpus(), ...) {

  biom <- as_rbiom(biom)
  validate_tree(null_ok = TRUE)

  params <- eval_envir(environment(), ...)
  cmd    <- current_cmd('bdiv_ord_table', cmd_params(params))
  schema_check('bdiv_ord_table', params)
  dots   <- params$.dots


  #________________________________________________________
  # Validate and restructure user's arguments.
  #________________________________________________________
  validate_ord(max = Inf)
  validate_bdiv(multiple = TRUE)

  validate_biom_field('md',       null_ok = TRUE, max = Inf)
  validate_biom_field('stat.by',  null_ok = TRUE)
  validate_biom_field('split.by', null_ok = TRUE, col_type = "cat", max = Inf)
  validate_biom_field('strata',   null_ok = TRUE, col_type = "cat")

  validate_var_choices('test', CHOICES$ord.test)
  validate_var_choices('p.adj', CHOICES$p.adj)
  validate_unc()
  validate_var_range('k', range = c(2, 3), n = 1, int = TRUE)
  validate_var_range('p.top', range = arg_range(ARGS$p.top), n = 1)
  validate_var_range('permutations', int = TRUE, range = arg_range(ARGS$permutations), n = 1)
  validate_seed()
  validate_alpha()
  validate_cpus()

  if (biom$n_samples < 3)
    cli_abort("At least three samples are needed for an ordination; `biom` has {biom$n_samples}.")


  #________________________________________________________
  # Drop samples with NA values in the grouping fields.
  #________________________________________________________
  fields <- unique(c(stat.by, split.by, strata))
  if (length(fields) > 0) {
    keep <- stats::complete.cases(biom$metadata[, fields, drop = FALSE])
    if (!any(keep)) cli_abort("All samples have NA values in {.field {fields}}.")
    if (!all(keep)) {
      cli::cli_inform(c('i' = "Dropped {sum(!keep)} of {length(keep)} sample{?s} with NA values in {.field {fields}}."))
      biom <- na.omit(biom, fields = fields)
    }
  }


  #________________________________________________________
  # Biplot - determine taxonomic rank
  #________________________________________________________
  if (is.character(taxa)) {
    rank <- names(which.max(lapply(biom$taxonomy, function (x) sum(x %in% taxa))))
  } else {
    validate_rank(max = Inf, null_ok = TRUE)
  }


  #________________________________________________________
  # See if this result is already in the cache.
  #________________________________________________________
  cache_file    <- get_cache_file('bdiv_ord_table', params)
  sample_coords <- get_cache_value(cache_file)

  if (is.null(sample_coords)) {

    #________________________________________________________
    # Everything for one stratum and one bdiv metric. The code texts
    # further down mirror this function.
    #________________________________________________________
    results <- blply(
      biom   = biom,
      vars   = split.by,
      iters  = list(bdiv = bdiv),
      prefix = TRUE,
      FUN    = function (b, bdiv) {

        dm <- bdiv_distmat(b, bdiv = bdiv, tree = tree, alpha = alpha, cpus = cpus)


        #________________________________________________________
        # Sample Ordination(s). x/y/z coordinates.
        #________________________________________________________
        coords <- do.call(distmat_ord_table, c(list(dm = dm, ord = ord, k = k, seed = seed), dots))
        res    <- list(coords = coords, axes = attr(coords, 'axes', exact = TRUE))


        #________________________________________________________
        # Statistics: Sample Groups vs Distance Matrix
        #________________________________________________________
        if (!is.null(stat.by) && test != "none") {

          res$stats <- tryCatch(
            error = function (e) { warning(conditionMessage(e), call. = FALSE); tibble(.p.val = NA_real_) },
            expr  = distmat_stats(
              dm           = dm,
              groups       = pull(b, stat.by),
              strata       = if (!is.null(strata)) pull(b, strata),
              test         = test,
              seed         = seed,
              permutations = permutations ))

          # The permutation scheme, for the table header.
          if (hasName(res$stats, '.n'))
            res$plan <- perm_plan(pull(b, stat.by), if (!is.null(strata)) pull(b, strata), permutations)
        }


        #________________________________________________________
        # Biplot: each taxon's abundance vector fitted onto the
        # ordination with vegan::envfit().
        #________________________________________________________
        if (!is.null(rank)) {

          fits <- dplyr::bind_rows(lapply(rank, function (rank) {

            taxa_mtx <- taxa_matrix(b, rank = rank, taxa = taxa, unc = unc, transform = "percent")

            dplyr::bind_rows(lapply(levels(coords$.ord), function (ord) {

              xy  <- as.matrix(coords[coords$.ord == ord, intersect(c('.x', '.y', '.z'), names(coords))])
              rownames(xy) <- as.character(coords$.sample[coords$.ord == ord])
              env <- t(taxa_mtx)[rownames(xy), , drop = FALSE]

              perms <- if (is.null(strata)) permutations else
                permute::how(nperm = permutations, blocks = factor(pull(b, strata)[rownames(xy)]))

              fit <- with_seed(seed, vegan::envfit(xy, env, permutations = perms))

              tibble(
                .rank  = rank,
                .ord   = ord,
                .taxa  = names(fit$vectors$r),
                .r.sqr = unname(fit$vectors$r),
                .p.val = unname(fit$vectors$pvals),
                .size  = unname(rowMeans(taxa_mtx)[names(fit$vectors$r)]),
                biplot_arrows(fit, xy) )
            }))
          }))

          coord_cols      <- grep("^\\.[xyz]0?$", names(fits), value = TRUE)
          res$taxa_stats  <- fits[, c('.rank', '.ord', '.taxa', '.r.sqr', '.p.val')]
          res$taxa_coords <- fits[, c('.rank', '.ord', '.taxa', '.size', coord_cols)]
        }

        return (res)
      })


    #________________________________________________________
    # Combine the per-stratum results, prefixed by the stratum labels
    # (split.by fields and .bdiv).
    #________________________________________________________
    labels <- attr(results, 'split_labels', exact = TRUE)

    combine <- function (what) {
      rows <- lapply(seq_along(results), function (i) {
        x <- results[[i]][[what]]
        if (is.null(x) || nrow(x) == 0) return (NULL)
        cbind(labels[rep(i, nrow(x)), , drop = FALSE], as.data.frame(x))
      })
      res <- dplyr::bind_rows(rows)
      if (nrow(res) == 0) return (NULL)
      as_rbiom_tbl(res)
    }

    sample_coords <- combine('coords')
    axes          <- combine('axes')
    sample_stats  <- combine('stats')
    taxa_stats    <- combine('taxa_stats')
    taxa_coords   <- combine('taxa_coords')

    sample_coords[['.bdiv']] %<>% factor(levels = bdiv)
    axes[['.bdiv']]          %<>% factor(levels = bdiv)

    for (i in unique(c(split.by, stat.by, md)))
      if (!hasName(sample_coords, i))
        sample_coords[[i]] <- unname(pull(biom, i)[as.character(sample_coords[['.sample']])])


    #________________________________________________________
    # Code shared by the stats tables.
    #________________________________________________________
    # The stats tables always carry `.bdiv`, so the code always iterates
    # over the metric(s) with bdply().
    dm_args <- paste(collapse = ", ", c(
      "bdiv = bdiv",
      if (!is.null(tree))  "tree = tree",
      if (!eq(alpha, 0.5)) sprintf("alpha = %s", deparse_vec(alpha)) ))

    preamble <- c(
      if (!is.null(tree)) sprintf("# tree <- %s   # not reproducible from code", obj_repr(tree)),
      if (length(fields) > 0 && !all(keep)) sprintf("biom <- na.omit(biom, fields = %s)", as.args(list(fields))) )


    #________________________________________________________
    # Sample statistics.
    #________________________________________________________
    if (!is.null(sample_stats)) {

      sample_stats[['.adj.p']] <- p.adjust(sample_stats[['.p.val']], method = p.adj)

      inner <- c(
        sprintf("dm <- bdiv_distmat(biom, %s)", dm_args),
        "tryCatch(",
        "  error = function (e) { warning(conditionMessage(e), call. = FALSE); tibble::tibble(.p.val = NA_real_) },",
        "  expr  = distmat_stats(",
        "    dm           = dm,",
        sprintf("    groups       = pull(biom, %s),", double_quote(stat.by)),
        if (!is.null(strata)) sprintf("    strata       = pull(biom, %s),", double_quote(strata)),
        sprintf("    test         = %s,", double_quote(test)),
        sprintf("    seed         = %s,", deparse_vec(seed)),
        sprintf("    permutations = %s ))", deparse_vec(permutations)) )

      attr(sample_stats, 'code') <- paste(collapse = "\n", c(
        preamble,
        if (length(preamble) > 0) "",
        bdply_code(inner, split.by, bdiv, iterate = TRUE),
        sprintf("stats <- dplyr::mutate(stats, .adj.p = p.adjust(.p.val, %s))", double_quote(p.adj)) ))

      schemes <- do.call(rbind, lapply(results, function (r) if (!is.null(r$plan)) perm_info(r$plan)))

      attr(sample_stats, 'tbl_sum') <- c('Test' = sprintf(
        "%s ~ %s; %s.",
        switch(test, adonis2 = "PERMANOVA", mrpp = "MRPP (vegan::mrpp)", betadisper = "Dispersion (vegan::betadisper)"),
        coan(stat.by),
        if (is.null(schemes)) sprintf("%s permutations", permutations) else
          perm_scheme_summary(schemes, coan(strata %||% "strata"), permutations) ))

      sample_stats <- prov_attach(sample_stats, cmd = paste0(cmd, "$stats"), biom = biom)
    }


    #________________________________________________________
    # Taxa statistics: adjust over all tested taxa, then apply p.top
    # within each stratum / metric / ordination / rank.
    #________________________________________________________
    if (!is.null(taxa_stats)) {

      taxa_stats[['.adj.p']] <- p.adjust(taxa_stats[['.p.val']], method = p.adj)

      if (is.finite(p.top)) {

        by <- intersect(c(split.by, '.bdiv', '.ord', '.rank'), names(taxa_stats))

        taxa_stats <- split_apply(taxa_stats, by, function (x) {
          if (p.top >= 1) { x[head(order(x[['.p.val']]), p.top), setdiff(names(x), by), drop = FALSE]
          } else          { x[which(x[['.adj.p']] <= p.top), setdiff(names(x), by), drop = FALSE] }
        })

        keys        <- c(by, '.taxa')
        taxa_coords <- dplyr::semi_join(taxa_coords, taxa_stats, by = keys)
        taxa_stats  <- as_rbiom_tbl(taxa_stats)
      }


      dots_args <- if (length(dots) > 0) paste0(", ", as.args(dots))

      inner <- c(
        sprintf("dm     <- bdiv_distmat(biom, %s)", dm_args),
        sprintf("coords <- distmat_ord_table(dm, ord = %s, k = %s, seed = %s%s)", as.args(list(ord)), deparse_vec(k), deparse_vec(seed), if (is.null(dots_args)) "" else dots_args),
        sprintf("dplyr::bind_rows(lapply(%s, function (rank) {", as.args(list(rank))),
        sprintf("  taxa_mtx <- taxa_matrix(biom, rank = rank, taxa = %s, unc = %s, transform = 'percent')", as.args(list(taxa)), double_quote(unc)),
        "  dplyr::bind_rows(lapply(levels(coords$.ord), function (ord) {",
        "    xy  <- as.matrix(coords[coords$.ord == ord, intersect(c('.x', '.y', '.z'), names(coords))])",
        "    rownames(xy) <- as.character(coords$.sample[coords$.ord == ord])",
        "    env <- t(taxa_mtx)[rownames(xy), , drop = FALSE]",
        if (!is.null(strata)) sprintf("    blocks <- factor(pull(biom, %s)[rownames(xy)])", double_quote(strata)),
        sprintf("    set.seed(%s)", deparse_vec(seed)),
        sprintf("    fit <- vegan::envfit(xy, env, permutations = %s)", if (is.null(strata)) deparse_vec(permutations) else sprintf("permute::how(nperm = %s, blocks = blocks)", deparse_vec(permutations))),
        "    tibble::tibble(",
        "      .rank  = rank,",
        "      .ord   = ord,",
        "      .taxa  = names(fit$vectors$r),",
        "      .r.sqr = unname(fit$vectors$r),",
        "      .p.val = unname(fit$vectors$pvals) )",
        "  }))",
        "}))" )

      attr(taxa_stats, 'code') <- paste(collapse = "\n", c(
        preamble,
        if (length(preamble) > 0) "",
        bdply_code(inner, split.by, bdiv, iterate = TRUE),
        sprintf("stats <- dplyr::mutate(stats, .adj.p = p.adjust(.p.val, %s))", double_quote(p.adj)),
        if (p.top < 1)          sprintf("stats <- stats[stats$.adj.p <= %s, ]", deparse_vec(p.top)),
        if (p.top >= 1 && is.finite(p.top)) sprintf("# p.top = %s: the %s most significant taxa per stratum, metric, ordination, and rank were kept.", deparse_vec(p.top), deparse_vec(p.top)) ))

      attr(taxa_stats, 'tbl_sum') <- c('Test' = sprintf(
        "vegan::envfit ~ taxa; %s permutations%s.", permutations,
        if (is.null(strata)) "" else sprintf(" within %s", coan(strata)) ))

      taxa_stats <- prov_attach(taxa_stats, cmd = paste0(cmd, "$taxa_stats"), biom = biom)
    }


    attr(sample_coords, 'axes')        <- axes
    attr(sample_coords, 'stats')       <- sample_stats
    attr(sample_coords, 'taxa_coords') <- taxa_coords
    attr(sample_coords, 'taxa_stats')  <- taxa_stats

    set_cache_value(cache_file, sample_coords)
  }


  return (prov_attach(sample_coords, cmd, biom = params$biom))
}




#________________________________________________________
# Arrow coordinates for a biplot: envfit's vectors (unit direction
# scaled by sqrt(r2)) multiplied so that the longest arrow spans 75% of
# the sample cloud's extent from its centroid (after vegan::ordiArrowMul,
# but using the data range rather than a graphics device).
#________________________________________________________
biplot_arrows <- function (fit, xy) {

  arrows <- vegan::scores(fit, display = "vectors")
  arrows <- arrows[, seq_len(ncol(xy)), drop = FALSE]
  origin <- colMeans(xy)

  extent <- rbind(apply(xy, 2, min) - origin, apply(xy, 2, max) - origin)
  reach  <- rbind(apply(arrows, 2, min), apply(arrows, 2, max))
  ratios <- extent / reach
  ratios <- ratios[is.finite(ratios) & ratios > 0]
  mul    <- if (length(ratios) == 0) 1 else 0.75 * min(ratios)

  tips <- sweep(arrows * mul, 2, origin, `+`)

  axis_names <- c('.x', '.y', '.z')[seq_len(ncol(xy))]

  res <- as.data.frame(signif(tips, 6))
  names(res) <- axis_names
  for (i in seq_along(axis_names))
    res[[paste0(axis_names[[i]], '0')]] <- signif(origin[[i]], 6)

  rownames(res) <- NULL
  return (res)
}
