#________________________________________________________
# Methods text.
#
# `x$methods` on a table, statistics table, plot, or rbiom object gives a
# short paragraph, written for the methods section of a manuscript, that
# says how the result was produced, followed by the references it cites.
#
# Nothing is recomputed. The text is assembled from what the result
# already carries:
#
#   * `history`   - the code lines of the rbiom object's history (import,
#                   subsetting, rarefaction, ...);
#   * `cmd`       - the rbiom call that made the result. It is parsed and
#                   matched to the function's formals, so both the
#                   arguments given and the defaults are known;
#   * `data_cmd`  - the call that made the plotted / tested data;
#   * `tbl_sum`   - on statistics tables: the model or test that was run;
#   * `axes`      - on ordinations: the NMDS stress;
#   * `biom_info` - an attribute of `history`: sample and OTU counts and
#                   input_table_caption()'s description of the abundance
#                   table (attached by prov_attach()).
#
# methods_text() never fails: a result it can't describe gets the
# software sentence only.
#________________________________________________________


#________________________________________________________
# References, keyed by method. `cite` is the in-text form.
#________________________________________________________
CITATIONS <- list(

  # Ordination.
  gower1966 = list(
    cite = "Gower 1966",
    ref  = "Gower JC (1966). Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika 53(3-4):325-338." ),
  kruskal1964 = list(
    cite = "Kruskal 1964",
    ref  = "Kruskal JB (1964). Multidimensional scaling by optimizing goodness of fit to a nonmetric hypothesis. Psychometrika 29(1):1-27." ),
  maaten2008 = list(
    cite = "van der Maaten and Hinton 2008",
    ref  = "van der Maaten L, Hinton G (2008). Visualizing data using t-SNE. Journal of Machine Learning Research 9:2579-2605." ),
  mcinnes2018 = list(
    cite = "McInnes et al. 2018",
    ref  = "McInnes L, Healy J, Melville J (2018). UMAP: Uniform Manifold Approximation and Projection for dimension reduction. arXiv:1802.03426." ),

  # Tests on distance matrices.
  anderson2001 = list(
    cite = "Anderson 2001",
    ref  = "Anderson MJ (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecology 26(1):32-46." ),
  anderson2006 = list(
    cite = "Anderson 2006",
    ref  = "Anderson MJ (2006). Distance-based tests for homogeneity of multivariate dispersions. Biometrics 62(1):245-253." ),
  mielke2001 = list(
    cite = "Mielke and Berry 2001",
    ref  = "Mielke PW, Berry KJ (2001). Permutation Methods: A Distance Function Approach. Springer, New York." ),
  mantel1967 = list(
    cite = "Mantel 1967",
    ref  = "Mantel N (1967). The detection of disease clustering and a generalized regression approach. Cancer Research 27(2):209-220." ),

  # Tests on one response.
  kruskal1952 = list(
    cite = "Kruskal and Wallis 1952",
    ref  = "Kruskal WH, Wallis WA (1952). Use of ranks in one-criterion variance analysis. Journal of the American Statistical Association 47(260):583-621." ),
  wilcoxon1945 = list(
    cite = "Wilcoxon 1945",
    ref  = "Wilcoxon F (1945). Individual comparisons by ranking methods. Biometrics Bulletin 1(6):80-83." ),
  mann1947 = list(
    cite = "Mann and Whitney 1947",
    ref  = "Mann HB, Whitney DR (1947). On a test of whether one of two random variables is stochastically larger than the other. Annals of Mathematical Statistics 18(1):50-60." ),
  wood2011 = list(
    cite = "Wood 2011",
    ref  = "Wood SN (2011). Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. Journal of the Royal Statistical Society: Series B 73(1):3-36." ),

  # Differential abundance.
  zhou2022 = list(
    cite = "Zhou et al. 2022",
    ref  = "Zhou H, He K, Chen J, Zhang X (2022). LinDA: linear models for differential abundance analysis of microbiome compositional data. Genome Biology 23:95." ),
  kuznetsova2017 = list(
    cite = "Kuznetsova et al. 2017",
    ref  = "Kuznetsova A, Brockhoff PB, Christensen RHB (2017). lmerTest package: tests in linear mixed effects models. Journal of Statistical Software 82(13):1-26." ),

  # Multiple comparisons.
  benjamini1995 = list(
    cite = "Benjamini and Hochberg 1995",
    ref  = "Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society: Series B 57(1):289-300." ),
  benjamini2001 = list(
    cite = "Benjamini and Yekutieli 2001",
    ref  = "Benjamini Y, Yekutieli D (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29(4):1165-1188." ),
  holm1979 = list(
    cite = "Holm 1979",
    ref  = "Holm S (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6(2):65-70." ),
  hochberg1988 = list(
    cite = "Hochberg 1988",
    ref  = "Hochberg Y (1988). A sharper Bonferroni procedure for multiple tests of significance. Biometrika 75(4):800-802." ),
  hommel1988 = list(
    cite = "Hommel 1988",
    ref  = "Hommel G (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika 75(2):383-386." ),

  # Graphics.
  lex2014 = list(
    cite = "Lex et al. 2014",
    ref  = "Lex A, Gehlenborg N, Strobelt H, Vuillemot R, Pfister H (2014). UpSet: visualization of intersecting sets. IEEE Transactions on Visualization and Computer Graphics 20(12):1983-1992." ),
  wickham2016 = list(
    cite = "Wickham 2016",
    ref  = "Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag, New York." )
)


# A reference for each diversity metric, keyed by ecodive's metric id:
# the first reference on the metric's page in ecodive's manual (its
# \references section), in the style of CITATIONS. Kept in step with
# the installed ecodive by a test in test-methods.r.
METRIC_REFS <- list(
  ace = list(
    cite = "Chao and Lee 1992",
    ref  = "Chao A, Lee SM (1992). Estimating the number of classes via sample coverage. Journal of the American Statistical Association 87(417):210-217. doi:10.1080/01621459.1992.10475194" ),
  aitchison = list(
    cite = "Aitchison 1986",
    ref  = "Aitchison J (1986). The statistical analysis of compositional data. Chapman and Hall. doi:10.1002/bimj.4710300705" ),
  berger = list(
    cite = "Berger and Parker 1970",
    ref  = "Berger WH, Parker FL (1970). Diversity of planktonic foraminifera in deep-sea sediments. Science 168(3937):1345-1347. doi:10.1126/science.168.3937.1345" ),
  bhattacharyya = list(
    cite = "Bhattacharyya 1943",
    ref  = "Bhattacharyya A (1943). On a measure of divergence between two statistical populations defined by their probability distributions. Bulletin of the Calcutta Mathematical Society 35:99-109." ),
  bray = list(
    cite = "Bray and Curtis 1957",
    ref  = "Bray JR, Curtis JT (1957). An ordination of the upland forest communities of southern Wisconsin. Ecological Monographs 27(4):325-349. doi:10.2307/1942268" ),
  brillouin = list(
    cite = "Brillouin 1956",
    ref  = "Brillouin L (1956). Science and information theory. Academic Press." ),
  canberra = list(
    cite = "Lance and Williams 1966",
    ref  = "Lance GN, Williams WT (1966). Computer programs for hierarchical polythetic classification (\"similarity analyses\"). The Computer Journal 9(1):60-64. doi:10.1093/comjnl/9.1.60" ),
  chao1 = list(
    cite = "Chao 1984",
    ref  = "Chao A (1984). Nonparametric estimation of the number of classes in a population. Scandinavian Journal of Statistics 11:265-270." ),
  chebyshev = list(
    cite = "Cantrell 2000",
    ref  = "Cantrell CD (2000). Modern mathematical methods for physicists and engineers. Cambridge University Press." ),
  chord = list(
    cite = "Orl\u00f3ci 1967",
    ref  = "Orl\u00f3ci L (1967). An agglomerative method for classification of plant communities. Journal of Ecology 55(1):193-206. doi:10.2307/2257725" ),
  clark = list(
    cite = "Clark 1952",
    ref  = "Clark PJ (1952). An extension of the coefficient of divergence for use with multiple characters. Copeia 1952(2):61-64. doi:10.2307/1438598" ),
  sorensen = list(
    cite = "S\u00f8rensen 1948",
    ref  = "S\u00f8rensen T (1948). A method of establishing groups of equal amplitude in plant sociology based on similarity of species content. Kongelige Danske Videnskabernes Selskab, Biologiske Skrifter 5:1-34." ),
  divergence = list(
    cite = "Cha 2007",
    ref  = "Cha S-H (2007). Comprehensive survey on distance/similarity measures between probability density functions. International Journal of Mathematical Models and Methods in Applied Sciences 1(4):300-307." ),
  euclidean = list(
    cite = "Legendre and Legendre 2012",
    ref  = "Legendre P, Legendre L (2012). Numerical ecology. Elsevier." ),
  faith = list(
    cite = "Faith 1992",
    ref  = "Faith DP (1992). Conservation evaluation and phylogenetic diversity. Biological Conservation 61(1):1-10. doi:10.1016/0006-3207(92)91201-3" ),
  fisher = list(
    cite = "Fisher et al. 1943",
    ref  = "Fisher RA, Corbet AS, Williams CB (1943). The relation between the number of species and the number of individuals in a random sample of an animal population. Journal of Animal Ecology 12:42-58. doi:10.2307/1411" ),
  generalized_unifrac = list(
    cite = "Chen et al. 2012",
    ref  = "Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD, ..., Li H (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. Bioinformatics 28(16):2106-2113. doi:10.1093/bioinformatics/bts342" ),
  simpson = list(
    cite = "Simpson 1949",
    ref  = "Simpson EH (1949). Measurement of diversity. Nature 163:688. doi:10.1038/163688a0" ),
  gower = list(
    cite = "Gower 1971",
    ref  = "Gower JC (1971). A general coefficient of similarity and some of its properties. Biometrics 27(4):857-871. doi:10.2307/2528823" ),
  hamming = list(
    cite = "Hamming 1950",
    ref  = "Hamming RW (1950). Error detecting and error correcting codes. Bell System Technical Journal 29(2):147-160. doi:10.1002/j.1538-7305.1950.tb00463.x" ),
  hellinger = list(
    cite = "Rao 1995",
    ref  = "Rao CR (1995). A review of canonical coordinates and an alternative to correspondence analysis using Hellinger distance. Q\u00fcestii\u00f3 19:23-63." ),
  horn = list(
    cite = "Horn 1966",
    ref  = "Horn HS (1966). Measurement of \"overlap\" in comparative ecological studies. The American Naturalist 100(914):419-424. doi:10.1086/282436" ),
  inv_simpson = list(
    cite = "Simpson 1949",
    ref  = "Simpson EH (1949). Measurement of diversity. Nature 163:688. doi:10.1038/163688a0" ),
  jaccard = list(
    cite = "Jaccard 1912",
    ref  = "Jaccard P (1912). The distribution of the flora in the alpine zone. New Phytologist 11(2):37-50. doi:10.1111/j.1469-8137.1912.tb05611.x" ),
  jensen = list(
    cite = "Endres and Schindelin 2003",
    ref  = "Endres DM, Schindelin JE (2003). A new metric for probability distributions. IEEE Transactions on Information Theory 49(7):1858-1860. doi:10.1109/TIT.2003.813506" ),
  jsd = list(
    cite = "Lin 1991",
    ref  = "Lin J (1991). Divergence measures based on the Shannon entropy. IEEE Transactions on Information Theory 37(1):145-151. doi:10.1109/18.61115" ),
  lorentzian = list(
    cite = "Cha 2007",
    ref  = "Cha S-H (2007). Comprehensive survey on distance/similarity measures between probability density functions. International Journal of Mathematical Models and Methods in Applied Sciences 1(4):300-307." ),
  manhattan = list(
    cite = "Krause 1987",
    ref  = "Krause EF (1987). Taxicab geometry: An adventure in non-Euclidean geometry. Dover Publications." ),
  margalef = list(
    cite = "Margalef 1958",
    ref  = "Margalef R (1958). Information theory in ecology. General Systems 3:36-71." ),
  matusita = list(
    cite = "Matusita 1955",
    ref  = "Matusita K (1955). Decision rules, based on the distance, for problems of fit, two samples, and estimation. The Annals of Mathematical Statistics 26(4):631-640. doi:10.1214/aoms/1177728422" ),
  mcintosh = list(
    cite = "McIntosh 1967",
    ref  = "McIntosh RP (1967). An index of diversity and the relation of certain concepts to diversity. Ecology 48(3):392-404. doi:10.2307/1932674" ),
  menhinick = list(
    cite = "Menhinick 1964",
    ref  = "Menhinick EF (1964). A comparison of some species-individuals diversity indices applied to samples of field insects. Ecology 45(4):859-861. doi:10.2307/1934933" ),
  minkowski = list(
    cite = "Deza and Deza 2009",
    ref  = "Deza MM, Deza E (2009). Encyclopedia of distances. Springer." ),
  morisita = list(
    cite = "Morisita 1959",
    ref  = "Morisita M (1959). Measuring of interspecific association and similarity between communities. Memoirs of the Faculty of Science, Kyushu University, Series E (Biology) 3:65-80." ),
  motyka = list(
    cite = "Motyka 1947",
    ref  = "Motyka J (1947). O celach i metodach badan geobotanicznych. Annales Universitatis Mariae Curie-Sklodowska, Sectio C 3:1-168." ),
  normalized_unifrac = list(
    cite = "Lozupone et al. 2007",
    ref  = "Lozupone CA, Hamady M, Kelley ST, Knight R (2007). Quantitative and qualitative beta diversity measures lead to different insights into factors that structure microbial communities. Applied and Environmental Microbiology 73(5):1576-1585. doi:10.1128/AEM.01996-06" ),
  ochiai = list(
    cite = "Ochiai 1957",
    ref  = "Ochiai A (1957). Zoogeographic studies on the soleoid fishes found in Japan and its neighbouring regions. Bulletin of the Japanese Society of Scientific Fisheries 22:526-530." ),
  psym_chisq = list(
    cite = "Cha 2007",
    ref  = "Cha S-H (2007). Comprehensive survey on distance/similarity measures between probability density functions. International Journal of Mathematical Models and Methods in Applied Sciences 1(4):300-307." ),
  robust_aitchison = list(
    cite = "Martino et al. 2019",
    ref  = "Martino C, Morton JT, Marotz CA, Thompson LR, Tripathi A, Knight R, Zengler K (2019). A novel sparse compositional technique reveals microbial perturbations. mSystems 4(1):e00016-19. doi:10.1128/mSystems.00016-19" ),
  shannon = list(
    cite = "Shannon 1948",
    ref  = "Shannon CE (1948). A mathematical theory of communication. Bell System Technical Journal 27:379-423." ),
  soergel = list(
    cite = "Cha 2007",
    ref  = "Cha S-H (2007). Comprehensive survey on distance/similarity measures between probability density functions. International Journal of Mathematical Models and Methods in Applied Sciences 1(4):300-307." ),
  squared_chisq = list(
    cite = "Cha 2007",
    ref  = "Cha S-H (2007). Comprehensive survey on distance/similarity measures between probability density functions. International Journal of Mathematical Models and Methods in Applied Sciences 1(4):300-307." ),
  squared_chord = list(
    cite = "Legendre and Legendre 2012",
    ref  = "Legendre P, Legendre L (2012). Numerical ecology. Elsevier." ),
  squared_euclidean = list(
    cite = "Legendre and Legendre 2012",
    ref  = "Legendre P, Legendre L (2012). Numerical ecology. Elsevier." ),
  squares = list(
    cite = "Alroy 2018",
    ref  = "Alroy J (2018). Limits to species richness estimates based on subsampling. Paleobiology 44(2):177-194. doi:10.1017/pab.2017.38" ),
  topsoe = list(
    cite = "Topsoe 2000",
    ref  = "Topsoe F (2000). Some inequalities for information divergence and related measures of discrimination. IEEE Transactions on Information Theory 46(4):1602-1609. doi:10.1109/18.850703" ),
  unweighted_unifrac = list(
    cite = "Lozupone and Knight 2005",
    ref  = "Lozupone C, Knight R (2005). UniFrac: a new phylogenetic method for comparing microbial communities. Applied and Environmental Microbiology 71(12):8228-8235. doi:10.1128/AEM.71.12.8228-8235.2005" ),
  variance_adjusted_unifrac = list(
    cite = "Chang et al. 2011",
    ref  = "Chang Q, Luan Y, Sun F (2011). Variance adjusted weighted UniFrac: a powerful beta diversity measure for comparing communities based on phylogeny. BMC Bioinformatics 12:118. doi:10.1186/1471-2105-12-118" ),
  wave_hedges = list(
    cite = "Cha 2007",
    ref  = "Cha S-H (2007). Comprehensive survey on distance/similarity measures between probability density functions. International Journal of Mathematical Models and Methods in Applied Sciences 1(4):300-307." ),
  weighted_unifrac = list(
    cite = "Lozupone et al. 2007",
    ref  = "Lozupone CA, Hamady M, Kelley ST, Knight R (2007). Quantitative and qualitative beta diversity measures lead to different insights into factors that structure microbial communities. Applied and Environmental Microbiology 73(5):1576-1585. doi:10.1128/AEM.01996-06" ) )

# Metric references are cited through CITATIONS as "metric_<id>".
CITATIONS <- c(CITATIONS, stats::setNames(METRIC_REFS, paste0("metric_", names(METRIC_REFS))))


# R packages named in the software sentence: authors and a fallback
# title. The title and version are read from the installed package.
# `key` is a paper in CITATIONS that is cited for the package instead.
METHODS_PKGS <- list(
  rbiom    = list(authors = "Smith DP",                                title = "Integrated Analysis and Visualization of Microbiome Data"),
  ecodive  = list(authors = "Smith DP",                                title = "Parallel and Memory-Efficient Ecological Diversity Metrics"),
  vegan    = list(authors = "Oksanen J, Simpson GL, Blanchet FG, et al", title = "Community Ecology Package"),
  emmeans  = list(authors = "Lenth RV",                                title = "Estimated Marginal Means, aka Least-Squares Means"),
  mgcv     = list(key = 'wood2011'),
  lmerTest = list(key = 'kuznetsova2017'),
  ggplot2  = list(key = 'wickham2016'),
  Rtsne    = list(authors = "Krijthe JH",                              title = "T-Distributed Stochastic Neighbor Embedding using a Barnes-Hut Implementation"),
  uwot     = list(authors = "Melville J",                              title = "The Uniform Manifold Approximation and Projection (UMAP) Method for Dimensionality Reduction") )



#________________________________________________________
# The methods text for `obj`: an 'rbiom_methods' string holding the
# paragraph and the reference list, with the two parts also available
# as attr(, 'text') and attr(, 'references').
#________________________________________________________
methods_text <- function (obj) {
  tryCatch(
    expr  = methods_build(obj),
    error = function (e) methods_build(NULL) )
}


methods_build <- function (obj) {

  keys <- character()
  pkgs <- list()

  ctx <- list(

    # "Author Year" for one or more references, which are then listed.
    cite = function (...) {
      k    <- c(...)
      keys <<- union(keys, k)
      paste(vapply(k, function (i) CITATIONS[[i]]$cite, ''), collapse = "; ")
    },

    # Note that a package did part of the work.
    use = function (pkg, purpose) {
      pkgs[[pkg]] <<- union(pkgs[[pkg]], purpose)
      invisible(NULL)
    } )

  info <- if (is.null(obj)) list() else methods_info(obj)

  body <- c(
    methods_prep(info, ctx),
    methods_analysis(info, ctx) )

  if (isTRUE(info$is_plot)) ctx$use('ggplot2', "graphics")

  soft <- methods_software(pkgs, ctx)
  refs <- c(soft$refs, vapply(keys, function (k) CITATIONS[[k]]$ref, ''))
  refs <- unname(unique(refs))
  refs <- refs[order(tolower(refs), method = "radix")]

  text <- paste(c(body, soft$text), collapse = " ")
  full <- paste0(text, "\n\nReferences\n", paste(refs, collapse = "\n"))

  structure(full, text = text, references = refs, class = c('rbiom_methods', 'character'))
}



#' Print methods text: the paragraph and each reference wrapped to the
#' console width.
#'
#' @param x   An object of class `rbiom_methods`.
#' @param ...   Not used.
#'
#' @noRd
#' @keywords internal
#' @export
print.rbiom_methods <- function (x, ...) {
  cat(format.rbiom_methods(x, ...), sep = "\n")
  return (invisible(x))
}


#' @noRd
#' @keywords internal
#' @export
format.rbiom_methods <- function (x, width = getOption("width", 80L), ...) {

  text <- attr(x, 'text',       exact = TRUE) %||% as.character(x)
  refs <- attr(x, 'references', exact = TRUE)

  width <- max(20L, as.integer(width) - 1L)

  c(
    strwrap(text, width = width),
    if (length(refs) > 0) c(
      "",
      "References",
      unlist(lapply(refs, strwrap, width = width, exdent = 4L)) ))
}



#________________________________________________________
# Everything the text is built from.
#________________________________________________________
methods_info <- function (obj) {

  if (inherits(obj, 'rbiom'))
    return (list(
      is_biom = TRUE,
      history = methods_history(prov_lines(obj)),
      biom    = prov_biom_info(obj) ))

  get <- function (x, nm) if (is.null(x)) NULL else attr(x, nm, exact = TRUE)

  is_plot <- inherits(obj, c('ggplot', 'patchwork', 'gg'))

  data <- tryCatch(obj[['data']], error = function (e) NULL)
  if (!is.data.frame(data)) data <- NULL

  tbl_sum <- get(obj, 'tbl_sum')
  stats   <- get(obj, 'stats')
  if (is.null(stats) && any(c('Model', 'Test', 'Method') %in% names(tbl_sum))) stats <- obj

  history <- get(obj, 'history') %||% get(data, 'history')

  list(
    is_biom    = FALSE,
    is_plot    = is_plot,
    obj        = obj,
    data       = data,
    call       = methods_call(get(obj, 'cmd')),
    data_call  = methods_call(get(obj, 'data_cmd') %||% get(data, 'cmd')),
    history    = methods_history(history),
    biom       = get(history, 'biom_info'),
    axes       = get(obj, 'axes') %||% get(data, 'axes'),
    stats      = stats,
    is_stats   = identical(stats, obj),
    taxa_stats = get(obj, 'taxa_stats') )
}


#________________________________________________________
# Parse a recorded command, e.g. `adiv_table(biom, "shannon")`, into
# list(fn, fun, args, part). `args` are matched to the function's
# formals; values that aren't literals (`biom`, `tree`, placeholders)
# stay as language objects. `part` is "stats" for `fn(...)$stats`.
#________________________________________________________
methods_call <- function (cmd) {

  if (is.null(cmd)) return (NULL)

  cmd  <- paste(as.character(cmd), collapse = "\n")
  part <- NULL

  if (grepl("\\$(stats|taxa_stats)\\s*$", cmd)) {
    part <- sub("^.*\\$(stats|taxa_stats)\\s*$", "\\1", cmd)
    cmd  <- sub("\\$(stats|taxa_stats)\\s*$", "", cmd)
  }

  expr <- tryCatch(str2lang(cmd), error = function (e) NULL)

  call <- methods_call_expr(expr)
  if (!is.null(call)) call$part <- part

  return (call)
}

methods_call_expr <- function (expr) {

  if (!is.call(expr) || !is.symbol(expr[[1]])) return (NULL)

  fn  <- as.character(expr[[1]])
  fun <- get0(fn, envir = asNamespace('rbiom'), mode = 'function', inherits = FALSE)

  if (is.function(fun) && !is.primitive(fun))
    expr <- tryCatch(match.call(fun, expr), error = function (e) expr)

  list(
    fn   = fn,
    fun  = if (is.function(fun) && !is.primitive(fun)) fun,
    args = lapply(as.list(expr)[-1], methods_literal) )
}


#________________________________________________________
# The value of an expression made only of constants and c() / list()
# / unary minus / `:`. Anything else is returned unevaluated: recorded
# commands are never run.
#________________________________________________________
methods_literal <- function (e) {

  if (!is.language(e)) return (e)
  if (!is.call(e) || !is.symbol(e[[1]])) return (e)

  head <- as.character(e[[1]])
  if (!head %in% c('c', 'list', '-', ':')) return (e)

  parts <- lapply(as.list(e)[-1], methods_literal)
  if (any(vapply(parts, is.language, logical(1)))) return (e)

  tryCatch(do.call(head, parts, envir = baseenv()), error = function (err) e)
}


#________________________________________________________
# An argument of a parsed call: the value given, else the function's
# default (evaluated with the other arguments in scope, so that
# `stat.by = x` resolves), else `default`. Values that can't be seen
# (objects) give `default` too.
#________________________________________________________
methods_arg <- function (call, name, default = NULL) {

  if (is.null(call)) return (default)

  if (hasName(call$args, name)) {
    v <- call$args[[name]]

  } else if (!is.null(call$fun)) {
    v <- default_value(call$fun, name, call$args)
    if (inherits(v, 'rbiom_no_default')) return (default)

  } else {
    return (default)
  }

  if (is.language(v) || is.function(v)) return (default)

  return (v)
}


#________________________________________________________
# Parse the history's code lines into calls. Each entry is
# list(fn, fun, args), where fn is the function called, ".dataset" for
# a bundled dataset, ".assign" for `biom$field <- value`, or
# ".modified" for an unrecorded change.
#________________________________________________________
methods_history <- function (lines) {

  lines <- as.character(lines)
  lines <- lines[nzchar(trimws(lines))]
  out   <- list()

  for (line in lines) {

    if (grepl("^# .* was subsequently modified", line)) {
      out <- c(out, list(list(fn = '.modified', args = list())))
      next
    }

    code <- gsub("(^|\n)# ?", "\\1", line)   # entries that can't be replayed are comments
    code <- if (startsWith(line, "#")) code else line
    expr <- tryCatch(str2lang(code), error = function (e) NULL)

    if (!is.call(expr) || !identical(expr[[1]], as.name('<-')) || length(expr) != 3) next

    lhs <- expr[[2]]
    rhs <- expr[[3]]

    if (is.call(lhs) && identical(lhs[[1]], as.name('$'))) {
      entry <- list(fn = '.assign', args = list(field = as.character(lhs[[3]])))

    } else if (is.symbol(rhs)) {
      entry <- list(fn = '.dataset', args = list(name = as.character(rhs)))

    } else {
      entry <- methods_call_expr(rhs)
    }

    if (!is.null(entry)) out <- c(out, list(entry))
  }

  return (out)
}



#________________________________________________________
# Small text helpers.
#________________________________________________________

# "a", "a and b", "a, b, and c".
methods_join <- function (x, word = "and") {
  x <- as.character(x)
  n <- length(x)
  if (n <= 1) return (paste(x, collapse = ""))
  if (n == 2) return (paste(x, collapse = paste0(" ", word, " ")))
  paste0(paste(x[-n], collapse = ", "), ", ", word, " ", x[[n]])
}

methods_int <- function (x) formatC(round(as.numeric(x)), format = "d", big.mark = ",")
methods_num <- function (x) format(as.numeric(x), trim = TRUE)
methods_pct <- function (x) paste0(format(as.numeric(x) * 100, trim = TRUE), "%")

# A field name as recorded, without `==` / `!=` prefixes.
methods_field <- function (x) {
  if (!is.character(x)) return (NULL)
  x <- unique(sub("^(==|!=)", "", x))
  x[nzchar(x)]
}

# Expressions for display, without the backticks deparse() adds.
methods_expr <- function (e) {
  txt <- gsub("`", "", paste(deparse(e, width.cutoff = 500L), collapse = " "), fixed = TRUE)
  gsub("\\s+", " ", txt)
}

methods_cap <- function (s) paste0(toupper(substr(s, 1, 1)), substring(s, 2))



#________________________________________________________
# Data preparation, from the history.
#________________________________________________________
methods_prep <- function (info, ctx) {

  out    <- character()
  stated <- FALSE   # has a step said what the abundance table holds?

  for (h in info$history) {
    out <- c(out, methods_prep_step(h, ctx))
    if (h$fn %in% c('rarefy', 'biom_relativize', 'biom_rescale', 'biom_inflate')) stated <- TRUE
    if (h$fn %in% c('.modified') || (h$fn == '.assign' && identical(h$args$field, 'counts'))) stated <- FALSE
  }

  b <- info$biom

  if (!is.null(b)) {

    if (is.numeric(b$samples) && is.numeric(b$otus))
      out <- c(out, sprintf(
        "The %s comprised %s sample%s and %s OTU%s.",
        if (isTRUE(info$is_biom)) "dataset" else "analyzed dataset",
        methods_int(b$samples), if (b$samples == 1) "" else "s",
        methods_int(b$otus),    if (b$otus    == 1) "" else "s" ))

    if (!stated) out <- c(out, methods_table_sentence(b$table))
  }

  return (out)
}


# input_table_caption()'s sentence, as a statement about the input.
methods_table_sentence <- function (s) {

  if (!is_scalar_character(s) || is.na(s) || s == "No samples.") return (NULL)

  if (startsWith(s, "Rarefied"))
    return (sub("^Rarefied", "The abundance table was rarefied", s))

  paste0("The abundance table held ", tolower(substr(s, 1, 1)), substring(s, 2))
}


methods_prep_step <- function (h, ctx) {

  A <- function (name, default = NULL) methods_arg(h, name, default)

  switch(
    EXPR = h$fn,

    read_biom = {
      src <- A('src') %||% (if (length(h$args) > 0 && is.character(h$args[[1]])) h$args[[1]])
      if (is_scalar_character(src) && !is.na(src) && !grepl("[{\n\t]", src) && nchar(src) < 300) {
        sprintf("Data were imported from %s.", basename(src))
      } else {
        "Data were imported with read_biom()."
      }
    },

    as_rbiom = "The dataset was assembled from R objects with as_rbiom().",

    .dataset = sprintf("The %s example dataset distributed with rbiom was used.", h$args$name),

    biom_merge = {
      n <- sum(vapply(h$args, function (a) is.symbol(a) && grepl("^biom[0-9]+$", as.character(a)), logical(1)))
      if (n > 1) sprintf("The %i datasets were merged.", n) else "Datasets were merged."
    },

    rarefy = {
      ctx$use('ecodive', "rarefaction")
      depth <- A('depth')
      seed  <- A('seed', 0L)
      sprintf(
        "%s were rarefied %s%s; samples with fewer reads were excluded.",
        if (isTRUE(A('inflate', FALSE))) "Abundances were scaled to integer counts and samples" else "Samples",
        if (is.numeric(depth) && length(depth) == 1) sprintf("to %s reads per sample", methods_int(depth)) else "to an even depth",
        if (is.numeric(seed) && length(seed) == 1) sprintf(" (random seed %s)", methods_num(seed)) else "" )
    },

    biom_relativize = "Counts were converted to relative abundances.",

    biom_inflate = "Abundances were scaled to integer counts with biom_inflate().",

    biom_rescale = {
      rng <- A('range', c(0, 1))
      if (is.numeric(rng) && length(rng) == 2) {
        sprintf("Abundances were rescaled within each sample to the range %s to %s.", methods_num(rng[[1]]), methods_num(rng[[2]]))
      } else {
        "Abundances were rescaled within each sample."
      }
    },

    subset = ,
    subset_taxa = {
      what <- if (h$fn == 'subset') "Samples" else "OTUs"
      txt  <- if (length(h$args) >= 2) methods_expr(h$args[[2]]) else ""
      if (nzchar(txt) && nchar(txt) <= 150 && !grepl("<", txt, fixed = TRUE)) {
        sprintf("%s were restricted to those satisfying %s.", what, txt)
      } else {
        sprintf("%s were filtered on %s.", what, if (h$fn == 'subset') "their metadata" else "their taxonomy")
      }
    },

    `[` = "A subset of samples and/or OTUs was selected.",

    na.omit = {
      fields <- A('fields')
      if (is.character(fields) && !identical(fields, ".all")) {
        sprintf("Samples with missing values for %s were excluded.", methods_join(fields, "or"))
      } else {
        "Samples with missing metadata values were excluded."
      }
    },

    slice_sample = "A random subset of samples was drawn.",

    slice      = ,
    slice_head = ,
    slice_tail = ,
    slice_min  = ,
    slice_max  = sprintf("A subset of samples was selected with %s().", h$fn),

    # Edits to the metadata don't change the analysis.
    mutate = ,
    rename = ,
    within = NULL,

    .assign = switch(
      EXPR     = h$args$field %||% "",
      counts   = "The abundance table was replaced by direct assignment, which was not recorded as code.",
      tree     = "A phylogenetic tree was assigned to the dataset.",
      taxonomy = "The taxonomy was assigned directly.",
      NULL ),

    .modified = "The dataset was subsequently modified in a way that was not recorded.",

    sprintf("The dataset was processed with %s().", h$fn) )
}



#________________________________________________________
# The analysis and its statistics.
#________________________________________________________
methods_analysis <- function (info, ctx) {

  if (isTRUE(info$is_biom) || is.null(info$call)) {
    # A statistics table without a recorded command still names its test.
    return (if (!is.null(info$stats)) methods_stats(info$stats, ctx))
  }

  call <- info$call
  fn   <- call$fn
  A    <- function (name, default = NULL) methods_arg(call, name, default)


  #________________________________________________________
  # A statistics table: what was measured, then the test.
  #________________________________________________________
  if (isTRUE(info$is_stats)) {

    src <- if (fn %in% c('stats_table', 'distmat_stats')) info$data_call else call

    return (c(
      methods_measure(src, info, ctx),
      if (!is.null(src)) methods_pairs(src),
      methods_stats(info$stats, ctx, resp = methods_resp(src)) ))
  }


  out <- methods_measure(call, info, ctx)


  #________________________________________________________
  # What each function adds to the measurement.
  #________________________________________________________
  extra <- switch(
    EXPR = fn,

    bdiv_table    = ,
    bdiv_boxplot  = ,
    bdiv_corrplot = methods_pairs(call),

    bdiv_heatmap = {
      clust <- A('clust', "complete")
      if (is.character(clust) && length(clust) >= 1 && !anyNA(clust)) {
        sprintf("Samples were ordered by '%s' hierarchical clustering of the plotted distances.", clust[[1]])
      } else if (!isFALSE(clust) && !is_na(clust)) {
        "Samples were ordered by hierarchical clustering."
      }
    },

    bdiv_ord_plot     = ,
    bdiv_ord_table    = ,
    distmat_ord_table = methods_ord(call, info, ctx),

    rare_corrplot  = ,
    rare_multiplot = NULL,   # in methods_measure()

    rare_stacked = {
      rline <- A('rline', TRUE)
      if (is.numeric(rline) && length(rline) == 1) {
        sprintf("Sequencing depth was summarized per sample, showing the reads that rarefying to %s reads per sample would keep and discard.", methods_int(rline))
      } else {
        "Sequencing depth was summarized per sample, showing the reads that rarefaction would keep and discard."
      }
    },

    taxa_heatmap = c(
      switch(
        EXPR = as.character(A('rescale', "none"))[[1]],
        rows = "Abundances were rescaled within each taxon for display.",
        cols = "Abundances were rescaled within each sample for display." ),
      methods_clust(A('clust', "complete"), A('dist', "euclidean"), "Taxa and samples") ),

    taxa_stacked = {
      order.by <- methods_field(A('order.by'))
      clust    <- methods_clust(A('clust', "complete"), A('dist', "euclidean"), "Samples")
      if (length(order.by) > 0) {
        sprintf("Samples were arranged by %s.", methods_join(order.by))
      } else {
        clust
      }
    },

    taxa_top = {
      by    <- methods_field(A('by'))
      n     <- A('n', 5)
      n_txt <- if (is.numeric(n) && length(n) == 1) sprintf("%s most abundant", methods_int(n)) else "most abundant"
      if (length(by) == 0) {
        sprintf("Taxa were ranked by their mean relative abundance across all samples, and the %s were listed.", n_txt)
      } else {
        sprintf(
          "Taxa were ranked by their mean relative abundance within each %s, and the %s were listed for each.",
          if ('.sample' %in% by) "sample" else sprintf("group of %s", methods_join(by)), n_txt )
      }
    },

    taxa_shared = methods_shared(call),

    taxa_upset = c(
      methods_shared(call),
      sprintf("The number of taxa in each combination of groups was displayed as an UpSet plot (%s).", ctx$cite('lex2014')) ),

    NULL )

  out <- c(out, extra)


  #________________________________________________________
  # Trend lines on scatter plots that weren't tested.
  #________________________________________________________
  if (grepl("_corrplot$|^rare_multiplot$", fn) && is.null(info$stats)) {
    fit <- A('fit', "gam")
    x   <- if (startsWith(fn, "rare_")) "rarefaction depth" else methods_field(A('x'))
    if (is_scalar_character(fit) && length(x) == 1)
      out <- c(out, sprintf(
        "Trend lines against %s were fitted with %s.", x,
        switch(
          EXPR = fit,
          lm   = "linear models",
          log  = sprintf("linear models on log(%s)", x),
          gam  = { ctx$use('mgcv', "generalized additive models")
                   "generalized additive models (GAM)" },
          fit ) ))
  }


  #________________________________________________________
  # Statistics attached to the result.
  #________________________________________________________
  out <- c(
    out,
    methods_stats(info$stats,      ctx, resp = methods_resp(call)),
    methods_stats(info$taxa_stats, ctx) )

  if (fn == 'taxa_forest' && !is.null(info$stats))
    out <- c(out, sprintf(
      "Taxa with an adjusted p-value of at most %s were considered significant.",
      methods_num(A('p.sig', 0.05)) ))

  return (out)
}


# What a function's response is called in the statistics sentences.
methods_resp <- function (call) {
  fn <- call$fn %||% ""
  if (grepl("^(adiv_|rare_)", fn)) return ("alpha diversity")
  if (grepl("^bdiv_", fn))         return ("beta diversity")
  if (grepl("^taxa_", fn))         return ("taxon abundance")
  "the response"
}



#________________________________________________________
# What was measured: diversity metrics or taxa abundances.
#________________________________________________________
methods_measure <- function (call, info, ctx) {

  fn <- call$fn %||% ""

  if (grepl("^rare_(corrplot|multiplot)$", fn)) return (methods_rare(call, info, ctx))
  if (grepl("^adiv_", fn))                      return (methods_adiv(call, ctx))
  if (grepl("^bdiv_", fn))                      return (methods_bdiv(call, ctx))

  if (fn %in% c('taxa_table', 'taxa_matrix', 'taxa_boxplot', 'taxa_corrplot', 'taxa_stats', 'taxa_heatmap', 'taxa_stacked', 'taxa_top', 'taxa_shared', 'taxa_upset', 'taxa_sums', 'taxa_means'))
    return (methods_taxa(call, info, ctx))

  NULL
}


# The names of metrics, each with its reference when there is one.
methods_metrics <- function (ids, div, ctx, alpha = NULL) {

  ctx$use('ecodive', "diversity metrics")

  vapply(as.character(ids), FUN.VALUE = '', USE.NAMES = FALSE, function (x) {

    m <- tryCatch(ecodive::match_metric(x, div = div), error = function (e) NULL)
    if (is.null(m)) return (x)

    name <- methods_metric_name(m$name)
    key  <- if (is.null(METRIC_REFS[[m$id]])) NA else paste0("metric_", m$id)
    note <- c(
      if (identical(m$id, 'generalized_unifrac') && is.numeric(alpha)) sprintf("alpha = %s", methods_num(alpha)),
      if (!is.na(key)) ctx$cite(unname(key)) )

    if (length(note) == 0) name else sprintf("%s (%s)", name, paste(note, collapse = "; "))
  })
}


# ecodive's display name, as it reads in a sentence: "the Shannon
# diversity index", "unweighted UniFrac", "Faith's phylogenetic diversity".
methods_metric_name <- function (name) {

  name <- sub("\\s*\\([^)]*\\)\\s*$", "", name)   # "(GUniFrac)", "(ACE)", ...

  generic <- c(
    "Diversity", "Index", "Distance", "Dissimilarity", "Divergence", "Richness",
    "Estimator", "Phylogenetic", "Observed", "Features", "Coverage", "Alpha",
    "Weighted", "Unweighted", "Normalized", "Generalized", "Variance-Adjusted",
    "Abundance-based", "Inverse", "Squared", "Robust", "Probabilistic", "Symmetric" )

  for (w in generic)
    name <- gsub(sprintf("\\b%s\\b", w), tolower(w), name)

  if (grepl("(index|distance|dissimilarity|divergence|estimator)$", name) && !grepl("'s ", name, fixed = TRUE))
    name <- paste("the", name)

  return (name)
}


methods_transform <- function (transform) {
  if (!is_scalar_character(transform)) return (NULL)
  switch(
    EXPR    = transform,
    rank    = "rank-transformed",
    log     = "log-transformed",
    log1p   = "log(1 + x)-transformed",
    sqrt    = "square-root-transformed",
    NULL )
}


methods_adiv <- function (call, ctx) {

  adiv <- methods_arg(call, 'adiv', "Shannon")
  if (!is.character(adiv)) return ("Alpha diversity was computed.")

  tf <- methods_transform(methods_arg(call, 'transform', "none"))

  sprintf(
    "Alpha diversity was measured with %s%s.",
    methods_join(methods_metrics(adiv, 'alpha', ctx)),
    if (is.null(tf)) "" else sprintf("; values were %s", tf) )
}


methods_bdiv <- function (call, ctx) {

  bdiv <- methods_arg(call, 'bdiv', "bray")
  if (!is.character(bdiv)) return ("Beta diversity between samples was computed.")

  tf <- methods_transform(methods_arg(call, 'transform', "none"))

  sprintf(
    "Beta diversity between samples was measured with %s%s.",
    methods_join(methods_metrics(bdiv, 'beta', ctx, alpha = methods_arg(call, 'alpha', 0.5))),
    if (is.null(tf)) "" else sprintf("; distances were %s", tf) )
}


# Which sample pairs were compared (within / between, `==` / `!=`).
methods_pairs <- function (call) {

  within  <- methods_arg(call, 'within')
  between <- methods_arg(call, 'between')

  for (i in c('x', 'stat.by', 'facet.by', 'md')) {
    v <- methods_arg(call, i)
    if (!is.character(v)) next
    within  <- c(within,  sub("^==", "", v[startsWith(v, "==")]))
    between <- c(between, sub("^!=", "", v[startsWith(v, "!=")]))
  }

  s <- within_between_caption(
    within  = if (is.character(within))  within,
    between = if (is.character(between)) between )

  if (is.null(s)) return (NULL)

  sub("^Comparing samples (.*)\\.$", "Only pairs of samples \\1 were compared.", s)
}


methods_rare <- function (call, info, ctx) {

  ctx$use('ecodive', "rarefaction")

  depths <- tryCatch(info$data[['.depth']], error = function (e) NULL)

  # rare_multiplot(): the curves are one of the patchwork's panels.
  if (is.null(depths))
    for (panel in tryCatch(info$obj[['patches']][['plots']], error = function (e) NULL)) {
      d <- tryCatch(panel[['data']][['.depth']], error = function (e) NULL)
      if (is.numeric(d)) { depths <- d; break }
    }

  depths <- depths[is.finite(depths)]

  first <- if (is.numeric(depths) && length(depths) > 0) {
    sub("^Rarefied", "Each sample was rarefied", input_table_caption(depth = depths))
  } else {
    "Each sample was rarefied to a series of depths."
  }

  c(first, sub(
    pattern     = "\\.$",
    replacement = " at each depth.",
    x           = methods_adiv(call, ctx) ))
}


# A phrase for the taxonomic rank(s): "the Genus rank", "OTUs", ...
methods_rank <- function (call, info) {

  rank <- methods_arg(call, 'rank', -1)

  if (!is.character(rank)) {

    # The result's `.rank` column names the rank a number selected.
    for (df in list(info$obj, info$data)) {
      r <- if (is.data.frame(df) && hasName(df, '.rank')) unique(as.character(df[['.rank']]))
      if (length(r) > 0) { rank <- r; break }
    }
  }

  if (is.character(rank)) {
    rank <- unique(rank)
    otu  <- tolower(rank) %in% c('.otu', 'otu')
    return (methods_join(c(
      if (any(!otu)) sprintf("the %s rank%s", methods_join(rank[!otu]), if (sum(!otu) > 1) "s" else ""),
      if (any(otu))  "the OTU level" )))
  }

  if (is.numeric(rank) && length(rank) == 1) {
    if (rank == -1) return ("the lowest taxonomic rank")
    if (rank ==  0) return ("the OTU level")
    if (rank ==  1) return ("the highest taxonomic rank")
  }

  "the selected taxonomic rank"
}


methods_taxa <- function (call, info, ctx) {

  fn   <- call$fn
  rank <- methods_rank(call, info)

  if (fn %in% c('taxa_shared', 'taxa_upset'))
    return (sprintf("Counts were summarized at %s.", rank))

  transform <- methods_arg(call, 'transform', "none")
  if (fn %in% c('taxa_stacked', 'taxa_heatmap')) transform <- "percent"

  units <- if (identical(transform, "percent")) "expressed as relative abundances" else methods_transform(transform)

  out <- sprintf(
    "Counts were summarized at %s%s.", rank,
    if (is.null(units)) "" else paste(" and", units) )

  if (fn == 'taxa_top') return (out)

  taxa  <- methods_arg(call, 'taxa')
  other <- methods_arg(call, 'other', FALSE)
  other <- isTRUE(other) || is_scalar_character(other)

  sel <- if (is.numeric(taxa) && length(taxa) == 1 && taxa >= 1) {
    sprintf("The %s most abundant taxa were retained", methods_int(taxa))
  } else if (is.numeric(taxa) && length(taxa) == 1) {
    sprintf("Taxa with a mean relative abundance of at least %s were retained", methods_pct(taxa))
  } else if (is.character(taxa)) {
    sprintf("%s named tax%s retained", methods_cap(methods_int(length(taxa))), if (length(taxa) == 1) "on was" else "a were")
  }

  if (!is.null(sel))
    out <- c(out, paste0(sel, if (other) "; the remaining taxa were pooled." else "."))

  return (out)
}


methods_shared <- function (call) {

  stat.by <- methods_field(methods_arg(call, 'stat.by'))
  detect  <- methods_arg(call, 'detect', 0)
  prev    <- methods_arg(call, 'prev',   0)
  group   <- if (length(stat.by) == 1) sprintf("a %s group", stat.by) else "the dataset"

  sprintf(
    "A taxon was counted as present in a sample when %s, and as belonging to %s when it was present in %s.",
    if (is.numeric(detect) && detect > 0) sprintf("its relative abundance exceeded %s", methods_pct(detect)) else "it was observed",
    group,
    if (is.numeric(prev) && prev > 0) sprintf("at least %s of its samples", methods_pct(prev)) else "at least one of its samples" )
}


# Ordering by hierarchical clustering, in clust_caption()'s words.
methods_clust <- function (clust, dist, what) {

  if (isFALSE(clust) || is_na(clust)) return (NULL)

  if (!is.character(clust) || anyNA(clust) || !is.character(dist))
    return (sprintf("%s were ordered by hierarchical clustering.", what))

  parts <- unique(mapply(clust_caption, clust, dist, USE.NAMES = FALSE))
  parts <- sub("clustering", "hierarchical clustering", parts, fixed = TRUE)

  sprintf("%s were ordered by %s.", what, methods_join(parts))
}


methods_ord <- function (call, info, ctx) {

  axes <- info$axes

  ord <- if (is.data.frame(axes) && hasName(axes, '.ord')) unique(as.character(axes[['.ord']])) else methods_arg(call, 'ord', "PCoA")
  if (!is.character(ord)) return ("Samples were ordinated.")
  if (identical(ord, ".all")) ord <- CHOICES$ord

  phrases <- vapply(ord, FUN.VALUE = '', USE.NAMES = FALSE, function (o) {
    switch(
      EXPR = o,

      PCoA = sprintf("principal coordinates analysis (PCoA; %s)", ctx$cite('gower1966')),

      NMDS = {
        ctx$use('vegan', "NMDS")
        stress <- if (is.data.frame(axes) && hasName(axes, '.stress')) axes[['.stress']][as.character(axes[['.ord']]) == "NMDS"]
        stress <- unique(stress[is.finite(stress)])
        fmt    <- function (x) formatC(x, format = "f", digits = 3)
        sprintf(
          "non-metric multidimensional scaling (NMDS; %s%s)", ctx$cite('kruskal1964'),
          if (length(stress) == 1)     sprintf("; stress = %s", fmt(stress))
          else if (length(stress) > 1) sprintf("; stress from %s to %s", fmt(min(stress)), fmt(max(stress)))
          else "" )
      },

      tSNE = { ctx$use('Rtsne', "t-SNE")
               sprintf("t-distributed stochastic neighbor embedding (t-SNE; %s)", ctx$cite('maaten2008')) },

      UMAP = { ctx$use('uwot', "UMAP")
               sprintf("uniform manifold approximation and projection (UMAP; %s)", ctx$cite('mcinnes2018')) },

      o )
  })

  k    <- methods_arg(call, 'k', 2)
  seed <- methods_arg(call, 'seed', 0)
  rng  <- any(ord != "PCoA") && is.numeric(seed) && length(seed) == 1

  sprintf(
    "Samples were ordinated in %s dimensions by %s%s.",
    if (is.numeric(k) && length(k) == 1) methods_int(k) else "two",
    methods_join(phrases),
    if (rng) sprintf(", with random seed %s", methods_num(seed)) else "" )
}



#________________________________________________________
# Sentences for a statistics table, from its header (`tbl_sum`) and
# the arguments of the command that made it. `resp` names what was
# tested.
#________________________________________________________
methods_stats <- function (stats, ctx, resp = "the response") {

  if (!is.data.frame(stats)) return (NULL)

  ts <- attr(stats, 'tbl_sum', exact = TRUE)
  if (is.null(ts) || is.null(names(ts))) return (NULL)
  ts <- as.list(ts)

  call <- methods_call(attr(stats, 'cmd', exact = TRUE))
  A    <- function (name, default = NULL) methods_arg(call, name, default)

  out <- if (!is.null(ts[['Method']])) {  methods_stats_da(ts, A, ctx)
  } else if (!is.null(ts[['Test']]))   {  methods_stats_perm(stats, ts, A, ctx)
  } else if (!is.null(ts[['Model']]))  {  methods_stats_model(stats, ts, A, ctx, resp)
  } else                               {  return (NULL) }


  #________________________________________________________
  # Multiple comparisons. One test needs no adjustment.
  #________________________________________________________
  p.adj <- A('p.adj', "fdr")

  if (nrow(stats) > 1 && is_scalar_character(p.adj) && p.adj != "none") {

    label <- switch(
      EXPR       = p.adj,
      BH         = ,
      fdr        = sprintf("the Benjamini-Hochberg false discovery rate procedure (%s)", ctx$cite('benjamini1995')),
      BY         = sprintf("the Benjamini-Yekutieli false discovery rate procedure (%s)", ctx$cite('benjamini2001')),
      holm       = sprintf("Holm's method (%s)", ctx$cite('holm1979')),
      hochberg   = sprintf("Hochberg's method (%s)", ctx$cite('hochberg1988')),
      hommel     = sprintf("Hommel's method (%s)", ctx$cite('hommel1988')),
      bonferroni = "the Bonferroni correction",
      sprintf("the '%s' method", p.adj) )

    out <- c(out, sprintf(
      "P-values were adjusted for multiple comparisons with %s%s.", label,
      if (!is.null(ts[['Method']])) ", across taxa within each comparison" else "" ))
  }

  return (out)
}


# What the tests were run separately for.
methods_split <- function (stats, A) {

  split.by <- methods_field(A('split.by'))

  cols <- c('.adiv' = "alpha diversity metric", '.bdiv' = "beta diversity metric", '.taxa' = "taxon", '.rank' = "rank", '.ord' = "ordination")
  auto <- names(cols)[vapply(names(cols), function (i) hasName(stats, i) && length(unique(stats[[i]])) > 1, logical(1))]

  fields <- setdiff(split.by, names(cols))
  each   <- unique(c(
    cols[intersect(c(split.by, auto), names(cols))],
    if (length(fields) > 0) sprintf("level of %s", fields) ))
  if (length(each) == 0) return (NULL)

  sprintf("Tests were run separately for %s.", methods_join(paste("each", each)))
}


# Permutation tests on distance matrices.
methods_stats_perm <- function (stats, ts, A, ctx) {

  ctx$use('vegan', "permutation tests")

  txt    <- sub("\\.\\s*$", "", gsub("`", "", ts[['Test']], fixed = TRUE))
  label  <- sub(" ~ .*$", "", txt)
  rest   <- sub("^.*? ~ ", "", txt, perl = TRUE)
  term   <- trimws(sub(";.*$", "", rest))
  scheme <- if (grepl(";", rest, fixed = TRUE)) trimws(sub("^[^;]*;", "", rest)) else ""

  if (term == "groups") term <- "the sample groups"

  seed  <- A('seed', 0)
  notes <- c(
    if (nzchar(scheme)) scheme,
    if (is.numeric(seed) && length(seed) == 1) sprintf("random seed %s", methods_num(seed)) )
  notes <- if (length(notes) > 0) sprintf(" (%s)", paste(notes, collapse = "; ")) else ""

  first <- if (startsWith(label, "PERMANOVA")) {
    sprintf(
      "The effect of %s on community composition was tested by permutational multivariate analysis of variance (PERMANOVA; %s) with adonis2 from vegan%s.",
      term, ctx$cite('anderson2001'), notes )

  } else if (startsWith(label, "MRPP")) {
    sprintf(
      "Differences in community composition by %s were tested with the multi-response permutation procedure (MRPP; %s) with mrpp from vegan%s.",
      term, ctx$cite('mielke2001'), notes )

  } else if (startsWith(label, "Dispersion")) {
    sprintf(
      "Homogeneity of multivariate dispersion by %s was tested (%s) with betadisper and permutest from vegan%s.",
      term, ctx$cite('anderson2006'), notes )

  } else if (startsWith(label, "Mantel")) {
    sprintf(
      "The correlation between beta diversity and differences in %s was tested with a Mantel test (%s) on Spearman's rank correlation, with mantel from vegan%s.",
      term, ctx$cite('mantel1967'), notes )

  } else if (grepl("envfit", label, fixed = TRUE)) {
    return (sprintf("The abundance of each taxon was fitted onto the ordination with envfit from vegan%s.", notes))

  } else {
    sprintf("Statistics: %s.", txt)
  }

  c(first, methods_split(stats, A))
}


# Rank tests and regression models (stats_table() and its wrappers).
methods_stats_model <- function (stats, ts, A, ctx, resp) {

  model   <- ts[['Model']]
  stat.by <- methods_field(A('stat.by'))
  regr    <- methods_field(A('regr'))
  strata  <- methods_field(A('strata'))
  test    <- A('test', "emmeans")
  fit     <- A('fit',  "gam")
  level   <- A('level', 0.95)
  groups  <- if (length(stat.by) == 1) sprintf("%s groups", stat.by) else "groups"
  Resp    <- methods_cap(resp)

  ci <- if (is.numeric(level) && length(level) == 1) sprintf("Confidence intervals were computed at the %s level.", methods_pct(level))


  if (startsWith(model, "wilcox.test(.resp ~")) {

    first <- sprintf(
      "%s was compared between each pair of %s with Wilcoxon rank-sum (Mann-Whitney U) tests (%s).",
      Resp, groups, ctx$cite('mann1947', 'wilcoxon1945') )

  } else if (startsWith(model, "wilcox.test")) {

    first <- sprintf(
      "%s was compared with %s using a one-sample Wilcoxon signed-rank test (%s).",
      Resp, methods_num(A('mu', 0)), ctx$cite('wilcoxon1945') )

  } else if (startsWith(model, "kruskal.test")) {

    ci    <- NULL
    first <- sprintf(
      "%s was compared across %s with the Kruskal-Wallis rank-sum test (%s).",
      Resp, groups, ctx$cite('kruskal1952') )

  } else if (grepl("^(lm|gam)\\(", model)) {

    has_regr <- length(regr) == 1
    has_by   <- length(stat.by) == 1
    smooth   <- has_regr && identical(fit, "gam")
    is_gam   <- startsWith(model, "gam(")

    if (is_gam) ctx$use('mgcv', if (smooth) "generalized additive models" else "mixed models")

    covs <- if (!is.null(ts[['Covariates']]))
      gsub("`", "", gsub("\\.cov[0-9]+ = ", "", strsplit(ts[['Covariates']], ", ", fixed = TRUE)[[1]]), fixed = TRUE)

    x <- if (has_regr && identical(fit, "log")) sprintf("log(%s)", regr) else regr

    re <- if (length(strata) == 1) sprintf("a random intercept for each %s", strata)

    kind <- if (smooth) {
      sprintf(
        "a generalized additive model (GAM; %s) with a smooth term for %s%s, fitted by REML",
        ctx$cite('wood2011'), regr, if (is.null(re)) "" else paste(" and", re) )
    } else if (!is.null(re)) {
      sprintf("a linear mixed model with %s, fitted by REML with mgcv (%s)", re, ctx$cite('wood2011'))
    } else {
      "a linear model"
    }

    terms <- c(
      if (has_regr && !smooth) x,
      if (has_by)              stat.by,
      covs )

    at <- A('at')
    at <- if (is.numeric(at)) sprintf("%s = %s", regr, methods_join(methods_num(at))) else sprintf("the median of %s", regr)

    emm <- TRUE

    quantity <- if (identical(test, "emtrends")) {
      slope <- sprintf("the slope of %s against %s%s", resp, x, if (smooth) sprintf(" at %s", at) else "")
      if (has_by) sprintf("pairwise differences between %s in %s were", groups, slope) else sprintf("%s was", slope)

    } else if (has_by && has_regr) {
      sprintf("pairwise differences between the means of %s at %s were", groups, at)

    } else if (has_by) {
      sprintf("pairwise differences between the means of %s were", groups)

    } else if (has_regr) {
      emm <- FALSE
      if (smooth) sprintf("the significance of the smooth term for %s was", regr) else sprintf("the slope of %s against %s was", resp, x)

    } else {
      sprintf("the difference between the mean of %s and %s was", resp, methods_num(A('mu', 0)))
    }

    how <- if (identical(test, "emtrends")) "estimated marginal trends" else "estimated marginal means"
    if (emm) ctx$use('emmeans', how)

    first <- c(
      sprintf(
        "%s was modeled%s with %s%s.",
        Resp,
        if (length(terms) > 0) sprintf(" as a function of %s", methods_join(terms)) else "",
        kind,
        if (has_by) sprintf("; a separate model was fitted to each pair of %s", groups) else "" ),
      sprintf(
        "%s tested%s.",
        methods_cap(quantity),
        if (emm) sprintf(" with %s from the emmeans package", how) else "" ))

  } else {
    first <- sprintf("Statistics were computed with %s.", model)
  }

  c(first, methods_split(stats, A), ci)
}


# Differential abundance.
methods_stats_da <- function (ts, A, ctx) {

  rank     <- A('rank')
  strata   <- methods_field(A('strata'))
  prev.cut <- A('prev.cut', 0.1)
  level    <- A('level', 0.95)
  model    <- sub(",? per taxon$", "", gsub("`", "", ts[['Model']] %||% "", fixed = TRUE))

  if (!startsWith(ts[['Method']], "LinDA"))
    return (sprintf("Differential abundance was tested with %s.", ts[['Method']]))

  if (length(strata) == 1) ctx$use('lmerTest', "linear mixed models")

  c(
    sprintf(
      "Differential abundance%s was tested with LinDA (%s) as implemented in rbiom: for each taxon, centered log-ratio transformed counts were fitted with a linear %smodel%s, and each coefficient was corrected by its mode across taxa.",
      if (is.character(rank)) sprintf(" at the %s rank%s", methods_join(rank), if (length(rank) > 1) "s" else "") else "",
      ctx$cite('zhou2022'),
      if (length(strata) == 1) "mixed " else "",
      if (nzchar(model)) sprintf(" (%s)", model) else "" ),

    if (is.numeric(prev.cut) && length(prev.cut) == 1 && prev.cut > 0)
      sprintf("Taxa observed in fewer than %s of the samples were not tested.", methods_pct(prev.cut)),

    if (is.numeric(level) && length(level) == 1)
      sprintf("Confidence intervals were computed at the %s level.", methods_pct(level)) )
}



#________________________________________________________
# "Analyses were performed with rbiom 4.0.0 in R 4.6.1 (R Core Team
# 2026), using ..." and
# the references for the software named.
#________________________________________________________
methods_software <- function (pkgs, ctx) {

  version <- function (pkg) {
    v <- tryCatch(utils::packageDescription(pkg, fields = "Version"), error = function (e) NA, warning = function (w) NA)
    if (is_scalar_character(v) && !is.na(v)) return (v)
    tryCatch(as.character(utils::packageVersion(pkg)), error = function (e) NULL)
  }
  named   <- function (pkg) paste(c(pkg, version(pkg)), collapse = " ")

  ref <- function (pkg) {
    def <- METHODS_PKGS[[pkg]]
    if (is.null(def) || !is.null(def$key)) return (NULL)
    title <- tryCatch(utils::packageDescription(pkg)$Title, error = function (e) NULL, warning = function (w) NULL)
    if (!is_scalar_character(title) || is.na(title)) title <- def$title
    sprintf(
      "%s. %s: %s. R package%s.", def$authors, pkg, gsub("\\s+", " ", title),
      if (is.null(v <- version(pkg))) "" else paste(" version", v) )
  }

  # Most specific first; graphics last.
  order <- c('ecodive', 'vegan', 'Rtsne', 'uwot', 'mgcv', 'emmeans', 'lmerTest', 'ggplot2')
  used  <- c(intersect(order, names(pkgs)), setdiff(names(pkgs), order))

  parts <- vapply(used, FUN.VALUE = '', function (pkg) {
    key <- METHODS_PKGS[[pkg]]$key
    sprintf(
      "%s%s for %s", named(pkg),
      if (is.null(key)) "" else sprintf(" (%s)", ctx$cite(key)),
      methods_join(pkgs[[pkg]]) )
  })

  r_ver  <- paste(R.version$major, R.version$minor, sep = ".")
  r_cite <- sprintf("R Core Team %s", R.version$year)

  text <- sprintf(
    "Analyses were performed with %s in R %s (%s)%s.", named('rbiom'), r_ver, r_cite,
    if (length(parts) > 0) paste(", using", methods_join(parts)) else "" )

  refs <- c(
    unlist(lapply(c('rbiom', used), ref)),
    sprintf("R Core Team (%s). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria.", R.version$year) )

  list(text = text, refs = refs)
}
