# Describe the arguments of rbiom functions.

A machine-readable description of what each argument of an rbiom
function accepts: its kind, its fixed choices, numeric limits, whether
`NULL` is allowed and what it means, and where dataset-dependent choices
come from. Another program (such as a web interface) can use it to draw
input controls and to check a request without calling the function. The
same description is what the functions use to check their own arguments,
so the two cannot disagree.

## Usage

``` r
fn_schema(fn = NULL, json = FALSE)
```

## Arguments

- fn:

  The name of a function, e.g. `"adiv_boxplot"`. `NULL` describes every
  covered function. Default: `NULL`

- json:

  Return the description as JSON text instead of an R list. Default:
  `FALSE`

## Value

With `fn`, one entry (a list). Without, a named list of entries, one per
covered function, so `names(fn_schema())` lists the covered functions.
With `json = TRUE`, the same as a JSON string.

## Entries

Each entry has:

- `fn`:

  The function's name.

- `returns`:

  What the function returns: `type` is one of `"plot"`, `"table"`,
  `"factor"`, `"matrix"`, `"dist"`, `"biom"`, or `"file"`; `has` lists
  the `$` fields the result carries (e.g. `"data"`, `"code"`,
  `"stats"`).

- `rules`:

  (Optional) Constraints between arguments, each one of: `one_of`
  (exactly one of these arguments must be non-`NULL`); `any_of` (at
  least one must be); or `when` (conditions, as for an argument's
  `when`) with `requires` (these arguments must then be non-`NULL`) or
  `forbids` (they must then be `NULL`).

- `params`:

  One entry per argument, in the order of
  [`formals()`](https://rdrr.io/r/base/formals.html).

## Arguments

Every argument has `name`, `type`, and `desc` (the first sentence of its
help text). The other fields are present only when they apply.

- `type`:

  One of: `"logical"`, `"integer"`, `"number"`, `"string"`; `"choice"`
  (one or more of `choices`); `"layers"` (plot layers: `choices` lists
  each `letter` and `name`; letters may be run together, as in `"xd"`);
  `"field"` (metadata field names: see `class`); `"rank"` (taxonomic
  ranks of the dataset, by name or by position: `1` is the highest rank,
  `-1` the lowest, `0` the OTUs); `"taxa"` (see `forms`); `"metric"`
  (`div` is `"alpha"` or `"beta"`; `choices` lists each metric's `id`,
  `name`, whether it is `phylo`genetic, needing a tree, and whether it
  needs integer counts, `int_only`); `"aesthetic"` (colors, shapes, or
  patterns: `aes`, and `by`, the argument whose values they are mapped
  to); `"palette"`; `"tracks"`; `"biom"`, `"tree"`, `"dist"`,
  `"vector"`, `"file"`, `"dots"` (R objects, a file path, and `...`: see
  `r_only`).

- `required`:

  `TRUE` when the argument has no default.

- `default`:

  The default value, from
  [`formals()`](https://rdrr.io/r/base/formals.html).

- `default_expr`, `default_note`:

  For a default that is an expression, such as `stat.by = x`: the
  expression (not evaluated) and what it means. `default` is then
  `NULL`.

- `null`:

  Present when `NULL` is accepted: what it means.

- `choices`:

  The fixed options. Matching ignores case and accepts unambiguous
  abbreviations.

- `choices_when`:

  Options that depend on other arguments: a list of `when` conditions
  (as below) with the `choices` that apply while they hold. The first
  that holds applies; otherwise `choices`.

- `min`, `max`:

  Numeric limits. `"integer"` values must be whole.

- `max_n`:

  How many values may be given. Absent: any number.

- `per_side`:

  One value for both rows and columns of a heat map, or two: for rows,
  then columns.

- `class`:

  For fields: `"categorical"`, `"numeric"`, or `"any"`.

- `coerce`:

  For categorical fields: numeric fields are accepted too, and treated
  as categories.

- `ops`:

  For fields: field names may be prefixed with `"=="` (compare within
  groups) or `"!="` (between groups).

- `extra`:

  For fields: values accepted besides the dataset's fields, such as
  `".taxa"`.

- `forms`:

  The forms accepted, when there are several: `"logical"` (`TRUE` or
  `FALSE`); `"true"` (`TRUE` only); `"integer"`; `"string"`; `"choice"`
  (one of `choices`); `"palette"` (a palette name from
  [`list_palettes()`](https://cmmr.github.io/rbiom/reference/list_palettes.md),
  optionally prefixed with `"-"` to reverse it); `"values"` (a vector of
  colors, shapes, or patterns); `"named_values"` (a vector named by the
  values of `by`); `"count"` (the n most abundant taxa); `"fraction"` (a
  number below 1: taxa with at least that mean relative abundance);
  `"names"` (taxa names at the chosen rank); `"fields"` (metadata field
  names); `"spec"` (a list of options; see the function's help).

- `palettes`:

  The palette kinds suited to the argument.

- `needs`:

  An optional package the argument requires.

- `when`:

  Conditions under which the argument has any effect; all must hold.
  Each names a `param` and one of: `is` (its value is one of these),
  `is_not`, `has_any` (it includes one of these), `is_null`, `not_null`.
  Values are canonical: layer names, not letters, and metric ids.
  Absent: the argument may always matter.

- `deprecated`:

  The version that deprecated the argument.

- `r_only`:

  Only an R caller can set it: an R object, a file path, a host setting
  such as `cpus`, or `...`.
  [`fn_check()`](https://cmmr.github.io/rbiom/reference/fn_check.md)
  refuses it.

- `note`:

  Anything else a caller needs to know.

In JSON, `NULL` is `null`, numbers are written in full, `Inf` is the
string `"Inf"`, and the list fields (`choices`, `forms`, `ops`,
`palettes`, `has`, and the values in `when`) are always arrays.

## See also

[`fn_check()`](https://cmmr.github.io/rbiom/reference/fn_check.md) to
check arguments against a dataset;
[`list_palettes()`](https://cmmr.github.io/rbiom/reference/list_palettes.md)
and
[`list_metrics()`](https://cmmr.github.io/ecodive/reference/metrics.html).

## Examples

``` r
    library(rbiom)

    head(names(fn_schema()))
#> [1] "adiv_boxplot"  "adiv_corrplot" "adiv_stats"    "adiv_matrix"  
#> [5] "adiv_table"    "bdiv_boxplot" 

    x <- fn_schema("adiv_boxplot")
    x$returns
#> $type
#> [1] "plot"
#> 
#> $has
#> [1] "data"  "code"  "stats"
#> 
    str(x$params[[10]])
#> List of 6
#>  $ name   : chr "test"
#>  $ type   : chr "choice"
#>  $ choices: chr [1:5] "auto" "none" "wilcox" "emmeans" ...
#>  $ max_n  : num 1
#>  $ default: chr "auto"
#>  $ desc   : chr "Method for computing p-values."

    fn_schema("rarefy", json = TRUE)
#> {"fn":"rarefy","returns":{"type":"biom","has":["code"]},"params":[{"name":"biom","type":"biom","required":true,"r_only":true,"desc":"An rbiom object, or any value accepted by `as_rbiom()`."},{"name":"depth","type":"integer","min":1,"max_n":1,"null":"Chosen automatically, as by suggest_rarefy_depth().","default":null,"desc":"The number of observations to keep per sample."},{"name":"seed","type":"integer","min":0,"max_n":1,"default":0,"desc":"Random seed for selecting which observations to keep."},{"name":"inflate","type":"logical","max_n":1,"default":false,"desc":"Logical. Handling for non-integer data (e.g. relative abundances)."},{"name":"clone","type":"logical","r_only":true,"max_n":1,"default":true,"desc":"Create a copy of `biom` before modifying."},{"name":"cpus","type":"integer","r_only":true,"min":1,"max":256,"max_n":1,"null":"All available cores.","default":null,"default_expr":"n_cpus()","default_note":"All available cores.","desc":"The number of CPUs to use."}]} 
```
