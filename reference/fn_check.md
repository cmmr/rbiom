# Check arguments for an rbiom function against a dataset.

Applies the checks the function itself would apply on entry, using its
[`fn_schema()`](https://cmmr.github.io/rbiom/reference/fn_schema.md)
entry, without running it. The dataset can be an rbiom object or the
summary returned by `fn_check_data()`, so that a program that only keeps
a summary can check requests before passing them on.

## Usage

``` r
fn_check(fn, args = list(), data)

fn_check_data(biom)
```

## Arguments

- fn:

  The name of a function covered by
  [`fn_schema()`](https://cmmr.github.io/rbiom/reference/fn_schema.md).

- args:

  A named list of the arguments to check, as they would be passed to
  `fn`. Arguments not given take their defaults. Arguments marked
  `r_only` in the schema can't be given.

- data:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  the list returned by `fn_check_data()` (also after a round trip
  through JSON).

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

## Value

`fn_check()`: the arguments in `args`, converted as the function would
convert them: abbreviations expanded, rank numbers replaced by rank
names, layer letters by layer names, metric names by metric ids. A field
given with a `"=="` or `"!="` prefix is returned without it, and added
to `within` or `between`, which are then returned too. An invalid
argument is an error, with the message the function would give.

        `fn_check_data()`: a list with `fields` (one entry per
        metadata field, giving its `class`: `"numeric"` or
        `"categorical"`), `ranks` (the rank names), and `tree` (whether
        `biom` has a phylogenetic tree).

## Details

Checks cover the type, choices, limits, and number of each value;
metadata field names and their classes; ranks; metric names (and whether
a phylogenetic metric has a tree); and the constraints in the entry's
`rules`. Checks that need computed data are left to the function: names
in `colors`, `shapes`, and `patterns` against the groups, palette and
color names, the parts of a `tracks` or `grid` list, and taxa names.

## See also

[`fn_schema()`](https://cmmr.github.io/rbiom/reference/fn_schema.md)

## Examples

``` r
    library(rbiom)

    fn_check("adiv_boxplot", list(x = "body", adiv = "shan", layers = "xd"), hmp50)
#> $x
#> [1] "Body Site"
#> 
#> $adiv
#> [1] "shannon"
#> 
#> $layers
#> [1] "box" "dot"
#> 

    data <- fn_check_data(hmp50)
    str(data)
#> List of 3
#>  $ fields:List of 5
#>   ..$ .sample  :List of 1
#>   .. ..$ class: chr "categorical"
#>   ..$ Age      :List of 1
#>   .. ..$ class: chr "numeric"
#>   ..$ BMI      :List of 1
#>   .. ..$ class: chr "numeric"
#>   ..$ Body Site:List of 1
#>   .. ..$ class: chr "categorical"
#>   ..$ Sex      :List of 1
#>   .. ..$ class: chr "categorical"
#>  $ ranks : chr [1:7] ".otu" "Kingdom" "Phylum" "Class" ...
#>  $ tree  : logi TRUE

    fn_check("taxa_boxplot", list(rank = 2, taxa = 3), data)
#> $rank
#> [1] "Phylum"
#> 
#> $taxa
#> [1] 3
#> 

    try(fn_check("adiv_boxplot", list(ci = "zz"), data))
#> Error in validate_var_choices(var = name, choices = spec$choices, env = env,  : 
#>   Invalid `ci` option(s): zz
#> Choices are: ci, range, sd, se, mad
```
