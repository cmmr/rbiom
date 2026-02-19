# Suggest Rarefaction Depth

Calculates a rarefaction depth that balances retaining samples against
retaining total observations.

## Usage

``` r
suggest_rarefy_depth(biom)
```

## Arguments

- biom:

  An [rbiom
  object](https://cmmr.github.io/rbiom/reference/rbiom_objects.md), or
  any value accepted by
  [`as_rbiom()`](https://cmmr.github.io/rbiom/reference/as_rbiom.md).

## Value

A single integer representing the suggested rarefaction depth.

## Heuristic

This function selects a depth by analyzing the trade-off between
dropping samples (to increase depth) and lowering depth (to keep
samples).

1.  **Calculate Yields:** For every distinct sample depth in the
    dataset, calculate the total number of observations that would
    remain if the dataset were rarefied to that level. \$\$Yield_d = d
    \times N\_{\ge d}\$\$ Where \\d\\ is the depth and \\N\_{\ge d}\\ is
    the number of samples with at least that many reads.

2.  **Define Threshold:** Calculate 10% of the total observations in the
    original un-rarefied dataset.

3.  **Select Depth:** Find the **lowest** depth \\d\\ where the
    \\Yield_d\\ exceeds this 10% threshold.

This approach prioritizes keeping as many samples as possible, provided
that doing so doesn't discard more than 90% of the dataset's total
information.

## See also

[`rarefy()`](https://cmmr.github.io/rbiom/reference/rarefy.md) which
uses this heuristic when `depth = NULL`.

## Examples

``` r
    library(rbiom)
    
    suggest_rarefy_depth(hmp50)
#> [1] 1183
```
