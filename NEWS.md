---
  output:
  html_document: default
pdf_document: default
---
  `MatchingFrontier` News and Updates
======

# MatchingFrontier (development version)

`MatchingFrontier` has been completely rewritten from scratch, providing more consistent syntax, new methods for constructing the frontier, and many other new features. All of these are breaking changes, meaning you should not expect result to agree with those from version 3.0.0 or below, and syntax from earlier versions may not work with this version. To install an older version of `MatchingFrontier` (e.g., to reproduce the results of an older analysis), you can install an archived version from CRAN. Here, we describe some of the changes made to the package and its functionality.

### Constructing the frontier

`makeFrontier()` has been completely rewritten. The biggest changes come in the new imbalance metrics and quantities of interest allowed and the changing of some past options. Factor variables can now be used as well.

* `QOI` can now be one of `"SATE"` (sample average treatment effect), `"FSATE"` (feasible sample average treatment effect), `"SATT"` (sample average treatment effect in the treated), or `"FSATT"` (feasible sample average treatment effect in the treated).

* The different types of frontiers are now available depending on the imbalance metric supplied to `metric`. These include pair distance-based frontiers, bin-based frontiers, and energy distance-based frontiers. For pair distance-based frontiers, `metric` can be one of `"Mahal"`, `"Euclid"`, or "`Custom"`. For bin-based frontiers, `metric` can be one `"L1"`, `"L2"`, `"L1median"`, or `"L2median"`. For energy distance-based frontiers, `metric` can be `"Energy"`. Not all `metric`s are available with all `QOI`s, but the ones available have expanded from previous versions.

* Bugs in how the balance statistics are computed have been fixed.

* Speed improvements for pair distance-based and bin-based frontiers.

* The `ratio` argument is now defunct and has been removed.

* The arguments now be specified either as a dataset and strings for the matching covariates and treatment (as in prior versions) or as a formula and dataset.

* A `verbose` argument has been added to control printing of results.

* A new `print()` method displays more metadata about the frontier.

* A new `summary()` method displays information about the number of units dropped and the imbalance measure at several points on the frontier.

* The `drop.order` component of a `matchFrontier` object now begins with an empty entry, indicating no units have been dropped, and ends just before the final unit(s) are dropped, to match with the `Xs` and `Ys` components.

### Plotting the frontier

`plot.matchFrontier()` has been rewritten to rely on `ggplot2` instead of the base R plotting functions, making it easier to customize the plots after producing them by simply using `ggplot2` syntax on the resulting object. `plot.matchFrontier()` also subsumes `plotMeans()` and has even more functionality.

* The resulting plot now has two x-axes: one displaying the number of units dropped, and one displaying the number of units remaining.

* The plot is now a line plot rather than a dot plot.

* The `covs` argument can be specified to produce balance plots for specific covariates; these include plots of (standardized) mean differences, (standardized) means, or Kolmogorov-Smirnov statistics depending on the arguments supplied to `stat`.

* There is now a help page for `plot.matchFrontier()`.

### Extracting data from one point on the frontier

A new function, `frontier_to_matchit()`, has been added to convert the solution from one point on the frontier to a `matchit` object to be further processed by `MatchIt` functions, like `summary.matchit()` and `plot.matchit()` for assessing balance. The utility `generateDataset()` returns a dataset that is identical to that produced if `MatchIt::match.data()` or `MatchIt::get_matches()` were used on the resulting `matchit` object.

* `generateDataset()` has a new `Ndrop` argument to request a dataset after dropping `Ndrop` units. This can be used instead of the `N` argument (which corresponds to the units *remaining* in the matched dataset).

* `generateDataset()` has a new `dup` argument to control whether units that have been matched more than once when using  pair distance-based frontier should appear more than one time in the matched dataset (`dup = TRUE`) or should only appear once but be weighted to reflect their reuse (`dup = FALSE`). Using `dup = TRUE` is equivalent to using `MatchIt::get_matches()` on the `matchit` object resulting from `frontier_to_matchit()`, and using `dup = FALSE` is equivalent to using `MatchIt::match.data()` on the `matchit` object resulting from `frontier_to_matchit()`. The `dup` argument is ignored for other types of frontiers.

* To be in line with `MatchIt::match.data()` or `MatchIt::get_matches()`, `generateDataset()` accepts three new arguments, `weights`, `subclass`, and `id`, that control the names of these components in the resulting dataset. `subclass` and `id` are only used when `dup = TRUE` to rename the subclass (i.e., pair membership) variable and ID (i.e., unit ID) variables, respectively. Weights are always produced (though they are often all 1), so the `weights` argument can always be used. The argument supplied will be stored in corresponding attributes in the resulting object.

### Assessing model dependence

The `modelDependence()` function for assessing model dependence in a single dataset has not changed drastically. The order of arguments has changed to more closely resemble traditional modeling functions (e.g., by allowing a formula to be supplied as the first argument).

* `modelDependence()` now accepts a `weights` argument to supply sampling or matching weights to the fit models.

* A `method` argument has been added to allow users to select the model dependence method rather than guessing based on the other supplied inputs.

* The output is now always a pair of mode dependence bounds; previously, it would return a pair of values for the extreme bounds procedure and a single value for the Athey-Imbens procedure.

* When using the extreme bounds procedure, the 2.5th and 97.5th percentiles of the resulting distribution of effect estimates are returned rather than the minimum and maximum in order to encourage the use of larger number of model specifications.

* A new `plot()` method can be used to visualize the resulting model dependence estimates.

### Estimating effects across the frontier

`estimateEffects()` has been updated to estimate effects in accordance with modern recommendations, e.g., by using robust standard errors and corrected formulas for computing the matching weights. Some of the argument names have changed to be more consistent with other functions.

* The new `method` argument controls whether model dependence bounds are computed or not and which method to use if so. Model dependence bounds can be suppressed by setting `method = "none"`.

* Athey-Imbens model dependence bounds can no longer be used with pair distance-based frontiers. The extreme bounds procedure can be used with any frontier.

* The `alpha` argument now refers to the alpha level used to select the critical test statistic for the confidence intervals rather than the confidence level.

* The new `cutpoint.method` argument controls how covariate cut points when using the Athey-Imbens model dependence procedure are calculated, and can be specified either as `"segmented"` to use segmented regression or `"mean"` or `"median"` to use the covariate mean or median, respectively.

* A new `print()` method displays some metadata about the estimated effects.

### Plotting estimated effects

`plot.frontierEstimates()` has been rewritten to rely on `ggplot2` instead of the base R plotting functions, making it easier to customize the plots after producing them by simply using `ggplot2` syntax on the resulting object.

* Like `plot.matchFrontier()`, the resulting plot now has two x-axes: one displaying the number of units dropped, and one displaying the number of units remaining.

* The new `band` argument controls which type of band is displayed around the estimated effects. Can be `"none"` for no band, `"confidence"` for confidence bands, or `"model-dependence"` for model dependence bands if previously requested.

* New arguments have been added to control the appearance of the bands.

### Documentation

* All help pages have been completely rewritten.

* The vignette has been rewritten to address the recent changes.

* A README has been added to the package.

### Other changes

* The functions `parallelPlot()`, `plotMeans()`, and `plotPrunedMeans()` have been removed. The functionality of `plotMeans()` is subsumed by `plot.matchFrontier()`, and diagnostics on individual matched datasets can be done on the `matchit` object resulting from `frontier_to_matchit()` instead of `parallelPlot()`.

* The `lalonde` dataset that used to come with the package has been removed. Now, the `lalonde` dataset that comes with `MatchIt` is used in the examples and vignettes. This is to reduce the number of `lalonde` datasets floating around that have different numbers of units and variables and correspond to different populations.

* `MatchingFrontier` now lists `MatchIt`, `ggplot2`, `scales`, `lmtest`, and `sandwich` in its `Imports` and `segmented` in its `Suggests`. Importantly, this removes the dependence on `igraph`, which is no longer used.
