# Changelog

## v2.3.1
*Released on 2023-07-22*

### Improvements
* Better help text and examples for CLI

## v2.3.0
*Released on 2023-07-21*

### New feature
* New method `D4xdata.save_D4x_correl()` to export a list of Δ<sub>4x</sub> values along with their SE and correlation matrix.

### Other changes
* The CLI now also calls `save_D4x_correl()`.

## v2.2.1
*Released on 2023-07-20*

### New feature
* The CLI now processes Δ48 as well as Δ47 data, thanks to the `--D48` option.


## v2.2.0
*Released on 2023-07-20*

### Command-line interface (CLI)
* Rejoice, you no longer need to know Python: it is now possible to process a multi-session Δ47 dataset with custom anchors and custom UID/sample exclusion list by simply calling `D47crunch -a anchors.csv -e exclude.csv -o outdir rawdata.csv`.

### New features
* Add `yspan` option to `D4xdata.plot_residuals()`
* Added `shuffle` option to `virtual_data()`.
* Added `filetype` option to `D4xdata.plot_sessions()`.
* Added `dpi` option to `D4xdata.plot_sessions()`, `D4xdata.plot_residuals()`.

### Bugfix
* Fix error in `D4xdata.plot_residuals()` when `hist = False` and `kde = False`.

## v2.1.1
*Released on 2023-05-16*

### Bugfix
* `D4xdata.compute_r()` uses an improved computation for degrees of freedom for arbitrary subsets of sessions and/or samples, yielding more estimates of analytical repeatabilities for Δ47 and Δ48.

### New feature
* Added `kde` option to `D4xdata.plot_residuals()`

### Other changes
* Minor improvement to y axis tick labels) in `D4xdata.plot_residuals()`.

## v2.1.0
*Released on 2023-05-14*

### New feature
* `D4xdata.plot_bulk_compositions()` plots the dispersion of δ13C and δ18O values for each sample.

## v2.0.6
*Released on 2023-05-13*

### Bugfix
* Eliminate some spurious debugging messages in `_fullcovar()`

## v2.0.5
*Released on 2023-05-11*

### Changes
* Under the hood: constrained parameters in pooled standardization now get fully propagated variance and covariance, allowing for truly arbitrary constraints without having book-keeping problems further down the line.

## v2.0.4
*Released on 2023-05-11*

### Changes
* Graphically improved `D4xdata.plot_distribution_of_analyses()`

### Bugfix
* Fix `D4xdata.standardize()` when using weighted sessions

## v2.0.3
*Released on 2022-02-27*

### New feature
* `D4xdata.covar_table()` allows exporting the variance-covariance matrix or the correlation matrix for the Δ<sub>4x</sub> values of unknwon samples.

### Changes
* New `hist` keyword to `D4xdata.plot_residuals()`, which adds a histogram of residuals to the side of the plot.

## v2.0.2
*Released on 2021-08-16*

### Internals
* Remove HTML tags in all docstrings

## v2.0.1
*Released on 2021-08-08*

### Bugfix
* Fix silly mistake in readme.

## v2.0.0
*Released on 2021-08-08*

### New feature
* Support for Δ<sub>48</sub> standardization (cf section *Process paired Δ<sub>47</sub> and Δ<sub>48</sub> values* in the documentation).

### Changes
* Extensive changes to the documentation, with new sections (*Tutorial* and *How-to*).
* Documentation is now built with `pdoc` instead of `pdoc3`.
* `D47data.simulate()` replaced by `simulate_single_analysis()` and `virtual_data()`, with additional functionality.
* New method: `D4xdata.report()`
* `D4xdata.table_of_analyses()`, `D4xdata.table_of_sessions()`, and `D4xdata.table_of_samples()` have a new argument `output` controlling what the method should return.
* `SAMPLE_CONSTRAINING_WG_COMPOSITION` is gone, replaced by `D4xdata.Nominal_d13C_VPDB` and `D4xdata.Nominal_d18O_VPDB`.
* `lambda_17` is now replaced by `LAMBDA_17` everywhere.
* Additional tests

### Bugfixes
* Correct (or at least improve) computations of analytical repeatabilities in `D47data.repeatabilities()` and `D47data.compute_r()`.