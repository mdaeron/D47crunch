# Changelog

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