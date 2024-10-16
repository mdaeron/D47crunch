# 3. Command-Line Interface (CLI)

Instead of writing Python code, you may directly use the CLI to process raw Δ47 and Δ48 data using reasonable defaults. The simplest way is simply to call:

```txt
D47crunch rawdata.csv
```

This will create a directory named `output` and populate it by calling the following methods:

* `D47data.wg()`
* `D47data.crunch()`
* `D47data.standardize()`
* `D47data.summary()`
* `D47data.table_of_samples()`
* `D47data.table_of_sessions()`
* `D47data.plot_sessions()`
* `D47data.plot_residuals()`
* `D47data.table_of_analyses()`
* `D47data.plot_distribution_of_analyses()`
* `D47data.plot_bulk_compositions()`
* `D47data.save_D47_correl()`

You may specify a custom set of anchors instead of the default ones using the `--anchors` or `-a` option:

```txt
D47crunch -a anchors.csv rawdata.csv
```

In this case, the `anchors.csv` file (you may use any other file name) must have the following format:

```csv
Sample, d13C_VPDB, d18O_VPDB,    D47
 ETH-1,      2.02,     -2.19, 0.2052
 ETH-2,    -10.17,    -18.69, 0.2085
 ETH-3,      1.71,     -1.78, 0.6132
 ETH-4,          ,          , 0.4511
```
 
The samples with non-empty `d13C_VPDB`, `d18O_VPDB`, and `D47` values are used to standardize δ13C, δ18O, and Δ47 values respectively.

You may also provide a list of analyses and/or samples to exclude from the input. This is done with the `--exclude` or `-e` option:

```txt
D47crunch -e badbatch.csv rawdata.csv
```

In this case, the `badbatch.csv` file (again, you may use a different file name) must have the following format:

```csv
UID, Sample
A03
A09
B06
   , MYBADSAMPLE-1
   , MYBADSAMPLE-2
```

This will exclude (ignore) analyses with the UIDs `A03`, `A09`, and `B06`, and those of samples `MYBADSAMPLE-1` and `MYBADSAMPLE-2`. It is possible to have and exclude file with only the `UID` column, or only the `Sample` column, or both, in any order.

The `--output-dir` or `-o` option may be used to specify a custom directory name for the output. For example, in unix-like shells the following command will create a time-stamped output directory:

```txt
D47crunch -o `date "+%Y-%M-%d-%Hh%M"` rawdata.csv
```

To process Δ48 as well as Δ47, just add the `--D48` option.