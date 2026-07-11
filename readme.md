# D47crunch

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4314550.svg)](https://doi.org/10.5281/zenodo.4314550)
[![PyPI Version](https://img.shields.io/pypi/v/D47crunch.svg)](https://pypi.python.org/pypi/D47crunch)
![License](https://img.shields.io/pypi/l/D47crunch)
[![docs](https://shields.io/badge/docs-D47crunch-success?logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxOTIiIGhlaWdodD0iMTkyIj48cGF0aCBkPSJNNjcuODA3IDIzLjI0NGMtMTguMDMyIDUuNTIzLTI0LjI3MyA5LjQxMS0zMy4xOTMgMTYuMTY5UzE5LjQwNyA1MS4xMjUgOS4zNCA2OS43NkMzLjI3IDgwLjk5Ni44MjcgOTkuMjcxIDEgMTExLjY3OWMuMTE0IDguMTY5IDIuMTg4IDEzLjU1OSAyLjE4OCAxMy41NTlTNDIuMzI2IDI0Ny4yNyAyMC40MzUgMjc4LjEwNWMtMTAuNDg4IDE0Ljc3MyA4NC40MjIgMjcuMDA3IDgwLjA0NyAxMC40NS0xMC4zOTctMzkuMzQ4LTI3LjU5NS05Ny40MDktMzEuOTA4LTExMy42NjctMy4wNDctMTEuNDg0IDguMTYtMjEuOTk1IDE2LjI5NS0yNy4zNDMgMy4zODMtMi4yMjQgOC4wNTEtMy43ODYgOC4wNTEtMy43ODZsNTEuOTI0LTE1LjU5NXMyMy44MDktOC41NzkgMzQuMDYyLTI2LjU2NyAxMi42ODktMzcuNTg0IDEyLjY4OS0zNy41ODQuNDk1LTE2LjU4OC01LjkzOS0yMy4xNC0yLjczMi01Ljg2Ni0xOC42NDYtMTIuNDk4LTQ5LjEyMi0xMC4wNTMtNDkuMTIyLTEwLjA1My0zMi4wNDktLjYwMi01MC4wODEgNC45MjF6IiBmaWxsPSIjMjBhZDZjIi8+PGcgdHJhbnNmb3JtPSJyb3RhdGUoMzQzLjk5MykiPjxlbGxpcHNlIGN4PSI1Mi4yNzYiIGN5PSI5NC4wODkiIGZpbGw9IiNmZmYiIHJ4PSIyNy44NTIiIHJ5PSIyNy44MjMiLz48ZWxsaXBzZSBjeD0iNTIuNTM4IiBjeT0iOTQuNjg5IiByeD0iMTguNjg4IiByeT0iMTguMzMxIiBmaWxsPSIjMTA1YTQ4Ii8+PGVsbGlwc2UgY3g9IjYzLjMzMyIgY3k9Ijg2LjMiIGZpbGw9IiNmZmYiIHJ4PSI3LjU5NiIgcnk9IjcuNTg4Ii8+PC9nPjxnIGZpbGw9IiMxMDVhNDgiPjxlbGxpcHNlIGN4PSI5NC4xNTQiIGN5PSIxMzUuNzM0IiByeD0iNC41MDUiIHJ5PSI1LjI1NyIgdHJhbnNmb3JtPSJyb3RhdGUoMzI1LjAxNSkiLz48ZWxsaXBzZSBjeD0iMTY3Ljg2IiBjeT0iNTYuOTUxIiByeD0iNC40ODQiIHJ5PSI1LjIzMSIgdHJhbnNmb3JtPSJyb3RhdGUoMzU4Ljc3OSkiLz48L2c+PHBhdGggZD0iTTE4OC45NjIgNzcuNDU4bC0zNC4yNDggMTAuMTcxYy02Ljg3OCAxLjk3My05Ljk3NCAzLjM4OC0xNC44ODggMy42NnMtMTAuODY4LTEuODY5LTEwLjg2OC0xLjg2OWwtMi4zNDIgOC42MzNjMi4xOTcuOTQ0IDkuOTY5IDIuMzkgMTQuNzU2IDIuNTkyczE1LjU4Mi0yLjQ2IDE1LjU4Mi0yLjQ2bDI3LjQ3MS03LjUzNGMxLjk1MS00LjcwOSAzLjQzMi05LjIzOCA0LjUzNy0xMy4xOTN6IiBmaWxsPSIjZTc4MzYxIi8+PC9zdmc+)](https://mdaeron.github.io/D47crunch)

Python library for processing and standardizing carbonate clumped-isotope analyses, from low-level data out of a dual-inlet mass spectrometer to final, “absolute” Δ<sub>47</sub>, Δ<sub>48</sub> and  Δ<sub>49</sub> values with fully propagated analytical error estimates.

This also provides a command-line interface making it possible to process a multi-session Δ<sub>47</sub> (and potentially Δ<sub>48</sub>, Δ<sub>49</sub>) dataset by simply calling `D47crunch rawdata.csv`. See the [CLI documentation](https://mdaeron.github.io/D47crunch/#3-command-line-interface-cli) for more options.

> **Warning**
The upcoming version 3 of D47crunch will break backwards-compatibility. This is motivated by the need to accomodate different standardization approaches (least squares "pooled" standardization vs Bayesian standardization allowing for reference materials having nominal Δ4x values with non-zero uncertainties).

> **Note**
The new Bayesian standardization approach introduced in the upcoming version 3 of D47crunch works well, to the first order, but has not been extensively tested yet. To test this upcoming version without affecting your existing Python setup, use `uv` and start with the demo code available [here](https://github.com/mdaeron/D47crunch/tree/bayes/examples).

```
mkdir testD47crunch   # create new working directory
cd testD47crunch      # enter the working directory
uv init               # initialize a new uv project
uv add mylib==3.0.0a1 # install the new version of D47crunch
uv run bayes-demo.py  # run the demo code
```



## Documentation

For the full API and a short tutorial, see [https://mdaeron.github.io/D47crunch].

[https://mdaeron.github.io/D47crunch]: https://mdaeron.github.io/D47crunch

## Installation

This should do the trick:

```bash
pip install D47crunch
```

Alternatively:

1. download the [dev branch] or the [latest release] and unzip it
2. rename the resulting directory to `D47crunch`
3. chose one of one of the following options:
	+ move the `D47crunch` directory to somewhere in your Python path
	+ move the `D47crunch` directory to your desired working directory
	+ move the `D47crunch` directory to any other location (e.g., `/foo/bar`) and include the following code snippet in your scripts:

```py
import sys
sys.path.append('/foo/bar')
```
Having done any of the above you should now be able to `import D47crunch`, with the following requirements: [Python 3], [numpy], [matplotlib], [scipy], [lmfit], and [typer].

[Python 3]: https://www.python.org
[numpy]: https://numpy.org
[lmfit]: https://lmfit.github.io
[matplotlib]: https://matplotlib.org
[scipy]: https://www.scipy.org
[typer]: https://typer.tiangolo.com
[dev branch]: https://github.com/mdaeron/D47crunch/archive/dev.zip
[latest release]: https://github.com/mdaeron/D47crunch/releases/latest

## Contact

All questions and suggestions are welcome and should be directed at [Mathieu Daëron](mailto:daeron@lsce.ipsl.fr?subject=[D47crunch]).

