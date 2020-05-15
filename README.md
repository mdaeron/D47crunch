# D47crunch

Python library for processing and standardizing carbonate clumped-isotope analyses, from low-level data out of a dual-inlet mass spectrometer to final, “absolute” Δ<sub>47</sub> values with fully propagated analytical error estimates.

All questions and suggestions are welcome and should be directed at [Mathieu Daëron](mailto:daeron@lsce.ipsl.fr?subject=[D47crunch]).

## Requirements

[Python 3], [numpy], [matplotlib], and [lmfit]. For the first three we recommend installing the [Anaconda] distribution. Installing [lmfit] should be a simple as `pip install lmfit`.

[Python 3]: https://www.python.org
[numpy]: https://numpy.org
[lmfit]: https://lmfit.github.io
[Anaconda]: https://www.anaconda.com/distribution
[matplotlib]: https://matplotlib.org

## Installation

This should do the trick:

```bash
pip install D47crunch
```

Alternatively:

1. download the [current branch] or the [latest release] and unzip it
2. rename the resulting directory to `D47crunch`
3. chose one of one of the following options:
	+ move the `D47crunch` directory to somewhere in your Python path
	+ move the `D47crunch` directory to your desired working directory
	+ move the `D47crunch` directory to any other location (e.g., `/foo/bar`) and include the following code snippet in your scripts:

```py
import sys
sys.path.append('/foo/bar')
```
Having done any of the above you should now be able to `import D47crunch`.

[current branch]: https://github.com/mdaeron/D47crunch/archive/master.zip
[latest release]: https://github.com/mdaeron/D47crunch/archive/v0.1.zip

## Documentation

For a short tutorial and the full API, see [https://mdaeron.github.io/D47crunch].

[https://mdaeron.github.io/D47crunch]: https://mdaeron.github.io/D47crunch
