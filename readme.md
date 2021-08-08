# D47crunch

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4314550.svg)](https://doi.org/10.5281/zenodo.4314550)

Python library for processing and standardizing carbonate clumped-isotope analyses, from low-level data out of a dual-inlet mass spectrometer to final, “absolute” Δ<sub>47</sub> and Δ<sub>48</sub> values with fully propagated analytical error estimates.

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
Having done any of the above you should now be able to `import D47crunch`, with the following requirements: [Python 3], [numpy], [matplotlib], [scipy], and [lmfit].

[Python 3]: https://www.python.org
[numpy]: https://numpy.org
[lmfit]: https://lmfit.github.io
[matplotlib]: https://matplotlib.org
[scipy]: https://www.scipy.org
[dev branch]: https://github.com/mdaeron/D47crunch/archive/dev.zip
[latest release]: https://github.com/mdaeron/D47crunch/releases/latest

## Contact

All questions and suggestions are welcome and should be directed at [Mathieu Daëron](mailto:daeron@lsce.ipsl.fr?subject=[D47crunch]).

