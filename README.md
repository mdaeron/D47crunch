# D47crunch

Python library for processing and standardizing carbonate clumped-isotope analyses, from low-level data out of a dual-inlet mass spectrometer to final, “absolute” Δ<sub>47</sub> values with fully propagated analytical error estimates.

All questions and suggestions are welcome and should be directed at [Mathieu Daëron](mailto:daeron@lsce.ipsl.fr?subject=[D47crunch]).

## Requirements

Python 3, [numpy], [lmfit]. We recommend installing the [Anaconda] distribution.

[numpy]: https://numpy.org
[lmfit]: https://lmfit.github.io
[Anaconda]: https://www.anaconda.com/distribution 

## Installation

This should do the trick:

```bash
pip install D47crunch
```

Alternatively:

1. download [D47crunch-master.zip]
2. unzip it
3. rename the resulting directory to `D47crunch`
4. move the `D47crunch` directory to somewhere in your `PYTHONPATH` or to your current working directory

[D47crunch-master.zip]: https://github.com/mdaeron/D47crunch/archive/master.zip

## Documentation

For a short tutorial and the full API, see [https://github.com/mdaeron/D47crunch/docs/index.html].

[https://github.com/mdaeron/D47crunch/docs/index.html]: https://github.com/mdaeron/D47crunch/docs/index.html
