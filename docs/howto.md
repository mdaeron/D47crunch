## 2. How-to

### 2.1 Use a different set of anchors, change anchor nominal values, and/or change <sup>17</sup>O correction parameters

Nominal values for various carbonate standards are defined in four places:

* `D4xdata.Nominal_d13C_VPDB`
* `D4xdata.Nominal_d18O_VPDB`
* `D47data.Nominal_D4x` (also accessible through `D47data.Nominal_D47`)
* `D48data.Nominal_D4x` (also accessible through `D48data.Nominal_D48`)

<sup>17</sup>O correction parameters are defined by:

* `D4xdata.R13_VPDB`
* `D4xdata.R18_VSMOW`
* `D4xdata.R18_VPDB`
* `D4xdata.lambda_17`
* `D4xdata.R17_VSMOW`
* `D4xdata.R17_VPDB`

When creating a new instance of `D47data` or `D48data`, the current values of these variables are copied as properties of the new object. Applying custom values for, e.g., `R17_VSMOW` and `Nominal_D47` can thus be done in several ways:

**Option 1:** by redefining `D4xdata.R17_VSMOW` and `D47data.Nominal_D47` _before_ creating a `D47data` object:

```python
from D47crunch import D4xdata, D47data

# redefine R17_VSMOW:
D4xdata.R17_VSMOW = 0.00037 # new value

# redefine R17_VPDB for consistency:
D4xdata.R17_VPDB = D4xdata.R17_VSMOW * (D4xdata.R18_VPDB/D4xdata.R18_VSMOW) ** D4xdata.lambda_17

# edit Nominal_D47 to only include ETH-1/2/3:
D47data.Nominal_D4x = {
	a: D47data.Nominal_D4x[a]
	for a in ['ETH-1', 'ETH-2', 'ETH-3']
	}
# redefine ETH-3:
D47data.Nominal_D4x['ETH-3'] = 0.600

# only now create D47data object:
mydata = D47data()

# check the results:
print(mydata.R17_VSMOW, mydata.R17_VPDB)
print(mydata.Nominal_D47)
# NB: mydata.Nominal_D47 is just an alias for mydata.Nominal_D4x

# should print out:
# 0.00037 0.00037599710894149464
# {'ETH-1': 0.2052, 'ETH-2': 0.2085, 'ETH-3': 0.6}
```

**Option 2:** by redefining `R17_VSMOW` and `Nominal_D47` _after_ creating a `D47data` object:

```python
from D47crunch import D47data

# first create D47data object:
mydata = D47data()

# redefine R17_VSMOW:
mydata.R17_VSMOW = 0.00037 # new value

# redefine R17_VPDB for consistency:
mydata.R17_VPDB = mydata.R17_VSMOW * (mydata.R18_VPDB/mydata.R18_VSMOW) ** mydata.lambda_17

# edit Nominal_D47 to only include ETH-1/2/3:
mydata.Nominal_D47 = {
	a: mydata.Nominal_D47[a]
	for a in ['ETH-1', 'ETH-2', 'ETH-3']
	}
# redefine ETH-3:
mydata.Nominal_D47['ETH-3'] = 0.600

# check the results:
print(mydata.R17_VSMOW, mydata.R17_VPDB)
print(mydata.Nominal_D47)

# should print out:
# 0.00037 0.00037599710894149464
# {'ETH-1': 0.2052, 'ETH-2': 0.2085, 'ETH-3': 0.6}
```

The two options above are equivalent, but the latter provides a simple way to compare different data processing choices:

```python
from D47crunch import D47data

# create two D47data objects:
foo = D47data()
bar = D47data()

# modify foo in various ways:
foo.lambda_17 = 0.52
foo.R17_VSMOW = 0.00037 # new value
foo.R17_VPDB = foo.R17_VSMOW * (foo.R18_VPDB/foo.R18_VSMOW) ** foo.lambda_17
foo.Nominal_D47 = {
	'ETH-1': foo.Nominal_D47['ETH-1'],
	'ETH-2': foo.Nominal_D47['ETH-1'],
	'IAEA-C2': foo.Nominal_D47['IAEA-C2'],
	'INLAB_REF_MATERIAL': 0.666,
	}

# now import the same raw data into foo and bar:
foo.read('rawdata.csv')
foo.wg()          # compute δ13C, δ18O of working gas
foo.crunch()      # compute all δ13C, δ18O and raw Δ47 values
foo.standardize() # compute absolute Δ47 values

bar.read('rawdata.csv')
bar.wg()          # compute δ13C, δ18O of working gas
bar.crunch()      # compute all δ13C, δ18O and raw Δ47 values
bar.standardize() # compute absolute Δ47 values

# and compare the final results:
foo.table_of_samples(verbose = True, save_to_file = False)
bar.table_of_samples(verbose = True, save_to_file = False)
```

### 2.2 Simulate a virtual data set to play with

It is sometimes convenient to quickly build a virtual data set of analyses, for instance to assess the final analytical precision achievable for a given combination of anchor and unknown analyses (see also Fig. 6 of [Daëron, 2021]).

This can be achieved with `virtual_data()`. The example below creates a dataset with four sessions, each of which comprises four analyses of anchor ETH-1, five of ETH-2, six of ETH-3, and two analyses of an unknown sample named `FOO` with an arbitrarily defined isotopic composition. Analytical repeatabilities for Δ<sub>47</sub> and Δ<sub>48</sub> are also specified arbitrarily. See the `virtual_data()` documentation for additional configuration parameters.

```py
from D47crunch import *

args = dict(
	samples = [
		dict(Sample = 'ETH-1', N = 4),
		dict(Sample = 'ETH-2', N = 5),
		dict(Sample = 'ETH-3', N = 6),
		dict(
			Sample = 'FOO',
			N = 2,
			d13C_VPDB = -5.,
			d18O_VPDB = -10.,
			D47 = 0.3,
			D48 = 0.15
			),
		],
	rD47 = 0.010,
	rD48 = 0.030,
	)

session1 = virtual_data(session = 'Session_01', **args)
session2 = virtual_data(session = 'Session_02', **args)
session3 = virtual_data(session = 'Session_03', **args)
session4 = virtual_data(session = 'Session_04', **args)

D = D47data(session1 + session2 + session3 + session4)

D.crunch()
D.standardize()

D.table_of_sessions(verbose = True, save_to_file = False)
D.table_of_samples(verbose = True, save_to_file = False)
D.table_of_analyses(verbose = True, save_to_file = False)
```

### 2.3 Process paired Δ<sub>47</sub> and Δ<sub>48</sub> values

