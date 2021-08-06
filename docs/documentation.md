The **tutorial** section takes you through a series of simple steps to import/process data and print out the results.
The **how-to** section provides instructions applicable to various specific tasks.

## 1. Tutorial

Start by creating a file named `rawdata.csv` with the following contents:

```csv
UID,  Sample,           d45,       d46,        d47,        d48,       d49
A01,  ETH-1,        5.79502,  11.62767,   16.89351,   24.56708,   0.79486
A02,  MYSAMPLE-1,   6.21907,  11.49107,   17.27749,   24.58270,   1.56318
A03,  ETH-2,       -6.05868,  -4.81718,  -11.63506,  -10.32578,   0.61352
A04,  MYSAMPLE-2,  -3.86184,   4.94184,    0.60612,   10.52732,   0.57118
A05,  ETH-3,        5.54365,  12.05228,   17.40555,   25.96919,   0.74608
A06,  ETH-2,       -6.06706,  -4.87710,  -11.69927,  -10.64421,   1.61234
A07,  ETH-1,        5.78821,  11.55910,   16.80191,   24.56423,   1.47963
A08,  MYSAMPLE-2,  -3.87692,   4.86889,    0.52185,   10.40390,   1.07032
```

Then instantiate a `D47data` object which will store and process this data:

```python
import D47crunch
mydata = D47crunch.D47data()
```

For now, this object is empty:

```python
>>> print(mydata)
[]
```

To load the analyses saved in `rawdata.csv` into our `D47data` object and process the data:

```python
mydata.read('rawdata.csv')
mydata.wg()          # compute δ13C, δ18O of working gas
mydata.crunch()      # compute all δ13C, δ18O and raw Δ47 values
mydata.standardize() # compute absolute Δ47 values
```

We can now print a summary of the data processing:

```csv
>>> mydata.summary(verbose = True, save_to_file = False)
[summary]        
–––––––––––––––––––––––––––––––  –––––––––
N samples (anchors + unknowns)   5 (3 + 2)
N analyses (anchors + unknowns)  8 (5 + 3)
Repeatability of δ13C_VPDB         4.2 ppm
Repeatability of δ18O_VSMOW       47.5 ppm
Repeatability of Δ47 (anchors)    13.4 ppm
Repeatability of Δ47 (unknowns)    2.5 ppm
Repeatability of Δ47 (all)         9.6 ppm
Model degrees of freedom                 3
Student's 95% t-factor                3.18
Standardization method              pooled
–––––––––––––––––––––––––––––––  –––––––––
```

This tells us that our data set contains 5 different samples: 3 anchors (ETH-1, ETH-2, ETH-3) and 2 unknowns (MYSAMPLE-1, MYSAMPLE-2). The total number of analyses is 8, with 5 anchor analyses and 3 unknown analyses. We get an estimate of the analytical repeatability (i.e. the overall, pooled standard deviation) for δ<sup>13</sup>C, δ<sup>18</sup>O and Δ<sub>47</sub>, as well as the number of degrees of freedom (here, 3) that these estimated standard deviations are based on, along with the corresponding Student's t-factor (here, 3.18) for 95&nbsp;% confidence limits. Finally, the summary indicates that we used a “pooled” standardization approach (see [Daëron, 2021](https://dx.doi.org/10.1029/2020GC009592)).

To see the actual results:

```csv
>>> mydata.table_of_samples(verbose = True, save_to_file = False)
[table_of_samples] 
––––––––––  –  –––––––––  ––––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––
Sample      N  d13C_VPDB  d18O_VSMOW     D47      SE    95% CL      SD  p_Levene
––––––––––  –  –––––––––  ––––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––
ETH-1       2       2.01       37.01  0.2052                    0.0131          
ETH-2       2     -10.17       19.88  0.2085                    0.0026          
ETH-3       1       1.73       37.49  0.6132                                    
MYSAMPLE-1  1       2.48       36.90  0.2996  0.0091  ± 0.0291                  
MYSAMPLE-2  2      -8.17       30.05  0.6600  0.0115  ± 0.0366  0.0025          
––––––––––  –  –––––––––  ––––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––
```

This table lists, for each sample, the number of analytical replicates, average δ<sup>13</sup>C and δ<sup>18</sup>O values (for the analyte CO<sub>2</sub> , _not_ for the carbonate itself), the average Δ<sub>47</sub> value and the SD of Δ<sub>47</sub> for all replicates of this sample. For unknown samples, the SE and 95 % confidence limits for mean Δ<sub>47</sub> are also listed These 95 % CL take into account the number of degrees of freedom of the regression model, so that in large datasets the 95 % CL will tend to 1.96 times the SE, but in this case the applicable t-factor is much larger.

We can also generate a table of all analyses in the data set (again, note that `d18O_VSMOW` is the composition of the CO<sub>2</sub> analyte):

```csv
>>> mydata.table_of_analyses(verbose = True, save_to_file = False)
[table_of_analyses] 
–––  –––––––––  ––––––––––  –––––––––––  ––––––––––––  –––––––––  –––––––––  ––––––––––  ––––––––––  ––––––––  ––––––––––  ––––––––––  –––––––––  –––––––––  ––––––––––  ––––––––
UID    Session      Sample  d13Cwg_VPDB  d18Owg_VSMOW        d45        d46         d47         d48       d49   d13C_VPDB  d18O_VSMOW     D47raw     D48raw      D49raw       D47
–––  –––––––––  ––––––––––  –––––––––––  ––––––––––––  –––––––––  –––––––––  ––––––––––  ––––––––––  ––––––––  ––––––––––  ––––––––––  –––––––––  –––––––––  ––––––––––  ––––––––
A01  mySession       ETH-1       -3.807        24.921   5.795020  11.627670   16.893510   24.567080  0.794860    2.014086   37.041843  -0.574686   1.149684  -27.690250  0.214454
A02  mySession  MYSAMPLE-1       -3.807        24.921   6.219070  11.491070   17.277490   24.582700  1.563180    2.476827   36.898281  -0.499264   1.435380  -27.122614  0.299589
A03  mySession       ETH-2       -3.807        24.921  -6.058680  -4.817180  -11.635060  -10.325780  0.613520  -10.166796   19.907706  -0.685979  -0.721617   16.716901  0.206693
A04  mySession  MYSAMPLE-2       -3.807        24.921  -3.861840   4.941840    0.606120   10.527320  0.571180   -8.159927   30.087230  -0.248531   0.613099   -4.979413  0.658270
A05  mySession       ETH-3       -3.807        24.921   5.543650  12.052280   17.405550   25.969190  0.746080    1.727029   37.485567  -0.226150   1.678699  -28.280301  0.613200
A06  mySession       ETH-2       -3.807        24.921  -6.067060  -4.877100  -11.699270  -10.644210  1.612340  -10.173599   19.845192  -0.683054  -0.922832   17.861363  0.210328
A07  mySession       ETH-1       -3.807        24.921   5.788210  11.559100   16.801910   24.564230  1.479630    2.009281   36.970298  -0.591129   1.282632  -26.888335  0.195926
A08  mySession  MYSAMPLE-2       -3.807        24.921  -3.876920   4.868890    0.521850   10.403900  1.070320   -8.173486   30.011134  -0.245768   0.636159   -4.324964  0.661803
–––  –––––––––  ––––––––––  –––––––––––  ––––––––––––  –––––––––  –––––––––  ––––––––––  ––––––––––  ––––––––  ––––––––––  ––––––––––  –––––––––  –––––––––  ––––––––––  ––––––––
```

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

It is sometimes convenient to quickly build a virtual data set of analyses, for instance to assess the final analytical precision achievable for a given combination of anchor and unknown analyses (see also Fig. 6 in
[Daëron, 2021](https://dx.doi.org/10.1029/2020GC009592)).

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

## 3. Discussions
