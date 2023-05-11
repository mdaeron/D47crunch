# 1. Tutorial

## 1.1 Installation

The easy option is to use `pip`; open a shell terminal and simply type:

```
python -m pip install D47crunch
```

For those wishing to experiment with the bleeding-edge development version, this can be done through the following steps:

1. Download the `dev` branch source code [here](https://raw.githubusercontent.com/mdaeron/D47crunch/dev/D47crunch/__init__.py) and rename it to `D47crunch.py`.
2. Do any of the following:
    * copy `D47crunch.py` to somewhere in your Python path
    * copy `D47crunch.py` to a working directory (`import D47crunch` will only work if called within that directory)
    * copy `D47crunch.py` to any other location (e.g., `/foo/bar`) and then use the following code snippet in your own code to import `D47crunch`:

```py
import sys
sys.path.append('/foo/bar')
import D47crunch
```

Documentation for the development version can be downloaded [here](https://github.com/mdaeron/D47crunch/raw/dev/docs/index.html) (save html file and open it locally).

## 1.2 Usage

Start by creating a file named `rawdata.csv` with the following contents:

```html
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

```py
import D47crunch
mydata = D47crunch.D47data()
```

For now, this object is empty:

```html
>>> print(mydata)
[]
```

To load the analyses saved in `rawdata.csv` into our `D47data` object and process the data:

```py
mydata.read('rawdata.csv')

# compute δ13C, δ18O of working gas:
mydata.wg()

# compute δ13C, δ18O, raw Δ47 values for each analysis:
mydata.crunch()

# compute absolute Δ47 values for each analysis
# as well as average Δ47 values for each sample:
mydata.standardize()
```

We can now print a summary of the data processing:

```html
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

This tells us that our data set contains 5 different samples: 3 anchors (ETH-1, ETH-2, ETH-3) and 2 unknowns (MYSAMPLE-1, MYSAMPLE-2). The total number of analyses is 8, with 5 anchor analyses and 3 unknown analyses. We get an estimate of the analytical repeatability (i.e. the overall, pooled standard deviation) for δ13C, δ18O and Δ47, as well as the number of degrees of freedom (here, 3) that these estimated standard deviations are based on, along with the corresponding Student's t-factor (here, 3.18) for 95&nbsp;% confidence limits. Finally, the summary indicates that we used a “pooled” standardization approach (see [Daëron, 2021]).

To see the actual results:

```html
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

This table lists, for each sample, the number of analytical replicates, average δ13C and δ18O values (for the analyte CO2 , *not* for the carbonate itself), the average Δ47 value and the SD of Δ47 for all replicates of this sample. For unknown samples, the SE and 95 % confidence limits for mean Δ47 are also listed These 95 % CL take into account the number of degrees of freedom of the regression model, so that in large datasets the 95 % CL will tend to 1.96 times the SE, but in this case the applicable t-factor is much larger.

We can also generate a table of all analyses in the data set (again, note that `d18O_VSMOW` is the composition of the CO2 analyte):

```html
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

