# D47crunch

Python library for processing and standardizing carbonate clumped-isotope analyses, from low-level data out of a dual-inlet mass spectrometer to final, “absolute” Δ<sub>47</sub> values with fully propagated analytical error estimates.

All questions and suggestions are welcome and should be directed at [Mathieu Daëron](mailto:daeron@lsce.ipsl.fr?subject=[D47crunch]).

## 1. Requirements

Python 3, [numpy], [lmfit]. We recommend installing the [Anaconda] distribution.

[numpy]: https://numpy.org
[lmfit]: https://lmfit.github.io
[Anaconda]: https://www.anaconda.com/distribution 

## 2. Installation

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

## 3. Documentation

For the full API documentation, see [https://github.com/mdaeron/D47crunch/docs/index.html].

For a short tutorial see below.

[https://github.com/mdaeron/D47crunch/docs/index.html]: https://github.com/mdaeron/D47crunch/docs/index.html

## 4. Usage

### 4.1 Import data

Start with some raw data stored as CSV in a file named `rawdata.csv` (spaces after commas are optional). Each line corresponds to a single analysis.

The only required fields are a sample identifier (`Sample`), and the working-gas delta values `d45`, `d46`, `d47`. If no session information is provided, all analuses will be treated as belonging to a single analytical session. Alternatively, to group analyses into sessions, provide session identifiers in a `Session` field. If not specified by the user, a unique identifier (`UID`) will be assigned automatically to each analysis. Independently known oxygen-17 anomalies may be provided as `D17O` (in ‰ relative to VSMOW, with λ equal to `D47data.lambda_17`), and are assumed to be zero otherwise. Working-gas deltas `d48` and `d49` may also be provided, and are otherwise treated as `NaN`.

Example `rawdata.csv` file:

```
UID,  Session,  Sample,       d45,      d46,       d47,       d48,      d49
A01, Session1,   ETH-1,   5.79502, 11.62767,  16.89351,  24.56708,  0.79486
A02, Session1, IAEA-C1,   6.21907, 11.49107,  17.27749,  24.58270,  1.56318
A03, Session1,   ETH-2,  -6.05868, -4.81718, -11.63506, -10.32578,  0.61352
A04, Session1, IAEA-C2,  -3.86184,  4.94184,   0.60612,  10.52732,  0.57118
A05, Session1,   ETH-3,   5.54365, 12.05228,  17.40555,  25.96919,  0.74608
A06, Session1,   ETH-2,  -6.06706, -4.87710, -11.69927, -10.64421,  1.61234
A07, Session1,   ETH-1,   5.78821, 11.55910,  16.80191,  24.56423,  1.47963
A08, Session1, IAEA-C2,  -3.87692,  4.86889,   0.52185,  10.40390,  1.07032
A09, Session1,   ETH-3,   5.53984, 12.01344,  17.36863,  25.77145,  0.53264
A10, Session1, IAEA-C1,   6.21905, 11.44785,  17.23428,  24.30975,  1.05702
A11, Session2,   ETH-1,   5.79958, 11.63130,  16.91766,  25.12232,  1.25904
A12, Session2, IAEA-C1,   6.22514, 11.51264,  17.33588,  24.92770,  2.54331
A13, Session2,   ETH-2,  -6.03042, -4.74644, -11.52551, -10.55907,  0.04024
A14, Session2, IAEA-C2,  -3.83702,  4.99278,   0.67529,  10.73885,  0.70929
A15, Session2,   ETH-3,   5.53700, 12.04892,  17.42023,  26.21793,  2.16400
A16, Session2,   ETH-2,  -6.06820, -4.84004, -11.68630, -10.72563,  0.04653
A17, Session2,   ETH-1,   5.78263, 11.57182,  16.83519,  25.09964,  1.26283
A18, Session2, IAEA-C2,  -3.85355,  4.91943,   0.58463,  10.56221,  0.71245
A19, Session2,   ETH-3,   5.52227, 12.01174,  17.36841,  26.19829,  1.03740
A20, Session2, IAEA-C1,   6.21937, 11.44701,  17.26426,  24.84678,  0.76866
```

First create a `D47data` object named `foo` and import `rawdata.csv`:

```python
import D47crunch
 
foo = D47crunch.D47data()
foo.read('rawdata.csv')
    
print('foo contains:')
print(f'{len(foo)} analyses')
print(f'{len({r["Sample"] for r in foo})} samples')
print(f'{len({r["Session"] for r in foo})} sessions')

# output:
# foo contains:
# 20 analyses
# 5 samples
# 2 sessions
```

We can inspect the elements of `foo`:

```python
r = foo[0]
for k in r:
    print(f'r["{k}"] = {repr(r[k])}')

# output:
# r["UID"] = 'A01'
# r["Session"] = 'Session1'
# r["Sample"] = 'ETH-1'
# r["d45"] = 5.79502
# r["d46"] = 11.62767
# r["d47"] = 16.89351
# r["d48"] = 24.56708
# r["d49"] = 0.79486
```

### 4.2 Working gas composition

There are two ways to define the isotpic composition of the working gas.

#### 4.2.1 Option 1: explicit definition

Directly writing to fields `d13Cwg_VPDB` and `d18Owg_VSMOW`:

```python
for r in foo:
    if r['Session'] == 'Session1':
        r['d13Cwg_VPDB'] = -3.75
        r['d18Owg_VSMOW'] = 25.14
    elif r['Session'] == 'Session2':
        r['d13Cwg_VPDB'] = -3.74
        r['d18Owg_VSMOW'] = 25.17
```

#### 4.2.2 Option 2: based on the known composition of a sample:

```python
# The 2 code lines below are the default settings. It is thus not
# necessary to include them unless you wish to use different values.

foo.SAMPLE_CONSTRAINING_WG_COMPOSITION = ('ETH-3', 1.71, -1.78)
foo.ALPHA_18O_ACID_REACTION = 1.00813 # (Kim et al., 2007), calcite at 90 °C

# Compute the WG composition for each session:
foo.wg()

```

### 4.3 Crunch the data

Now compute δ<sup>13</sup>C, δ<sup>18</sup>Ο, and raw Δ<sub>47</sub>, Δ<sub>48</sub>, Δ<sub>49</sub> values. Note that δ<sup>18</sup>Ο is the CO<sub>2</sub> composition. The user is responsible for any acid fractionation correction.

```python
foo.crunch()

r = foo[0]
for k in r:
    print(f'r["{k}"] = {r[k]}')

# output:
# r["UID"] = A01
# r["Session"] = Session1
# r["Sample"] = ETH-1
# r["d45"] = 5.79502
# r["d46"] = 11.62767
# r["d47"] = 16.89351
# r["d48"] = 24.56708
# r["d49"] = 0.79486
# r["d13Cwg_VPDB"] = -3.7555729459832765
# r["d18Owg_VSMOW"] = 25.1145492463934
# r["D17O"] = 0.0
# r["d13C_VPDB"] = 1.9948594073404546
# r["d18O_VSMOW"] = 37.03357105550355
# r["D47raw"] = -0.5746856128030498
# r["D48raw"] = 1.1496833191546596
# r["D49raw"] = -27.690248970251407
```

### 4.4 Oxygen-17 correction parameters

Note that this crunching step uses the IUPAC oxygen-17 correction parameters, as recommended by [Daëron et al. (2016)](https://dx.doi.org/10.1016/j.chemgeo.2016.08.014) and  [Schauer et al. (2016)](https://dx.doi.org/10.1002/rcm.7743):

```python
R13_VPDB = 0.01118  # (Chang & Li, 1990)
R18_VSMOW = 0.0020052  # (Baertschi, 1976)
lambda_17 = 0.528  # (Barkan & Luz, 2005)
R17_VSMOW = 0.00038475  # (Assonov & Brenninkmeijer, 2003, rescaled to R13_VPDB)
R18_VPDB = R18_VSMOW * 1.03092
R17_VPDB = R17_VSMOW * 1.03092 ** lambda_17
```

To use different numerical values for these parameters, change them before performing `foo.crunch()`:

```python
# to change the lambda value to 0.5164,
# leaving the other parameters unchanged:
foo.lambda_17 = 0.5164
```

### 4.5 Reference frame

The nominal Δ<sub>47</sub> values assigned to the anchor samples are defined in `foo.Nominal_D47`, which may be redefined arbitrarily:

```python
print(foo.Nominal_D47) # default values from Bernasconi et al. (2018)
# output:
# {'ETH-1': 0.258, 'ETH-2': 0.256, 'ETH-3': 0.691}

foo.Nominal_D47 = {
    "Foo-1": 0.232,
    "Foo-2": 0.289,
    "Foo-3": 0.455,
    "Foo-4": 0.704,
    }

print(foo.Nominal_D47)
# output:
# {'Foo-1': 0.232, 'Foo-2': 0.289, 'Foo-3': 0.455, 'Foo-4': 0.704}
```

### 4.6 Standardization (`pooled`)

### 4.6.1 Default method (`pooled`)

The default standardization approach computes the best-fit standardization parameters (a,b,c) for each session, along with the best-fit Δ<sub>47</sub> values of unknown samples, using a pooled regression model taking into account the relative mapping of all samples (anchors and unknowns) in (δ<sub>47</sub>, Δ<sub>47</sub>) space.

```python
foo.standardize()
```

The following text is output:

```
--------------------------------  -----------
N samples (anchors + unknowns)      5 (3 + 2)
N analyses (anchors + unknowns)   20 (12 + 8)
Repeatability of δ13C_VPDB           13.8 ppm
Repeatability of δ18O_VSMOW          41.9 ppm
Repeatability of Δ47 (anchors)       10.7 ppm
Repeatability of Δ47 (unknowns)       3.4 ppm
Repeatability of Δ47 (all)            8.6 ppm
Model degrees of freedom                   12
Student's 95% t-factor                   2.18
--------------------------------  -----------

--------  --  --  -----------  ------------  ------  ------  ------  -------------  -------------  --------------
Session   Na  Nu  d13Cwg_VPDB  d18Owg_VSMOW  r_d13C  r_d18O   r_D47         a ± SE   1e3 x b ± SE          c ± SE
--------  --  --  -----------  ------------  ------  ------  ------  -------------  -------------  --------------
Session1   6   4       -3.756        25.115  0.0035  0.0415  0.0066  0.838 ± 0.016  3.340 ± 0.247  -0.859 ± 0.007
Session2   6   4       -3.743        25.117  0.0174  0.0490  0.0119  0.815 ± 0.015  4.601 ± 0.246  -0.847 ± 0.007
--------  --  --  -----------  ------------  ------  ------  ------  -------------  -------------  --------------


-------  -  ---------  ----------  ------  ------  --------  ------  --------
Sample   N  d13C_VPDB  d18O_VSMOW     D47      SE    95% CL      SD  p_Levene
-------  -  ---------  ----------  ------  ------  --------  ------  --------
ETH-1    4       2.00       37.00  0.2580                    0.0096          
ETH-2    4     -10.03       20.18  0.2560                    0.0154          
ETH-3    4       1.71       37.45  0.6910                    0.0039          
IAEA-C1  4       2.46       36.88  0.3624  0.0061  ± 0.0133  0.0031     0.901
IAEA-C2  4      -8.04       30.18  0.7246  0.0082  ± 0.0178  0.0037     0.825
-------  -  ---------  ----------  ------  ------  --------  ------  --------
```

### 4.6.1 `D47data.sessions`

Under the hood, the normalization step does many things. It stores session information in `foo.sessions`:

```python
print([k for k in foo.sessions])
# output: ['Session1', 'Session2']

for k in foo.sessions['Session1']:
    if k == 'data':
        print(f"{k:>16}: [...] (too large to print)")
    else:
        print(f"{k:>16}: {foo.sessions['Session1'][k]}")
# output:
#             data: [...] (too large to print)
# scrambling_drift: False
#      slope_drift: False
#         wg_drift: False
#      d13Cwg_VPDB: -3.7555729459832765
#     d18Owg_VSMOW: 25.1145492463934
#               Na: 6
#               Nu: 4
#                a: 0.8381700022050721
#             SE_a: 0.015603758280720111
#                b: 0.0033401762623331823
#             SE_b: 0.00024740622188942793
#                c: -0.8586982120784981
#             SE_c: 0.006737855778815339
#               a2: 0.0
#               b2: 0.0
#               c2: 0.0
#      r_d13C_VPDB: 0.0035270933192414504
#     r_d18O_VSMOW: 0.04146499779427257
#            r_D47: 0.006638319347677248
```

each element of `foo.sessions` has the following attributes:

+ `data`: list of all the analyses in this session
+ `scrambling_drift`, `slope_drift`, `wg_drift`: whether parameters `a`, `b`,`c` are allowed to drift (change linearly with with time)
+ `d13Cwg_VPDB`, `d18Owg_VSMOW`: working gas composition
+ `Na`: number of anchor analyses in this session
+ `Nu`: number of unknown analyses in this session
+ `a`,`SE_a`: best-fit value and model SE of scrambling factor
+ `b`,`SE_b`: best-fit value and model SE of compositional slope
+ `c`,`SE_c`: best-fit value and model SE of working gas offset
+ `a2`,`b2`,`c2`: drift rates (per unit of `TimeTag`) of `a`,`b`, `c`. If `TimeTag` is one of the fields in the raw data, this will be used, otherwise `TimeTag` starts at 0 for each session and increases by 1 for each analysis, in the listed order (thus beware of datasets ordered by sample name).
+ `r_d13C_VPDB`, `r_d18O_VSMOW`, `r_D47`: repeatabilities for `d13C_VPDB`, `d18O_VSMOW`, `D47` in this session

### 4.6.2 `D47data.samples`, `D47data.anchors`, and `D47data.unknowns`

Additional information about the samples is stored in `foo.samples` (the same information can also be accessed via `foo.anchors` and `foo.unknowns`):

```python
print([k for k in foo.samples])
# output:
# ['ETH-1', 'ETH-2', 'ETH-3', 'IAEA-C1', 'IAEA-C2']

for k in foo.samples['IAEA-C1']:
    if k == 'data':
        print(f"{k:>12}: [...] (too large to print)")
    else:
        print(f"{k:>12}: {foo.samples['IAEA-C1'][k]}")
# output:
#         data: [...] (too large to print)
#            N: 4
#       SD_D47: 0.003120794222015294
#    d13C_VPDB: 2.4606390899379327
#   d18O_VSMOW: 36.87682448377142
#          D47: 0.36241877475632883
#       SE_D47: 0.006107113137661028
#     p_Levene: 0.9011524351870661
```

Each element of `foo.samples` has the following attributes:

+ `N`: total number of analyses in the whole data set
+ `SD_D47`: the sample SD of Δ<sub>47</sub> for this sample
+ `d13C_VPDB`, `d18O_VSMOW`: average δ<sup>13</sup>C, δ<sup>18</sup>Ο values for the analyte CO<sub>2</sub>.
+ `D47`, `SE_D47`: best-fit value and model SE for the Δ<sub>47</sub> of this sample
+ `p_Levene`: p-value for a [Levene's test](https://en.wikipedia.org/wiki/Levene%27s_test) of whether the observed Δ<sub>47</sub> variance for this sample is significantly larger than that for ETH-3 (to change the reference sample to compare with, e.g. to ETH-1: `foo.LEVENE_REF_SAMPLE = 'ETH-1'`).

### 4.6.3 `D47data.repeatability`

The overall analytical repeatabilities are now saved to `foo.repeatability`:

```python
for k in foo.repeatability:
    print(f"{k:>12}: {foo.repeatability[k]}")

# output:
#  r_d13C_VPDB: 0.013821704833171146
# r_d18O_VSMOW: 0.04191487414887982
#       r_D47a: 0.010690471302409636
#       r_D47u: 0.0034370447628642863
#        r_D47: 0.008561367687546161
```

+ `r_d13C_VPDB`: Analytical repeatability of δ<sup>13</sup>C for all samples
+ `r_d18O_VSMOW`: Analytical repeatability of δ<sup>18</sup>O for all samples (CO<sub>2</sub> values)
+ `r_D47a`: Analytical repeatability of Δ<sub>47</sub> for anchor samples only
+ `r_D47u`: Analytical repeatability of Δ<sub>47</sub> for unknown samples only
+ `r_D47`: Analytical repeatability of Δ<sub>47</sub> for all samples.

### 4.6.4 `D47data.result`

By default `foo.normalize()` uses the [`lmfit.Minimizer.leastsq()`](https://lmfit.github.io/lmfit-py/fitting.html#lmfit.minimizer.Minimizer) method, which returns an instance of [`lmfit.MinimizerResult`](https://lmfit.github.io/lmfit-py/fitting.html#lmfit.minimizer.MinimizerResult). This `MinimizerResult`instance is stored in `foo.result`. A detailed report may be printed using `foo.report()`

```python
print(type(foo.normalization))
# output:
# <class 'lmfit.minimizer.MinimizerResult'>
```

