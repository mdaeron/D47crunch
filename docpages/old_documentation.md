## Usage

### 1. Create `D47data` object

Start by creating a `D47data` object which will store and process your data.

```python
import D47crunch 
foo = D47crunch.D47data(verbose = True)
```

The `verbose` keyword specifies whether to print out extra information when calling `D47data` methods.
This property may be arbitrarily changed using the `verbose` attribute of the resulting object:

```python
foo.verbose = False
```

Even before importing any analyses, our `D47data` object has properties which may be inspected and/or edited:

#### 1.1 Nominal δ<sup>13</sup>C<sub>VPDB</sub>, δ<sup>18</sup>O<sub>VPDB</sub>, and Δ<sub>47</sub> values of carbonate standards

`foo.Nominal_d13C_VPDB` and `foo.Nominal_d18O_VPDB` are dictionaries storing the δ<sup>13</sup>C<sub>VPDB</sub> and δ<sup>18</sup>O<sub>VPDB</sub> values of carbonate standards. You may freely edit these values and/or which standards to consider:

```python
print(foo.Nominal_d13C_VPDB)
# output: {'ETH-1': 2.02, 'ETH-2': -10.17, 'ETH-3': 1.71}

print(foo.Nominal_d18O_VPDB)
# {'ETH-1': -2.19, 'ETH-2': -18.69, 'ETH-3': -1.78}

foo.Nominal_d13C_VPDB['ETH-4'] = -10.20
foo.Nominal_d18O_VPDB['ETH-4'] = -18.81

print(foo.Nominal_d13C_VPDB)
# output: {'ETH-1': 2.02, 'ETH-2': -10.17, 'ETH-3': 1.71, 'ETH-4': -10.2}

print(foo.Nominal_d18O_VPDB)
# {'ETH-1': -2.19, 'ETH-2': -18.69, 'ETH-3': -1.78, 'ETH-4': -18.81}
```

`foo.Nominal_D47` is another dictionary, storing the absolute Δ<sub>47</sub> values of the standards used to anchor your measurements to an absolute Δ<sub>47</sub> reference frame.
As above, you may feely edit these values:

```python
print(foo.Nominal_D47)
# output: {'ETH-1': 0.258, 'ETH-2': 0.256, 'ETH-3': 0.691}

foo.Nominal_D47['ETH-4'] = 0.507

print(foo.Nominal_D47)
# output: {'ETH-1': 0.258, 'ETH-2': 0.256, 'ETH-3': 0.691, 'ETH-4': 0.507}
```

#### 1.2 Oxygen-17 correction parameters

The oxygen-17 correction parameters used by `D47data.crunch()` (see below) are specified by `foo.R13_VPDB`, `foo.R17_VSMOW`, `foo.R18_VSMOW` and `foo.lambda_17`. Default values correspond to the IUPAC values as recommended by [Daëron et al. (2016)] and  [Schauer et al. (2016)].

[Daëron et al. (2016)]: https://dx.doi.org/10.1016/j.chemgeo.2016.08.014
[Schauer et al. (2016)]: https://dx.doi.org/10.1002/rcm.7743

```python
print(foo.R13_VPDB)  # -> 0.01118    (Chang & Li, 1990)
print(foo.R17_VSMOW) # -> 0.00038475 (Assonov & Brenninkmeijer, 2003, rescaled to R13_VPDB)
print(foo.R18_VSMOW) # -> 0.0020052  (Baertschi, 1976)
print(foo.lambda_17) # -> 0.528      (Barkan & Luz, 2005)
```

As above, the values for these parameters may be arbitrarily redefined:

```python
# to change the lambda value to 0.5164, leaving the other parameters unchanged:
foo.lambda_17 = 0.5164
```

#### 1.3 Default method for carbon-13 and oxygen-18 standardization

By default, bulk isotopic compositions are standardized using a “two-point” affine transformation (correcting for small offsets and stretching effects) based on the carbonate standards defined in `foo.Nominal_d13C_VPDB` and `foo.Nominal_d18O_VPDB` (see above).

Optionally, you may opt instead for a “single-point” standardization approach not correcting for strecthing effects, for instance if the cabonate standards in `foo.Nominal_d13C_VPDB` and `foo.Nominal_d18O_VPDB` cover only a small fraction of the full isotopic range of your measurements.

Finally, you may also opt to perform no _a posteriori_ standardization of bulk isotopic compositions, which implies that the quality of your final δ<sup>13</sup>C and δ<sup>18</sup>O values will depend strongly on the accuracy of your working gas composition and the linearity of your instrument.

Switching betwwen these three options can be achieved by setting `foo.d13C_STANDARDIZATION_METHOD` and `foo.d18O_STANDARDIZATION_METHOD` to `'2pt'`, `'1pt'`, and `'none'`, respectively. Note that you may later override this default behavior on a per-session basis.

#### 1.4 Oxygen-18 acid fractionation factor

`D47data` processing methods always return δ<sup>18</sup>O values of CO<sub>2</sub> analytes relative to VSMOW rather than carbonate δ<sup>18</sup>O<sub>VPDB</sub> values (which depend on sample mineralogies and acid reaction temperature). However, when using single-point or two-point δ<sup>18</sup>O standardization or when computing the bulk isotope composition of working gases based on carbonate standards (using `D47data.wg()`), it is necessary to specify the oxygen-18 fractionation factor associated with the phosphoric acid reaction, by setting the value of `foo.ALPHA_18O_ACID_REACTION`.

### 2. Import data

It's time to add some analyses to out `D47data` object.

Start with some raw data stored as CSV in a file named `rawdata.csv` (spaces after commas are optional). Each line corresponds to a single analysis.

The only required fields are a sample identifier (`Sample`), and the working-gas delta values `d45`, `d46`, `d47`. If no session information is provided, all analyses will be treated as belonging to a single analytical session. Alternatively, to group analyses into sessions, provide session identifiers in a `Session` field. If not specified by the user, a unique identifier (`UID`) will be assigned automatically to each analysis. Independently known oxygen-17 anomalies may be provided as `D17O` (in ‰ relative to VSMOW, with λ equal to `D47data.lambda_17`), and are assumed to be zero otherwise. Working-gas deltas `d48` and `d49` may also be provided, and are otherwise treated as `nan`.

Example `rawdata.csv` file:

```csv
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

Reading data from `rawdata.csv` can be done with `foo.read()`:

```python
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

At this stage, `foo` behaves like a `list` object. Yoy may slice it in the usual way (`foo[:10]` returns a list of the first 10 analyses) and use built-in methods in the expected way (e.g., `len(foo)` is equal to 20).

We can inspect the first record now stored in `foo`, corresponding to a single analysis:

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

#### 2.1 Sessions

After importing records from `rawdata.csv`, our `D47data` object now has a new dictionary attribute, `foo.sessions`:

```python
for session in foo.sessions:
	print(f"{session:>28}:")
	for k in foo.sessions[session]:
		if k == 'data':
			print(f"{k:>28}: [...] (too large to print)")
		else:
			print(f"{k:>28}: {foo.sessions[session][k]}")
	print()
# output:
#                     Session1:
#                         data: [...] (too large to print)
#             scrambling_drift: False
#                  slope_drift: False
#                     wg_drift: False
#  d13C_standardization_method: 2pt
#  d18O_standardization_method: 2pt
# 
#                     Session2:
#                         data: [...] (too large to print)
#             scrambling_drift: False
#                  slope_drift: False
#                     wg_drift: False
#  d13C_standardization_method: 2pt
#  d18O_standardization_method: 2pt
```

Each session in `foo.sessions` has the following attributes at this stage:

+ `data`: list of all the analyses in this session
+ `scrambling_drift`, `slope_drift`, `wg_drift`: whether parameters `a`, `b`,`c` of the Δ<sub>47</sub> standardization model are allowed to drift (change linearly with with time).
+ `d13C_standardization_method`, `d18O_standardization_method`: which method to use for this session.

You may arbitrarily edit the values of `d13C_standardization_method`, `d18O_standardization_method` for any session, which will affect the results of `foo.crunch()`.

Similarly, you may arbitrarily edit the values of `scrambling_drift`, `slope_drift`, `wg_drift` for any session, which will affect the results of `foo.standardize()`.

#### 2.2 Samples, anchors and unknowns

### 3. Working gas composition

There are two ways to define the isotpic composition of the working gas.

#### 3.1 Option 1: explicit definition

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

#### 3.2 Option 2: based on the known composition of a sample:

```python
# The 2 code lines below are the default settings. It is thus not
# necessary to include them unless you wish to use different values.

foo.SAMPLE_CONSTRAINING_WG_COMPOSITION = ('ETH-3', 1.71, -1.78)
foo.ALPHA_18O_ACID_REACTION = 1.00813 # (Kim et al., 2007), calcite at 90 °C

# Compute the WG composition for each session:
foo.wg()

```

### 4. Crunch the data

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

### 6. Standardization

#### 6.1 Default approach (<span style="text-transform:lowercase">`pooled`</span>)

The default standardization approach computes the best-fit standardization parameters (a,b,c) for each session, along with the best-fit Δ<sub>47</sub> values of unknown samples, using a pooled regression model taking into account the relative mapping of all samples (anchors and unknowns) in (δ<sub>47</sub>, Δ<sub>47</sub>) space.

```python
foo.standardize()
foo.table_of_sessions(verbose = True, save_to_file = False)
foo.table_of_samples(verbose = True, save_to_file = False)

```

The following text is output:

```csv
[table_of_sessions]
–––––––––––––––––––––––––––––––  –––––––––––
N samples (anchors + unknowns)     5 (3 + 2)
N analyses (anchors + unknowns)  20 (12 + 8)
Repeatability of δ13C_VPDB          13.8 ppm
Repeatability of δ18O_VSMOW         41.9 ppm
Repeatability of Δ47 (anchors)      13.1 ppm
Repeatability of Δ47 (unknowns)      3.4 ppm
Repeatability of Δ47 (all)           9.6 ppm
Model degrees of freedom                  12
Student's 95% t-factor                  2.18
Standardization method                pooled
–––––––––––––––––––––––––––––––  –––––––––––

[table_of_sessions]
––––––––  ––  ––  –––––––––––  ––––––––––––  ––––––  ––––––  ––––––  –––––––––––––  –––––––––––––  ––––––––––––––
Session   Na  Nu  d13Cwg_VPDB  d18Owg_VSMOW  r_d13C  r_d18O   r_D47         a ± SE   1e3 x b ± SE          c ± SE
––––––––  ––  ––  –––––––––––  ––––––––––––  ––––––  ––––––  ––––––  –––––––––––––  –––––––––––––  ––––––––––––––
Session1   6   4       -3.756        25.115  0.0035  0.0415  0.0066  0.838 ± 0.016  3.340 ± 0.247  -0.859 ± 0.007
Session2   6   4       -3.743        25.118  0.0174  0.0490  0.0119  0.815 ± 0.015  4.601 ± 0.246  -0.847 ± 0.007
––––––––  ––  ––  –––––––––––  ––––––––––––  ––––––  ––––––  ––––––  –––––––––––––  –––––––––––––  ––––––––––––––


[table_of_samples] 
–––––––  –  –––––––––  ––––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––
Sample   N  d13C_VPDB  d18O_VSMOW     D47      SE    95% CL      SD  p_Levene
–––––––  –  –––––––––  ––––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––
ETH-1    4       2.00       37.00  0.2580                    0.0096          
ETH-2    4     -10.03       20.18  0.2560                    0.0154          
ETH-3    4       1.71       37.45  0.6910                    0.0039          
IAEA-C1  4       2.46       36.88  0.3624  0.0061  ± 0.0133  0.0031     0.901
IAEA-C2  4      -8.04       30.19  0.7246  0.0082  ± 0.0178  0.0037     0.825
–––––––  –  –––––––––  ––––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––
```

#### 6.2 `D47data().sessions`

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
#            data: [...] (too large to print)
# scrambling_drift: False
#      slope_drift: False
#         wg_drift: False
#      d13Cwg_VPDB: -3.7555729339153743
#     d18Owg_VSMOW: 25.11497520475171
#               Na: 6
#               Nu: 4
#      r_d13C_VPDB: 0.0035270930676685897
#     r_d18O_VSMOW: 0.04146501520018946
#            r_D47: 0.006638319178058144
#               Np: 3
#                a: 0.8381700110925523
#             SE_a: 0.015603757788793743
#                b: 0.003340175397346955
#             SE_b: 0.0002474062198065805
#                c: -0.8586981978192628
#             SE_c: 0.006737855663518676
#               a2: 0.0
#            SE_a2: 0.0
#               b2: 0.0
#            SE_b2: 0.0
#               c2: 0.0
#            SE_c2: 0.0
#               CM: [...] (6x6 numpy.Array())
```

each element of `foo.sessions` has the following attributes:

+ `data`: list of all the analyses in this session
+ `scrambling_drift`, `slope_drift`, `wg_drift`: whether parameters `a`, `b`,`c` are allowed to drift (change linearly with with time)
+ `d13Cwg_VPDB`, `d18Owg_VSMOW`: working gas composition
+ `Na`: number of anchor analyses in this session
+ `Nu`: number of unknown analyses in this session
+ `r_d13C_VPDB`, `r_d18O_VSMOW`, `r_D47`: repeatabilities for `d13C_VPDB`, `d18O_VSMOW`, `D47` in this session
+ `a`,`SE_a`: best-fit value and model SE of scrambling factor
+ `b`,`SE_b`: best-fit value and model SE of compositional slope
+ `c`,`SE_c`: best-fit value and model SE of working gas offset
+ `a2`,`b2`,`c2`: drift rates (per unit of `TimeTag`) of `a`,`b`, `c`. If `TimeTag` is one of the fields in the raw data, this will be used, otherwise `TimeTag` starts at 0 for each session and increases by 1 for each analysis, in the listed order (thus beware of datasets ordered by sample name).
+ `CM`: the covariance matrix of (`a`, `b`, `c`, `a2`, `b2`, `c2`).

#### 6.3 `D47data().samples`, `D47data().anchors`, and `D47data().unknowns`

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
#       SD_D47: 0.0031207941052170305
#    d13C_VPDB: 2.460639104889639
#   d18O_VSMOW: 36.87725533010137
#     p_Levene: 0.901152441112675
#          D47: 0.3624187694150056
#       SE_D47: 0.00610711296513016
```

Each element of `foo.samples` has the following attributes:

+ `N`: total number of analyses in the whole data set
+ `SD_D47`: the sample SD of Δ<sub>47</sub> for this sample
+ `d13C_VPDB`, `d18O_VSMOW`: average δ<sup>13</sup>C, δ<sup>18</sup>Ο values for the analyte CO<sub>2</sub>.
+ `D47`, `SE_D47`: best-fit value and model SE for the Δ<sub>47</sub> of this sample
+ `p_Levene`: p-value for a [Levene test] of whether the observed Δ<sub>47</sub> variance for this sample is significantly larger than that for ETH-3 (to change the reference sample to compare with, e.g. to ETH-1: `foo.LEVENE_REF_SAMPLE = 'ETH-1'` before calling `foo.normalize()`).

[Levene test]: https://en.wikipedia.org/wiki/Levene%27s_test

#### 6.4 `D47data.()repeatability`

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

#### 6.5 `D47data.()result`

By default `foo.normalize()` uses the [`lmfit.Minimizer.leastsq()`](https://lmfit.github.io/lmfit-py/fitting.html#lmfit.minimizer.Minimizer) method, which returns an instance of [`lmfit.MinimizerResult`](https://lmfit.github.io/lmfit-py/fitting.html#lmfit.minimizer.MinimizerResult). This `MinimizerResult`instance is stored in `foo.result`. A detailed report may be printed using `foo.report()`

```python
print(type(foo.normalization))
# output:
# <class 'lmfit.minimizer.MinimizerResult'>
```

#### 6.6 Combining information from carbonate anchors and equilibrated gases

The `constraints` argument to `D47data.standardize()` in the pooled regression approach may be used to specify arbitrary constraints between regression model parameters. For instance, if a data set comprises two carbonate standards (`ETH-1` and `ETH-2`) and two gas standards (`HG-1000C` and `EG-25C`), it is possible to specify the Δ<sub>47</sub> difference between `HG-1000C` and `EG-25C` explicitly, essentially constraining the scrambling factor `a` based on the gas standards while constraining the other parameters based on `ETH-1` and `ETH-2`:

```python
from D47crunch import D47data, fCO2eqD47_Wang

rawdata = D47data()
rawdata.read('foo.csv') # foo.csv not provided in this example
rawdata.wg()
rawdata.crunch()
constr = {'D47_EG_25C': f'D47_HG_1000C + {fCO2eqD47_Wang(25)-fCO2eqD47_Wang(1000)}')
rawdata.standardize(constraints = constr)
rawdata.table_of_samples()

# outputs something like:
# 
# [table_of_samples] 
# ––––––––  –––  –––––––––  ––––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––
# Sample      N  d13C_VPDB  d18O_VSMOW     D47      SE    95% CL      SD  p_Levene
# ––––––––  –––  –––––––––  ––––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––
# ETH-1     205       2.03       37.03  0.2052                    0.0089          
# ETH-2     213     -10.17       19.88  0.2085                    0.0078          
# EG-25C    138     -18.43       40.35  0.9195  0.0000  ± 0.0000  0.0095     1.000
# HG-1000C  180      -7.95       26.50  0.0242  0.0007  ± 0.0014  0.0085     0.367
# ––––––––  –––  –––––––––  ––––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––
```

Note that in the example abobe, because `HG-1000C` was constrained as a function of `EG-25C`, the SE in its Δ<sub>47</sub> value is not reported. For now, it must instead be computed based on that of `EG-25C` (in this simple case, the two standard errors are identical).

#### 6.7 Legacy standardization approach (<span style="text-transform:lowercase">`indep_sessions`</span>)

Following a more traditional approach, `foo.standardize(method = 'indep_sessions')` computes the best-fit standardization parameters (a,b,c) for each session using independent regression models (one per session) only taking into account the anchor samples (samples defined in `foo.Nominal_D47`), then computes the Δ<sub>47</sub> value for each analysis and  the weighted average Δ<sub>47</sub> value for each unknown sample.

### 7. Viewing and saving the results

> under construction

***