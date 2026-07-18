'''
Standardization and analytical error propagation of Δ47 and Δ48 clumped-isotope measurements

Process and standardize carbonate and/or CO2 clumped-isotope analyses,
from low-level data out of a dual-inlet mass spectrometer to final, “absolute”
Δ47, Δ48 and Δ49 values with fully propagated analytical error estimates
([Daëron, 2021](https://doi.org/10.1029/2020GC009592)).

The **tutorial** section takes you through a series of simple steps to import/process data and print out the results.
The **how-to** section provides instructions applicable to various specific tasks.

.. include:: ../../docpages/tutorial.md
.. include:: ../../docpages/howto.md
.. include:: ../../docpages/cli.md

<h1>API Documentation</h1>
'''

from ._metadata import *

import os
import numpy as np
import typer
import warnings
import uncertainties
from typing_extensions import Annotated
from statistics import stdev
from scipy.stats import t as tstudent
# from scipy.stats import levene
from scipy.interpolate import interp1d
from numpy import linalg
from lmfit import Minimizer, Parameters, report_fit
from matplotlib import pyplot as ppl
from datetime import datetime as dt
from functools import wraps
from colorsys import hls_to_rgb
from matplotlib import rcParams
from typer import rich_utils

rich_utils.STYLE_HELPTEXT = ''

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Helvetica'
rcParams['font.size'] = 10
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'sans'
rcParams['mathtext.bf'] = 'sans:bold'
rcParams['mathtext.it'] = 'sans:italic'
rcParams['mathtext.cal'] = 'sans:italic'
rcParams['mathtext.default'] = 'rm'
rcParams['xtick.major.size'] = 4
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.size'] = 4
rcParams['ytick.major.width'] = 1
rcParams['axes.grid'] = False
rcParams['axes.linewidth'] = 1
rcParams['grid.linewidth'] = .75
rcParams['grid.linestyle'] = '-'
rcParams['grid.alpha'] = .15
rcParams['savefig.dpi'] = 150

Petersen_etal_CO2eqD47 = np.array([[-12, 1.147113572], [-11, 1.139961218], [-10, 1.132872856], [-9, 1.125847677], [-8, 1.118884889], [-7, 1.111983708], [-6, 1.105143366], [-5, 1.098363105], [-4, 1.091642182], [-3, 1.084979862], [-2, 1.078375423], [-1, 1.071828156], [0, 1.065337360], [1, 1.058902349], [2, 1.052522443], [3, 1.046196976], [4, 1.039925291], [5, 1.033706741], [6, 1.027540690], [7, 1.021426510], [8, 1.015363585], [9, 1.009351306], [10, 1.003389075], [11, 0.997476303], [12, 0.991612409], [13, 0.985796821], [14, 0.980028975], [15, 0.974308318], [16, 0.968634304], [17, 0.963006392], [18, 0.957424055], [19, 0.951886769], [20, 0.946394020], [21, 0.940945302], [22, 0.935540114], [23, 0.930177964], [24, 0.924858369], [25, 0.919580851], [26, 0.914344938], [27, 0.909150167], [28, 0.903996080], [29, 0.898882228], [30, 0.893808167], [31, 0.888773459], [32, 0.883777672], [33, 0.878820382], [34, 0.873901170], [35, 0.869019623], [36, 0.864175334], [37, 0.859367901], [38, 0.854596929], [39, 0.849862028], [40, 0.845162813], [41, 0.840498905], [42, 0.835869931], [43, 0.831275522], [44, 0.826715314], [45, 0.822188950], [46, 0.817696075], [47, 0.813236341], [48, 0.808809404], [49, 0.804414926], [50, 0.800052572], [51, 0.795722012], [52, 0.791422922], [53, 0.787154979], [54, 0.782917869], [55, 0.778711277], [56, 0.774534898], [57, 0.770388426], [58, 0.766271562], [59, 0.762184010], [60, 0.758125479], [61, 0.754095680], [62, 0.750094329], [63, 0.746121147], [64, 0.742175856], [65, 0.738258184], [66, 0.734367860], [67, 0.730504620], [68, 0.726668201], [69, 0.722858343], [70, 0.719074792], [71, 0.715317295], [72, 0.711585602], [73, 0.707879469], [74, 0.704198652], [75, 0.700542912], [76, 0.696912012], [77, 0.693305719], [78, 0.689723802], [79, 0.686166034], [80, 0.682632189], [81, 0.679122047], [82, 0.675635387], [83, 0.672171994], [84, 0.668731654], [85, 0.665314156], [86, 0.661919291], [87, 0.658546854], [88, 0.655196641], [89, 0.651868451], [90, 0.648562087], [91, 0.645277352], [92, 0.642014054], [93, 0.638771999], [94, 0.635551001], [95, 0.632350872], [96, 0.629171428], [97, 0.626012487], [98, 0.622873870], [99, 0.619755397], [100, 0.616656895], [102, 0.610519107], [104, 0.604459143], [106, 0.598475670], [108, 0.592567388], [110, 0.586733026], [112, 0.580971342], [114, 0.575281125], [116, 0.569661187], [118, 0.564110371], [120, 0.558627545], [122, 0.553211600], [124, 0.547861454], [126, 0.542576048], [128, 0.537354347], [130, 0.532195337], [132, 0.527098028], [134, 0.522061450], [136, 0.517084654], [138, 0.512166711], [140, 0.507306712], [142, 0.502503768], [144, 0.497757006], [146, 0.493065573], [148, 0.488428634], [150, 0.483845370], [152, 0.479314980], [154, 0.474836677], [156, 0.470409692], [158, 0.466033271], [160, 0.461706674], [162, 0.457429176], [164, 0.453200067], [166, 0.449018650], [168, 0.444884242], [170, 0.440796174], [172, 0.436753787], [174, 0.432756438], [176, 0.428803494], [178, 0.424894334], [180, 0.421028350], [182, 0.417204944], [184, 0.413423530], [186, 0.409683531], [188, 0.405984383], [190, 0.402325531], [192, 0.398706429], [194, 0.395126543], [196, 0.391585347], [198, 0.388082324], [200, 0.384616967], [202, 0.381188778], [204, 0.377797268], [206, 0.374441954], [208, 0.371122364], [210, 0.367838033], [212, 0.364588505], [214, 0.361373329], [216, 0.358192065], [218, 0.355044277], [220, 0.351929540], [222, 0.348847432], [224, 0.345797540], [226, 0.342779460], [228, 0.339792789], [230, 0.336837136], [232, 0.333912113], [234, 0.331017339], [236, 0.328152439], [238, 0.325317046], [240, 0.322510795], [242, 0.319733329], [244, 0.316984297], [246, 0.314263352], [248, 0.311570153], [250, 0.308904364], [252, 0.306265654], [254, 0.303653699], [256, 0.301068176], [258, 0.298508771], [260, 0.295975171], [262, 0.293467070], [264, 0.290984167], [266, 0.288526163], [268, 0.286092765], [270, 0.283683684], [272, 0.281298636], [274, 0.278937339], [276, 0.276599517], [278, 0.274284898], [280, 0.271993211], [282, 0.269724193], [284, 0.267477582], [286, 0.265253121], [288, 0.263050554], [290, 0.260869633], [292, 0.258710110], [294, 0.256571741], [296, 0.254454286], [298, 0.252357508], [300, 0.250281174], [302, 0.248225053], [304, 0.246188917], [306, 0.244172542], [308, 0.242175707], [310, 0.240198194], [312, 0.238239786], [314, 0.236300272], [316, 0.234379441], [318, 0.232477087], [320, 0.230593005], [322, 0.228726993], [324, 0.226878853], [326, 0.225048388], [328, 0.223235405], [330, 0.221439711], [332, 0.219661118], [334, 0.217899439], [336, 0.216154491], [338, 0.214426091], [340, 0.212714060], [342, 0.211018220], [344, 0.209338398], [346, 0.207674420], [348, 0.206026115], [350, 0.204393315], [355, 0.200378063], [360, 0.196456139], [365, 0.192625077], [370, 0.188882487], [375, 0.185226048], [380, 0.181653511], [385, 0.178162694], [390, 0.174751478], [395, 0.171417807], [400, 0.168159686], [405, 0.164975177], [410, 0.161862398], [415, 0.158819521], [420, 0.155844772], [425, 0.152936426], [430, 0.150092806], [435, 0.147312286], [440, 0.144593281], [445, 0.141934254], [450, 0.139333710], [455, 0.136790195], [460, 0.134302294], [465, 0.131868634], [470, 0.129487876], [475, 0.127158722], [480, 0.124879906], [485, 0.122650197], [490, 0.120468398], [495, 0.118333345], [500, 0.116243903], [505, 0.114198970], [510, 0.112197471], [515, 0.110238362], [520, 0.108320625], [525, 0.106443271], [530, 0.104605335], [535, 0.102805877], [540, 0.101043985], [545, 0.099318768], [550, 0.097629359], [555, 0.095974915], [560, 0.094354612], [565, 0.092767650], [570, 0.091213248], [575, 0.089690648], [580, 0.088199108], [585, 0.086737906], [590, 0.085306341], [595, 0.083903726], [600, 0.082529395], [605, 0.081182697], [610, 0.079862998], [615, 0.078569680], [620, 0.077302141], [625, 0.076059794], [630, 0.074842066], [635, 0.073648400], [640, 0.072478251], [645, 0.071331090], [650, 0.070206399], [655, 0.069103674], [660, 0.068022424], [665, 0.066962168], [670, 0.065922439], [675, 0.064902780], [680, 0.063902748], [685, 0.062921909], [690, 0.061959837], [695, 0.061016122], [700, 0.060090360], [705, 0.059182157], [710, 0.058291131], [715, 0.057416907], [720, 0.056559120], [725, 0.055717414], [730, 0.054891440], [735, 0.054080860], [740, 0.053285343], [745, 0.052504565], [750, 0.051738210], [755, 0.050985971], [760, 0.050247546], [765, 0.049522643], [770, 0.048810974], [775, 0.048112260], [780, 0.047426227], [785, 0.046752609], [790, 0.046091145], [795, 0.045441581], [800, 0.044803668], [805, 0.044177164], [810, 0.043561831], [815, 0.042957438], [820, 0.042363759], [825, 0.041780573], [830, 0.041207664], [835, 0.040644822], [840, 0.040091839], [845, 0.039548516], [850, 0.039014654], [855, 0.038490063], [860, 0.037974554], [865, 0.037467944], [870, 0.036970054], [875, 0.036480707], [880, 0.035999734], [885, 0.035526965], [890, 0.035062238], [895, 0.034605393], [900, 0.034156272], [905, 0.033714724], [910, 0.033280598], [915, 0.032853749], [920, 0.032434032], [925, 0.032021309], [930, 0.031615443], [935, 0.031216300], [940, 0.030823749], [945, 0.030437663], [950, 0.030057915], [955, 0.029684385], [960, 0.029316951], [965, 0.028955498], [970, 0.028599910], [975, 0.028250075], [980, 0.027905884], [985, 0.027567229], [990, 0.027234006], [995, 0.026906112], [1000, 0.026583445], [1005, 0.026265908], [1010, 0.025953405], [1015, 0.025645841], [1020, 0.025343124], [1025, 0.025045163], [1030, 0.024751871], [1035, 0.024463160], [1040, 0.024178947], [1045, 0.023899147], [1050, 0.023623680], [1055, 0.023352467], [1060, 0.023085429], [1065, 0.022822491], [1070, 0.022563577], [1075, 0.022308615], [1080, 0.022057533], [1085, 0.021810260], [1090, 0.021566729], [1095, 0.021326872], [1100, 0.021090622]])
_fCO2eqD47_Petersen = interp1d(Petersen_etal_CO2eqD47[:,0], Petersen_etal_CO2eqD47[:,1])
def fCO2eqD47_Petersen(T):
	'''
	CO2 equilibrium Δ47 value as a function of T (in degrees C)
	according to [Petersen et al. (2019)](https://doi.org/10.1029/2018GC008127).

	'''
	return float(_fCO2eqD47_Petersen(T))


Wang_etal_CO2eqD47 = np.array([[-83., 1.8954], [-73., 1.7530], [-63., 1.6261], [-53., 1.5126], [-43., 1.4104], [-33., 1.3182], [-23., 1.2345], [-13., 1.1584], [-3., 1.0888], [7., 1.0251], [17., 0.9665], [27., 0.9125], [37., 0.8626], [47., 0.8164], [57., 0.7734], [67., 0.7334], [87., 0.6612], [97., 0.6286], [107., 0.5980], [117., 0.5693], [127., 0.5423], [137., 0.5169], [147., 0.4930], [157., 0.4704], [167., 0.4491], [177., 0.4289], [187., 0.4098], [197., 0.3918], [207., 0.3747], [217., 0.3585], [227., 0.3431], [237., 0.3285], [247., 0.3147], [257., 0.3015], [267., 0.2890], [277., 0.2771], [287., 0.2657], [297., 0.2550], [307., 0.2447], [317., 0.2349], [327., 0.2256], [337., 0.2167], [347., 0.2083], [357., 0.2002], [367., 0.1925], [377., 0.1851], [387., 0.1781], [397., 0.1714], [407., 0.1650], [417., 0.1589], [427., 0.1530], [437., 0.1474], [447., 0.1421], [457., 0.1370], [467., 0.1321], [477., 0.1274], [487., 0.1229], [497., 0.1186], [507., 0.1145], [517., 0.1105], [527., 0.1068], [537., 0.1031], [547., 0.0997], [557., 0.0963], [567., 0.0931], [577., 0.0901], [587., 0.0871], [597., 0.0843], [607., 0.0816], [617., 0.0790], [627., 0.0765], [637., 0.0741], [647., 0.0718], [657., 0.0695], [667., 0.0674], [677., 0.0654], [687., 0.0634], [697., 0.0615], [707., 0.0597], [717., 0.0579], [727., 0.0562], [737., 0.0546], [747., 0.0530], [757., 0.0515], [767., 0.0500], [777., 0.0486], [787., 0.0472], [797., 0.0459], [807., 0.0447], [817., 0.0435], [827., 0.0423], [837., 0.0411], [847., 0.0400], [857., 0.0390], [867., 0.0380], [877., 0.0370], [887., 0.0360], [897., 0.0351], [907., 0.0342], [917., 0.0333], [927., 0.0325], [937., 0.0317], [947., 0.0309], [957., 0.0302], [967., 0.0294], [977., 0.0287], [987., 0.0281], [997., 0.0274], [1007., 0.0268], [1017., 0.0261], [1027., 0.0255], [1037., 0.0249], [1047., 0.0244], [1057., 0.0238], [1067., 0.0233], [1077., 0.0228], [1087., 0.0223], [1097., 0.0218]])
_fCO2eqD47_Wang = interp1d(Wang_etal_CO2eqD47[:,0] - 0.15, Wang_etal_CO2eqD47[:,1])
def fCO2eqD47_Wang(T):
	'''
	CO2 equilibrium Δ47 value as a function of `T` (in degrees C)
	according to [Wang et al. (2004)](https://doi.org/10.1016/j.gca.2004.05.039)
	(supplementary data of [Dennis et al., 2011](https://doi.org/10.1016/j.gca.2011.09.025)).
	'''
	return float(_fCO2eqD47_Wang(T))


def correlated_sum(X, C, w = None):
	'''
	Compute covariance-aware linear combinations

	**Parameters**

	+ `X`: list or 1-D array of values to sum
	+ `C`: covariance matrix for the elements of `X`
	+ `w`: list or 1-D array of weights to apply to the elements of `X`
	       (all equal to 1 by default)

	Return the sum (and its SE) of the elements of `X`, with optional weights equal
	to the elements of `w`, accounting for covariances between the elements of `X`.
	'''
	if w is None:
		w = [1 for x in X]
	return np.dot(w,X), (np.dot(w,np.dot(C,w)))**.5


def make_csv(x, hsep = ',', vsep = '\n'):
	'''
	Formats a list of lists of strings as a CSV

	**Parameters**

	+ `x`: the list of lists of strings to format
	+ `hsep`: the field separator (`,` by default)
	+ `vsep`: the line-ending convention to use (`\\n` by default)

	**Example**

	```py
	print(make_csv([['a', 'b', 'c'], ['d', 'e', 'f']]))
	```

	outputs:

	```py
	a,b,c
	d,e,f
	```
	'''
	return vsep.join([hsep.join(l) for l in x])


def pf(txt):
	'''
	Modify string `txt` to follow `lmfit.Parameter()` naming rules.
	'''
	return txt.replace('-','_').replace('.','_').replace(' ','_')


def smart_type(x):
	'''
	Tries to convert string `x` to a float if it includes a decimal point, or
	to an integer if it does not. If both attempts fail, return the original
	string unchanged.
	'''
	try:
		y = float(x)
	except ValueError:
		return x
	if '.' not in x:
		return int(y)
	return y

class _Defaults():
	def __init__(self):
		pass

D47crunch_defaults = _Defaults()
D47crunch_defaults.PRETTY_TABLE_VSEP = '—'

def pretty_table(x, header = 1, hsep = '  ', vsep = None, align = '<'):
	'''
	Reads a list of lists of strings and outputs an ascii table

	**Parameters**

	+ `x`: a list of lists of strings
	+ `header`: the number of lines to treat as header lines
	+ `hsep`: the horizontal separator between columns
	+ `vsep`: the character to use as vertical separator
	+ `align`: string of left (`<`) or right (`>`) alignment characters.

	**Example**

	```py
	print(pretty_table([
		['A', 'B', 'C'],
		['1', '1.9999', 'foo'],
		['10', 'x', 'bar'],
	]))
	```
	yields:
	```
	——  ——————  ———
	A        B    C
	——  ——————  ———
	1   1.9999  foo
	10       x  bar
	——  ——————  ———
	```

	To change the default `vsep` globally, redefine `D47crunch_defaults.PRETTY_TABLE_VSEP`:

	```py
	D47crunch_defaults.PRETTY_TABLE_VSEP = '='
	print(pretty_table([
		['A', 'B', 'C'],
		['1', '1.9999', 'foo'],
		['10', 'x', 'bar'],
	]))
	```
	yields:
	```
	==  ======  ===
	A        B    C
	==  ======  ===
	1   1.9999  foo
	10       x  bar
	==  ======  ===
	```
	'''

	if vsep is None:
		vsep = D47crunch_defaults.PRETTY_TABLE_VSEP

	txt = []
	widths = [np.max([len(e) for e in c]) for c in zip(*x)]

	if len(widths) > len(align):
		align += '>' * (len(widths)-len(align))
	sepline = hsep.join([vsep*w for w in widths])
	txt += [sepline]
	for k,l in enumerate(x):
		if k and k == header:
			txt += [sepline]
		txt += [hsep.join([f'{e:{a}{w}}' for e, w, a in zip(l, widths, align)])]
	txt += [sepline]
	txt += ['']
	return '\n'.join(txt)


def transpose_table(x):
	'''
	Transpose a list if lists

	**Parameters**

	+ `x`: a list of lists

	**Example**

	```py
	x = [[1, 2], [3, 4]]
	print(transpose_table(x)) # yields: [[1, 3], [2, 4]]
	```
	'''
	return [[e for e in c] for c in zip(*x)]


def w_avg(X, sX) :
	'''
	Compute variance-weighted average

	Returns the value and SE of the weighted average of the elements of `X`,
	with relative weights equal to their inverse variances (`1/sX**2`).

	**Parameters**

	+ `X`: array-like of elements to average
	+ `sX`: array-like of the corresponding SE values

	**Tip**

	If `X` and `sX` are initially arranged as a list of `(x, sx)` doublets,
	they may be rearranged using `zip()`:

	```python
	foo = [(0, 1), (1, 0.5), (2, 0.5)]
	print(w_avg(*zip(*foo))) # yields: (1.3333333333333333, 0.3333333333333333)
	```
	'''
	X = [ x for x in X ]
	sX = [ sx for sx in sX ]
	W = [ sx**-2 for sx in sX ]
	W = [ w/sum(W) for w in W ]
	Xavg = sum([ w*x for w,x in zip(W,X) ])
	sXavg = sum([ w**2*sx**2 for w,sx in zip(W,sX) ])**.5
	return Xavg, sXavg


def read_csv(filename, sep = ''):
	'''
	Read contents of `filename` in csv format and return a list of dictionaries.

	In the csv string, spaces before and after field separators (`','` by default)
	are optional.

	**Parameters**

	+ `filename`: the csv file to read
	+ `sep`: csv separator delimiting the fields. By default, use `,`, `;`, or `\t`,
	whichever appers most often in the contents of `filename`.
	'''
	with open(filename) as fid:
		txt = fid.read()

	if sep == '':
		sep = sorted(',;\t', key = lambda x: - txt.count(x))[0]
	txt = [[x.strip() for x in l.split(sep)] for l in txt.splitlines() if l.strip()]
	return [{k: smart_type(v) for k,v in zip(txt[0], l) if v} for l in txt[1:]]


def simulate_single_analysis(
	sample = 'MYSAMPLE',
	d13Cwg_VPDB = -4., d18Owg_VSMOW = 26.,
	d13C_VPDB = None, d18O_VPDB = None,
	D47 = None, D48 = None, D49 = 0., D17O = 0.,
	a47 = 1., b47 = 0., c47 = -0.9,
	a48 = 1., b48 = 0., c48 = -0.45,
	D47_RMs = None,
	D48_RMs = None,
	Nominal_d13C_VPDB = None,
	Nominal_d18O_VPDB = None,
	ALPHA_18O_ACID_REACTION = None,
	R13_VPDB = None,
	R17_VSMOW = None,
	R18_VSMOW = None,
	LAMBDA_17 = None,
	R18_VPDB = None,
	):
	'''
	Compute working-gas delta values for a single analysis, assuming a stochastic working
	gas and a “perfect” measurement (i.e. raw Δ values are identical to absolute values).

	**Parameters**

	+ `sample`: sample name
	+ `d13Cwg_VPDB`, `d18Owg_VSMOW`: bulk composition of the working gas
		(respectively –4 and +26 ‰ by default)
	+ `d13C_VPDB`, `d18O_VPDB`: bulk composition of the carbonate sample
	+ `D47`, `D48`, `D49`, `D17O`: clumped-isotope and oxygen-17 anomalies
		of the carbonate sample
	+ `D47_RMs`, `D48_RMs`: where to lookup Δ47 and
		Δ48 values if `D47` or `D48` are not specified
	+ `Nominal_d13C_VPDB`, `Nominal_d18O_VPDB`: where to lookup δ13C and
		δ18O values if `d13C_VPDB` or `d18O_VPDB` are not specified
	+ `ALPHA_18O_ACID_REACTION`: 18O/16O acid fractionation factor
	+ `R13_VPDB`, `R17_VSMOW`, `R18_VSMOW`, `LAMBDA_17`, `R18_VPDB`: oxygen-17
		correction parameters (by default equal to the `D4xdata` default values)

	Returns a dictionary with fields
	`['Sample', 'D17O', 'd13Cwg_VPDB', 'd18Owg_VSMOW', 'd45', 'd46', 'd47', 'd48', 'd49']`.
	'''

	if Nominal_d13C_VPDB is None:
		Nominal_d13C_VPDB = D4xdata().Nominal_d13C_VPDB

	if Nominal_d18O_VPDB is None:
		Nominal_d18O_VPDB = D4xdata().Nominal_d18O_VPDB

	if ALPHA_18O_ACID_REACTION is None:
		ALPHA_18O_ACID_REACTION = D4xdata().ALPHA_18O_ACID_REACTION

	if R13_VPDB is None:
		R13_VPDB = D4xdata().R13_VPDB

	if R17_VSMOW is None:
		R17_VSMOW = D4xdata().R17_VSMOW

	if R18_VSMOW is None:
		R18_VSMOW = D4xdata().R18_VSMOW

	if LAMBDA_17 is None:
		LAMBDA_17 = D4xdata().LAMBDA_17

	if R18_VPDB is None:
		R18_VPDB = D4xdata().R18_VPDB

	R17_VPDB = R17_VSMOW * (R18_VPDB / R18_VSMOW) ** LAMBDA_17

	if D47_RMs is None:
		D47_RMs = D47data().RMs

	if D48_RMs is None:
		D48_RMs = D48data().RMs

	if d13C_VPDB is None:
		if sample in Nominal_d13C_VPDB:
			d13C_VPDB = Nominal_d13C_VPDB[sample]
		else:
			raise KeyError(f"Sample {sample} is missing d13C_VPDB value, and it is not defined in Nominal_d13C_VPDB.")

	if d18O_VPDB is None:
		if sample in Nominal_d18O_VPDB:
			d18O_VPDB = Nominal_d18O_VPDB[sample]
		else:
			raise KeyError(f"Sample {sample} is missing d18O_VPDB value, and it is not defined in Nominal_d18O_VPDB.")

	if D47 is None:
		if sample in D47_RMs:
			D47 = D47_RMs[sample]
		else:
			raise KeyError(f"Sample {sample} is missing D47 value, and it is not defined in D47_RMs.")

	if D48 is None:
		if sample in D48_RMs:
			D48 = D48_RMs[sample]
		else:
			raise KeyError(f"Sample {sample} is missing D48 value, and it is not defined in D48_RMs.")

	X = D4xdata()
	X.R13_VPDB = R13_VPDB
	X.R17_VSMOW = R17_VSMOW
	X.R18_VSMOW = R18_VSMOW
	X.LAMBDA_17 = LAMBDA_17
	X.R18_VPDB = R18_VPDB
	X.R17_VPDB = R17_VSMOW * (R18_VPDB / R18_VSMOW)**LAMBDA_17

	R45wg, R46wg, R47wg, R48wg, R49wg = X.compute_isobar_ratios(
		R13 = R13_VPDB * (1 + d13Cwg_VPDB/1000),
		R18 = R18_VSMOW * (1 + d18Owg_VSMOW/1000),
		)
	R45, R46, R47, R48, R49 = X.compute_isobar_ratios(
		R13 = R13_VPDB * (1 + d13C_VPDB/1000),
		R18 = R18_VPDB * (1 + d18O_VPDB/1000) * ALPHA_18O_ACID_REACTION,
		D17O=D17O, D47=D47, D48=D48, D49=D49,
		)
	R45stoch, R46stoch, R47stoch, R48stoch, R49stoch = X.compute_isobar_ratios(
		R13 = R13_VPDB * (1 + d13C_VPDB/1000),
		R18 = R18_VPDB * (1 + d18O_VPDB/1000) * ALPHA_18O_ACID_REACTION,
		D17O=D17O,
		)

	d45 = 1000 * (R45/R45wg - 1)
	d46 = 1000 * (R46/R46wg - 1)
	d47 = 1000 * (R47/R47wg - 1)
	d48 = 1000 * (R48/R48wg - 1)
	d49 = 1000 * (R49/R49wg - 1)

	for k in range(3): # dumb iteration to adjust for small changes in d47
		R47raw = (1 + (a47 * D47 + b47 * d47 + c47)/1000) * R47stoch
		R48raw = (1 + (a48 * D48 + b48 * d48 + c48)/1000) * R48stoch
		d47 = 1000 * (R47raw/R47wg - 1)
		d48 = 1000 * (R48raw/R48wg - 1)

	return dict(
		Sample = sample,
		D17O = D17O,
		d13Cwg_VPDB = d13Cwg_VPDB,
		d18Owg_VSMOW = d18Owg_VSMOW,
		d45 = d45,
		d46 = d46,
		d47 = d47,
		d48 = d48,
		d49 = d49,
		)


def virtual_data(
	samples = [],
	a47 = 1., b47 = 0., c47 = -0.9,
	a48 = 1., b48 = 0., c48 = -0.45,
	rd45 = 0.020, rd46 = 0.060,
	rD47 = 0.015, rD48 = 0.045,
	d13Cwg_VPDB = None, d18Owg_VSMOW = None,
	session = None,
	D47_RMs = None, D48_RMs = None,
	Nominal_d13C_VPDB = None, Nominal_d18O_VPDB = None,
	ALPHA_18O_ACID_REACTION = None,
	R13_VPDB = None,
	R17_VSMOW = None,
	R18_VSMOW = None,
	LAMBDA_17 = None,
	R18_VPDB = None,
	seed = 0,
	shuffle = True,
	):
	'''
	Return list with simulated analyses from a single session.

	**Parameters**

	+ `samples`: a list of entries; each entry is a dictionary with the following fields:
	    * `Sample`: the name of the sample
	    * `d13C_VPDB`, `d18O_VPDB`: bulk composition of the carbonate sample
	    * `D47`, `D48`, `D49`, `D17O` (all optional): clumped-isotope and oxygen-17 anomalies of the carbonate sample
	    * `N`: how many analyses to generate for this sample
	+ `a47`: scrambling factor for Δ47
	+ `b47`: compositional nonlinearity for Δ47
	+ `c47`: working gas offset for Δ47
	+ `a48`: scrambling factor for Δ48
	+ `b48`: compositional nonlinearity for Δ48
	+ `c48`: working gas offset for Δ48
	+ `rd45`: analytical repeatability of δ45
	+ `rd46`: analytical repeatability of δ46
	+ `rD47`: analytical repeatability of Δ47
	+ `rD48`: analytical repeatability of Δ48
	+ `d13Cwg_VPDB`, `d18Owg_VSMOW`: bulk composition of the working gas
		(by default equal to the `simulate_single_analysis` default values)
	+ `session`: name of the session (no name by default)
	+ `D47_RMs`, `D48_RMs`: where to lookup Δ47 and Δ48 values
		if `D47` or `D48` are not specified (by default equal to the `simulate_single_analysis` defaults)
	+ `Nominal_d13C_VPDB`, `Nominal_d18O_VPDB`: where to lookup δ13C and
		δ18O values if `d13C_VPDB` or `d18O_VPDB` are not specified
		(by default equal to the `simulate_single_analysis` defaults)
	+ `ALPHA_18O_ACID_REACTION`: 18O/16O acid fractionation factor
		(by default equal to the `simulate_single_analysis` defaults)
	+ `R13_VPDB`, `R17_VSMOW`, `R18_VSMOW`, `LAMBDA_17`, `R18_VPDB`: oxygen-17
		correction parameters (by default equal to the `simulate_single_analysis` default)
	+ `seed`: explicitly set to a non-zero value to achieve random but repeatable simulations
	+ `shuffle`: randomly reorder the sequence of analyses


	Here is an example of using this method to generate an arbitrary combination of
	anchors and unknowns for a bunch of sessions:

	```py
	.. include:: ../../code_examples/virtual_data/example.py
	```

	This should output something like:

	```
	.. include:: ../../code_examples/virtual_data/output.txt
	```
	'''

	kwargs = locals().copy()

	from numpy import random as nprandom
	if seed:
		nprandom.seed(seed)
		rng = nprandom.default_rng(seed)
	else:
		rng = nprandom.default_rng()

	N = sum([s['N'] for s in samples])
	errors45 = rng.normal(loc = 0, scale = 1, size = N) # generate random measurement errors
	errors45 *= rd45 / stdev(errors45) # scale errors to rd45
	errors46 = rng.normal(loc = 0, scale = 1, size = N) # generate random measurement errors
	errors46 *= rd46 / stdev(errors46) # scale errors to rd46
	errors47 = rng.normal(loc = 0, scale = 1, size = N) # generate random measurement errors
	errors47 *= rD47 / stdev(errors47) # scale errors to rD47
	errors48 = rng.normal(loc = 0, scale = 1, size = N) # generate random measurement errors
	errors48 *= rD48 / stdev(errors48) # scale errors to rD48

	k = 0
	out = []
	for s in samples:
		kw = {}
		kw['sample'] = s['Sample']
		kw = {
			**kw,
			**{var: kwargs[var]
				for var in [
					'd13Cwg_VPDB', 'd18Owg_VSMOW', 'ALPHA_18O_ACID_REACTION',
					'D47_RMs', 'D48_RMs', 'Nominal_d13C_VPDB', 'Nominal_d18O_VPDB',
					'R13_VPDB', 'R17_VSMOW', 'R18_VSMOW', 'LAMBDA_17', 'R18_VPDB',
					'a47', 'b47', 'c47', 'a48', 'b48', 'c48',
					]
				if kwargs[var] is not None},
			**{var: s[var]
				for var in ['d13C_VPDB', 'd18O_VPDB', 'D47', 'D48', 'D49', 'D17O']
				if var in s},
			}

		sN = s['N']
		while sN:
			out.append(simulate_single_analysis(**kw))
			out[-1]['d45'] += errors45[k]
			out[-1]['d46'] += errors46[k]
			out[-1]['d47'] += (errors45[k] + errors46[k] + errors47[k]) * a47
			out[-1]['d48'] += (2*errors46[k] + errors48[k]) * a48
			sN -= 1
			k += 1

		if session is not None:
			for r in out:
				r['Session'] = session

		if shuffle:
			nprandom.shuffle(out)

	return out

def table_of_samples(
	data47 = None,
	data48 = None,
	dir = 'output',
	filename = None,
	save_to_file = True,
	print_out = True,
	output = None,
	):
	'''
	Print out, save to disk and/or return a combined table of samples
	for a pair of `D47data` and `D48data` objects.

	**Parameters**

	+ `data47`: `D47data` instance
	+ `data48`: `D48data` instance
	+ `dir`: the directory in which to save the table
	+ `filename`: the name to the csv file to write to
	+ `save_to_file`: whether to save the table to disk
	+ `print_out`: whether to print out the table
	+ `output`: if set to `'pretty'`: return a pretty text table (see `pretty_table()`);
		if set to `'raw'`: return a list of list of strings
		(e.g., `[['header1', 'header2'], ['0.1', '0.2']]`)
	'''
	if data47 is None:
		if data48 is None:
			raise TypeError("Arguments must include at least one D47data() or D48data() instance.")
		else:
			return data48.table_of_samples(
				dir = dir,
				filename = filename,
				save_to_file = save_to_file,
				print_out = print_out,
				output = output
				)
	else:
		if data48 is None:
			return data47.table_of_samples(
				dir = dir,
				filename = filename,
				save_to_file = save_to_file,
				print_out = print_out,
				output = output
				)
		else:
			samples = (
				sorted([a for a in data47.anchors if a in data48.anchors])
				+ sorted([a for a in data47.anchors if a not in data48.anchors])
				+ sorted([a for a in data48.anchors if a not in data47.anchors])
				+ sorted([a for a in data47.unknowns if a in data48.unknowns])
			)

			out47 = data47.table_of_samples(save_to_file = False, print_out = False, output = 'raw')
			out48 = data48.table_of_samples(save_to_file = False, print_out = False, output = 'raw')

			out47 = {l[0]: l for l in out47}
			out48 = {l[0]: l for l in out48}

			out = [out47['Sample'] + out48['Sample'][4:]]
			for s in samples:
				out.append(out47[s] + out48[s][4:])

			if save_to_file:
				if not os.path.exists(dir):
					os.makedirs(dir)
				if filename is None:
					filename = f'D47D48_samples.csv'
				with open(f'{dir}/{filename}', 'w') as fid:
					fid.write(make_csv(out))
			if print_out:
				print('\n'+pretty_table(out))
			if output == 'raw':
				return out
			elif output == 'pretty':
				return pretty_table(out)


def table_of_sessions(
	data47 = None,
	data48 = None,
	dir = 'output',
	filename = None,
	save_to_file = True,
	print_out = True,
	output = None,
	target = 'latest',
	):
	'''
	Print out, save to disk and/or return a combined table of sessions
	for a pair of `D47data` and `D48data` objects.
	***Only applicable if the sessions in `data47` and those in `data48`
	consist of the exact same sets of analyses.***

	**Parameters**

	+ `data47`: `D47data` instance
	+ `data48`: `D48data` instance
	+ `dir`: the directory in which to save the table
	+ `filename`: the name to the csv file to write to
	+ `save_to_file`: whether to save the table to disk
	+ `print_out`: whether to print out the table
	+ `output`: if set to `'pretty'`: return a pretty text table (see `pretty_table()`);
		if set to `'raw'`: return a list of list of strings
		(e.g., `[['header1', 'header2'], ['0.1', '0.2']]`)
	'''

	if data47 is None:
		if data48 is None:
			raise TypeError("Arguments must include at least one D47data() or D48data() instance.")
		else:
			return data48.table_of_sessions(
				dir = dir,
				filename = filename,
				save_to_file = save_to_file,
				print_out = print_out,
				output = output,
				target = target,
				)
	else:
		if data48 is None:
			return data47.table_of_sessions(
				dir = dir,
				filename = filename,
				save_to_file = save_to_file,
				print_out = print_out,
				output = output,
				target = target,
				)
		else:
			out47 = data47.table_of_sessions(save_to_file = False, print_out = False, output = 'raw', target = target)
			out48 = data48.table_of_sessions(save_to_file = False, print_out = False, output = 'raw', target = target)
			for k,x in enumerate(out47[0]):
				if k>7:
					out47[0][k] = out47[0][k].replace('a', 'a_47').replace('b', 'b_47').replace('c', 'c_47')
					out48[0][k] = out48[0][k].replace('a', 'a_48').replace('b', 'b_48').replace('c', 'c_48')
			out = transpose_table(transpose_table(out47) + transpose_table(out48)[7:])

			if save_to_file:
				if not os.path.exists(dir):
					os.makedirs(dir)
				if filename is None:
					filename = f'D47D48_sessions.csv'
				with open(f'{dir}/{filename}', 'w') as fid:
					fid.write(make_csv(out))
			if print_out:
				print('\n'+pretty_table(out))
			if output == 'raw':
				return out
			elif output == 'pretty':
				return pretty_table(out)


def table_of_analyses(
	data47 = None,
	data48 = None,
	dir = 'output',
	filename = None,
	save_to_file = True,
	print_out = True,
	output = None,
	):
	'''
	Print out, save to disk and/or return a combined table of analyses
	for a pair of `D47data` and `D48data` objects.

	If the sessions in `data47` and those in `data48` do not consist of
	the exact same sets of analyses, the table will have two columns
	`Session_47` and `Session_48` instead of a single `Session` column.

	**Parameters**

	+ `data47`: `D47data` instance
	+ `data48`: `D48data` instance
	+ `dir`: the directory in which to save the table
	+ `filename`: the name to the csv file to write to
	+ `save_to_file`: whether to save the table to disk
	+ `print_out`: whether to print out the table
	+ `output`: if set to `'pretty'`: return a pretty text table (see `pretty_table()`);
		if set to `'raw'`: return a list of list of strings
		(e.g., `[['header1', 'header2'], ['0.1', '0.2']]`)
	'''
	if data47 is None:
		if data48 is None:
			raise TypeError("Arguments must include at least one D47data() or D48data() instance.")
		else:
			return data48.table_of_analyses(
				dir = dir,
				filename = filename,
				save_to_file = save_to_file,
				print_out = print_out,
				output = output
				)
	else:
		if data48 is None:
			return data47.table_of_analyses(
				dir = dir,
				filename = filename,
				save_to_file = save_to_file,
				print_out = print_out,
				output = output
				)
		else:
			out47 = data47.table_of_analyses(save_to_file = False, print_out = False, output = 'raw')
			out48 = data48.table_of_analyses(save_to_file = False, print_out = False, output = 'raw')

			if [l[1] for l in out47[1:]] == [l[1] for l in out48[1:]]: # if sessions are identical
				out = transpose_table(transpose_table(out47) + transpose_table(out48)[-1:])
			else:
				out47[0][1] = 'Session_47'
				out48[0][1] = 'Session_48'
				out47 = transpose_table(out47)
				out48 = transpose_table(out48)
				out = transpose_table(out47[:2] + out48[1:2] + out47[2:] + out48[-1:])

			if save_to_file:
				if not os.path.exists(dir):
					os.makedirs(dir)
				if filename is None:
					filename = f'D47D48_analyses.csv'
				with open(f'{dir}/{filename}', 'w') as fid:
					fid.write(make_csv(out))
			if print_out:
				print('\n'+pretty_table(out))
			if output == 'raw':
				return out
			elif output == 'pretty':
				return pretty_table(out)


def _fullcovar(minresult, epsilon = 0.01, named = False):
	'''
	Construct full covariance matrix in the case of constrained parameters
	'''

	import asteval

	def f(values):
		interp = asteval.Interpreter()
		for n,v in zip(minresult.var_names, values):
			interp(f'{n} = {v}')
		for q in minresult.params:
			if minresult.params[q].expr:
				interp(f'{q} = {minresult.params[q].expr}')
		return np.array([interp.symtable[q] for q in minresult.params])

	# construct Jacobian
	J = np.zeros((minresult.nvarys, len(minresult.params)))
	X = np.array([minresult.params[p].value for p in minresult.var_names])
	sX = np.array([minresult.params[p].stderr for p in minresult.var_names])

	for j in range(minresult.nvarys):
		x1 = [_ for _ in X]
		x1[j] += epsilon * sX[j]
		x2 = [_ for _ in X]
		x2[j] -= epsilon * sX[j]
		J[j,:] = (f(x1) - f(x2)) / (2 * epsilon * sX[j])

	_names = [q for q in minresult.params]
	_covar = J.T @ minresult.covar @ J
	_se = np.diag(_covar)**.5
	_correl = _covar.copy()
	for k,s in enumerate(_se):
		if s:
			_correl[k,:] /= s
			_correl[:,k] /= s

	if named:
		_covar = {i: {j:_covar[i,j] for j in minresult.params} for i in minresult.params}
		_se = {i: _se[i] for i in minresult.params}
		_correl = {i: {j:_correl[i,j] for j in minresult.params} for i in minresult.params}

	return _names, _covar, _se, _correl


class D4xdata(list):
	'''
	Store and process data for a large set of Δ47 and/or Δ48
	analyses, usually comprising more than one analytical session.
	'''

	### 17O CORRECTION PARAMETERS
	R13_VPDB = 0.01118  # (Chang & Li, 1990)
	'''
	Absolute (13C/12C) ratio of VPDB.
	By default equal to 0.01118 ([Chang & Li, 1990](http://www.cnki.com.cn/Article/CJFDTotal-JXTW199004006.htm))
	'''

	R18_VSMOW = 0.0020052  # (Baertschi, 1976)
	'''
	Absolute (18O/16C) ratio of VSMOW.
	By default equal to 0.0020052 ([Baertschi, 1976](https://doi.org/10.1016/0012-821X(76)90115-1))
	'''

	LAMBDA_17 = 0.528  # (Barkan & Luz, 2005)
	'''
	Mass-dependent exponent for triple oxygen isotopes.
	By default equal to 0.528 ([Barkan & Luz, 2005](https://doi.org/10.1002/rcm.2250))
	'''

	R17_VSMOW = 0.00038475  # (Assonov & Brenninkmeijer, 2003, rescaled to R13_VPDB)
	'''
	Absolute (17O/16C) ratio of VSMOW.
	By default equal to 0.00038475
	([Assonov & Brenninkmeijer, 2003](https://dx.doi.org/10.1002/rcm.1011),
	rescaled to `R13_VPDB`)
	'''

	R18_VPDB = R18_VSMOW * 1.03092
	'''
	Absolute (18O/16C) ratio of VPDB.
	By definition equal to `R18_VSMOW * 1.03092`.
	'''

	R17_VPDB = R17_VSMOW * 1.03092 ** LAMBDA_17
	'''
	Absolute (17O/16C) ratio of VPDB.
	By definition equal to `R17_VSMOW * 1.03092 ** LAMBDA_17`.
	'''

	# LEVENE_REF_SAMPLE = 'ETH-3'
	# '''
	# After the Δ4x standardization step, each sample is tested to
	# assess whether the Δ4x variance within all analyses for that
	# sample differs significantly from that observed for a given reference
	# sample (using [Levene's test](https://en.wikipedia.org/wiki/Levene%27s_test),
	# which yields a p-value corresponding to the null hypothesis that the
	# underlying variances are equal).

	# `LEVENE_REF_SAMPLE` (by default equal to `'ETH-3'`) specifies which
	# sample should be used as a reference for this test.
	# '''

	ALPHA_18O_ACID_REACTION = round(np.exp(3.59 / (90 + 273.15) - 1.79e-3), 6)  # (Kim et al., 2007, calcite)
	'''
	Specifies the 18O/16O fractionation factor generally applicable
	to acid reactions in the dataset. Currently used by `D4xdata.wg()`,
	`D4xdata.standardize_d13C`, and `D4xdata.standardize_d18O`.

	By default equal to 1.008129 (calcite reacted at 90 °C,
	[Kim et al., 2007](https://dx.doi.org/10.1016/j.chemgeo.2007.08.005)).
	'''

	Nominal_d13C_VPDB = {
		'ETH-1': 2.02,
		'ETH-2': -10.17,
		'ETH-3': 1.71,
		}	# (Bernasconi et al., 2018)
	'''
	Nominal δ13C_VPDB values assigned to carbonate standards, used by
	`D4xdata.standardize_d13C()`.

	By default equal to `{'ETH-1': 2.02, 'ETH-2': -10.17, 'ETH-3': 1.71}` after
	[Bernasconi et al. (2018)](https://doi.org/10.1029/2017GC007385).
	'''

	Nominal_d18O_VPDB = {
		'ETH-1': -2.19,
		'ETH-2': -18.69,
		'ETH-3': -1.78,
		}	# (Bernasconi et al., 2018)
	'''
	Nominal δ18O_VPDB values assigned to carbonate standards, used by
	`D4xdata.standardize_d18O()`.

	By default equal to `{'ETH-1': -2.19, 'ETH-2': -18.69, 'ETH-3': -1.78}` after
	[Bernasconi et al. (2018)](https://doi.org/10.1029/2017GC007385).
	'''

	d13C_STANDARDIZATION_METHOD = '2pt'
	'''
	Method by which to standardize δ13C values:

	+ `none`: do not apply any δ13C standardization.
	+ `'1pt'`: within each session, offset all initial δ13C values so as to
	minimize the difference between final δ13C_VPDB values and
	`Nominal_d13C_VPDB` (averaged over all analyses for which `Nominal_d13C_VPDB` is defined).
	+ `'2pt'`: within each session, apply a affine trasformation to all δ13C
	values so as to minimize the difference between final δ13C_VPDB
	values and `Nominal_d13C_VPDB` (averaged over all analyses for which `Nominal_d13C_VPDB`
	is defined).
	'''

	d18O_STANDARDIZATION_METHOD = '2pt'
	'''
	Method by which to standardize δ18O values:

	+ `none`: do not apply any δ18O standardization.
	+ `'1pt'`: within each session, offset all initial δ18O values so as to
	minimize the difference between final δ18O_VPDB values and
	`Nominal_d18O_VPDB` (averaged over all analyses for which `Nominal_d18O_VPDB` is defined).
	+ `'2pt'`: within each session, apply a affine trasformation to all δ18O
	values so as to minimize the difference between final δ18O_VPDB
	values and `Nominal_d18O_VPDB` (averaged over all analyses for which `Nominal_d18O_VPDB`
	is defined).
	'''

	def __init__(self, l = [], mass = '47', logfile = '', session = 'mySession', verbose = False):
		'''
		**Parameters**

		+ `l`: a list of dictionaries, with each dictionary including at least the keys
		`Sample`, `d45`, `d46`, and `d47` or `d48`.
		+ `mass`: `'47'` or `'48'`
		+ `logfile`: if specified, write detailed logs to this file path when calling `D4xdata` methods.
		+ `session`: define session name for analyses without a `Session` key
		+ `verbose`: if `True`, print out detailed logs when calling `D4xdata` methods.

		Returns a `D4xdata` object derived from `list`.
		'''
		self._4x = mass
		self.verbose = verbose
		self.prefix = 'D4xdata'
		self.logfile = logfile
		list.__init__(self, l)
		self.repeatability = {}
		self.standardization = {}
		self.refresh(session = session)

	@property
	def fixed_RMs(self):
		return {
			k: v for k,v in self.RMs.items()
			if not isinstance(v, tuple)
		}

	@property
	def loose_RMs(self):
		return {
			k: v for k,v in self.RMs.items()
			if isinstance(v, tuple)
		}

	@property
	def bare_RMs(self):
		return {
			k: v[0] if isinstance(v, tuple) else v
			for k,v in self.RMs.items()
		}

	def make_verbal(oldfun):
		'''
		Decorator: allow temporarily changing `self.prefix` and overriding `self.verbose`.
		'''
		@wraps(oldfun)
		def newfun(*args, verbose = '', **kwargs):
			myself = args[0]
			oldprefix = myself.prefix
			myself.prefix = oldfun.__name__
			if verbose != '':
				oldverbose = myself.verbose
				myself.verbose = verbose
			out = oldfun(*args, **kwargs)
			myself.prefix = oldprefix
			if verbose != '':
				myself.verbose = oldverbose
			return out
		return newfun


	def msg(self, txt):
		'''
		Log a message to `self.logfile`, and print it out if `verbose = True`
		'''
		self.log(txt)
		if self.verbose:
			print(f'{f"[{self.prefix}]":<16} {txt}')


	def vmsg(self, txt):
		'''
		Log a message to `self.logfile` and print it out
		'''
		self.log(txt)
		print(txt)


	def log(self, *txts):
		'''
		Log a message to `self.logfile`
		'''
		if self.logfile:
			with open(self.logfile, 'a') as fid:
				for txt in txts:
					fid.write(f'\n{dt.now().strftime("%Y-%m-%d %H:%M:%S")} {f"[{self.prefix}]":<16} {txt}')


	def refresh(self, session = 'mySession'):
		'''
		Update `self.sessions`, `self.samples`, `self.anchors`, and `self.unknowns`.
		'''
		self.fill_in_missing_info(session = session)
		self.refresh_sessions()
		self.refresh_samples()


	def refresh_sessions(self):
		'''
		Update `self.sessions` and set `scrambling_drift`, `slope_drift`, and `wg_drift`
		to `False` for all sessions.
		'''
		self.sessions = {
			s: {'data': [r for r in self if r['Session'] == s]}
			for s in sorted({r['Session'] for r in self})
			}
		for s in self.sessions:
			self.sessions[s]['N'] = len(self.sessions[s]['data'])
			self.sessions[s]['scrambling_drift'] = False
			self.sessions[s]['slope_drift'] = False
			self.sessions[s]['wg_drift'] = False
			self.sessions[s]['d13C_standardization_method'] = self.d13C_STANDARDIZATION_METHOD
			self.sessions[s]['d18O_standardization_method'] = self.d18O_STANDARDIZATION_METHOD


	def refresh_samples(self):
		'''
		Define `self.samples`, `self.anchors`, and `self.unknowns`.
		'''
		self.samples = {
			s: {'data': [r for r in self if r['Sample'] == s]}
			for s in sorted({r['Sample'] for r in self})
			}
		self.anchors = {s: self.samples[s] for s in self.samples if s in self.RMs}
		self.fixed_anchors = {s: self.samples[s] for s in self.samples if s in self.fixed_RMs}
		self.loose_anchors = {s: self.samples[s] for s in self.samples if s in self.loose_RMs}
		self.unknowns = {s: self.samples[s] for s in self.samples if s not in self.RMs}

	def read(self, filename, sep = '', session = ''):
		'''
		Read file in csv format to load data into a `D47data` object.

		In the csv file, spaces before and after field separators (`','` by default)
		are optional. Each line corresponds to a single analysis.

		The required fields are:

		+ `UID`: a unique identifier
		+ `Session`: an identifier for the analytical session
		+ `Sample`: a sample identifier
		+ `d45`, `d46`, and at least one of `d47` or `d48`: the working-gas delta values

		Independently known oxygen-17 anomalies may be provided as `D17O` (in ‰ relative to
		VSMOW, λ = `self.LAMBDA_17`), and are otherwise assumed to be zero. Working-gas deltas `d47`, `d48`
		and `d49` are optional, and set to NaN by default.

		**Parameters**

		+ `fileneme`: the path of the file to read
		+ `sep`: csv separator delimiting the fields
		+ `session`: set `Session` field to this string for all analyses
		'''
		with open(filename) as fid:
			self.input(fid.read(), sep = sep, session = session)


	def input(self, txt, sep = '', session = ''):
		'''
		Read `txt` string in csv format to load analysis data into a `D47data` object.

		In the csv string, spaces before and after field separators (`','` by default)
		are optional. Each line corresponds to a single analysis.

		The required fields are:

		+ `UID`: a unique identifier
		+ `Session`: an identifier for the analytical session
		+ `Sample`: a sample identifier
		+ `d45`, `d46`, and at least one of `d47` or `d48`: the working-gas delta values

		Independently known oxygen-17 anomalies may be provided as `D17O` (in ‰ relative to
		VSMOW, λ = `self.LAMBDA_17`), and are otherwise assumed to be zero. Working-gas deltas `d47`, `d48`
		and `d49` are optional, and set to NaN by default.

		**Parameters**

		+ `txt`: the csv string to read
		+ `sep`: csv separator delimiting the fields. By default, use `,`, `;`, or `\t`,
		whichever appers most often in `txt`.
		+ `session`: set `Session` field to this string for all analyses
		'''
		if sep == '':
			sep = sorted(',;\t', key = lambda x: - txt.count(x))[0]
		txt = [[x.strip() for x in l.split(sep)] for l in txt.splitlines() if l.strip()]
		data = [{k: v if k in ['UID', 'Session', 'Sample'] else smart_type(v) for k,v in zip(txt[0], l) if v != ''} for l in txt[1:]]

		if session != '':
			for r in data:
				r['Session'] = session

		self += data
		self.refresh()


	@make_verbal
	def wg(self,
		samples = None,
		session_groups = None,
	):
		'''
		Compute bulk composition of the working gas for each session based (by default)
		on the carbonate standards defined in both `self.Nominal_d13C_VPDB` and
		`self.Nominal_d18O_VPDB`.

		**Parameters**

		+ `samples`: A list of samples specifying the subset of samples (defined in both
		`self.Nominal_d13C_VPDB` and `self.Nominal_d18O_VPDB`) which will be considered
		when computing the working gas. By default, use all samples defined both in
		`self.Nominal_d13C_VPDB` and `self.Nominal_d18O_VPDB`.
		+ `session_groups`: a list of lists of sessions
		(e.g., `[['session1', 'session2'], ['session3', 'session4', 'session5']]`)
		specifying which sessions groups, if any, have the exact same WG composition.
		If set to `'all'`, force all sessions to have the same WG composition (use with
		caution and on short time scales, since the WG may drift slowly a long time scales).
		'''

		self.msg('Computing WG composition:')

		a18_acid = self.ALPHA_18O_ACID_REACTION

		if samples is None:
			samples = [s for s in self.Nominal_d13C_VPDB if s in self.Nominal_d18O_VPDB]
		if session_groups is None:
			session_groups = [[s] for s in self.sessions]
		elif session_groups == 'all':
			session_groups = [[s for s in self.sessions]]

		samples = [s for s in samples if s in self.Nominal_d13C_VPDB and s in self.Nominal_d18O_VPDB]
		R45R46_standards = {}
		for sample in samples:
			d13C_vpdb = self.Nominal_d13C_VPDB[sample]
			d18O_vpdb = self.Nominal_d18O_VPDB[sample]
			R13_s = self.R13_VPDB * (1 + d13C_vpdb / 1000)
			R17_s = self.R17_VPDB * ((1 + d18O_vpdb / 1000) * a18_acid) ** self.LAMBDA_17
			R18_s = self.R18_VPDB * (1 + d18O_vpdb / 1000) * a18_acid

			C12_s = 1 / (1 + R13_s)
			C13_s = R13_s / (1 + R13_s)
			C16_s = 1 / (1 + R17_s + R18_s)
			C17_s = R17_s / (1 + R17_s + R18_s)
			C18_s = R18_s / (1 + R17_s + R18_s)

			C626_s = C12_s * C16_s ** 2
			C627_s = 2 * C12_s * C16_s * C17_s
			C628_s = 2 * C12_s * C16_s * C18_s
			C636_s = C13_s * C16_s ** 2
			C637_s = 2 * C13_s * C16_s * C17_s
			C727_s = C12_s * C17_s ** 2

			R45_s = (C627_s + C636_s) / C626_s
			R46_s = (C628_s + C637_s + C727_s) / C626_s
			R45R46_standards[sample] = (R45_s, R46_s)

		for sg in session_groups:
			db = [r for s in sg for r in self.sessions[s]['data'] if r['Sample'] in samples]
			assert db, f'No sample from {samples} found in session group {sg}.'

			X = [r['d45'] for r in db]
			Y = [R45R46_standards[r['Sample']][0] for r in db]
			x1, x2 = np.min(X), np.max(X)

			if x1 < x2:
				wgcoord = x1/(x1-x2)
			else:
				wgcoord = 999

			if wgcoord < -.5 or wgcoord > 1.5:
				# unreasonable to extrapolate to d45 = 0
				R45_wg = np.mean([y/(1+x/1000) for x,y in zip(X,Y)])
			else :
				# d45 = 0 is reasonably well bracketed
				R45_wg = np.polyfit(X, Y, 1)[1]

			X = [r['d46'] for r in db]
			Y = [R45R46_standards[r['Sample']][1] for r in db]
			x1, x2 = np.min(X), np.max(X)

			if x1 < x2:
				wgcoord = x1/(x1-x2)
			else:
				wgcoord = 999

			if wgcoord < -.5 or wgcoord > 1.5:
				# unreasonable to extrapolate to d46 = 0
				R46_wg = np.mean([y/(1+x/1000) for x,y in zip(X,Y)])
			else :
				# d46 = 0 is reasonably well bracketed
				R46_wg = np.polyfit(X, Y, 1)[1]

			d13Cwg_VPDB, d18Owg_VSMOW = self.compute_bulk_delta(R45_wg, R46_wg)

			for s in sg:
				self.msg(f'Sessions {s} WG:   δ13C_VPDB = {d13Cwg_VPDB:.3f}   δ18O_VSMOW = {d18Owg_VSMOW:.3f}')

				self.sessions[s]['d13Cwg_VPDB'] = d13Cwg_VPDB
				self.sessions[s]['d18Owg_VSMOW'] = d18Owg_VSMOW
				for r in self.sessions[s]['data']:
					r['d13Cwg_VPDB'] = d13Cwg_VPDB
					r['d18Owg_VSMOW'] = d18Owg_VSMOW


	def compute_bulk_delta(self, R45, R46, D17O = 0):
		'''
		Compute δ13C_VPDB and δ18O_VSMOW,
		by solving the generalized form of equation (17) from
		[Brand et al. (2010)](https://doi.org/10.1351/PAC-REP-09-01-05),
		assuming that δ18O_VSMOW is not too big (0 ± 50 ‰) and
		solving the corresponding second-order Taylor polynomial.
		(Appendix A of [Daëron et al., 2016](https://doi.org/10.1016/j.chemgeo.2016.08.014))
		'''

		K = np.exp(D17O / 1000) * self.R17_VSMOW * self.R18_VSMOW ** -self.LAMBDA_17

		A = -3 * K ** 2 * self.R18_VSMOW ** (2 * self.LAMBDA_17)
		B = 2 * K * R45 * self.R18_VSMOW ** self.LAMBDA_17
		C = 2 * self.R18_VSMOW
		D = -R46

		aa = A * self.LAMBDA_17 * (2 * self.LAMBDA_17 - 1) + B * self.LAMBDA_17 * (self.LAMBDA_17 - 1) / 2
		bb = 2 * A * self.LAMBDA_17 + B * self.LAMBDA_17 + C
		cc = A + B + C + D

		d18O_VSMOW = 1000 * (-bb + (bb ** 2 - 4 * aa * cc) ** .5) / (2 * aa)

		R18 = (1 + d18O_VSMOW / 1000) * self.R18_VSMOW
		R17 = K * R18 ** self.LAMBDA_17
		R13 = R45 - 2 * R17

		d13C_VPDB = 1000 * (R13 / self.R13_VPDB - 1)

		return d13C_VPDB, d18O_VSMOW


	@make_verbal
	def crunch(self, verbose = ''):
		'''
		Compute bulk composition and raw clumped isotope anomalies for all analyses.
		'''
		for r in self:
			self.compute_bulk_and_clumping_deltas(r)
		self.standardize_d13C()
		self.standardize_d18O()
		self.msg(f"Crunched {len(self)} analyses.")


	def fill_in_missing_info(self, session = 'mySession'):
		'''
		Fill in optional fields with default values
		'''
		for i,r in enumerate(self):
			if 'D17O' not in r:
				r['D17O'] = 0.
			if 'UID' not in r:
				r['UID'] = f'{i+1}'
			if 'Session' not in r:
				r['Session'] = session
			for k in ['d47', 'd48', 'd49']:
				if k not in r:
					r[k] = np.nan


	def standardize_d13C(self):
		'''
		Perform δ13C standadization within each session `s` according to
		`self.sessions[s]['d13C_standardization_method']`, which is defined by default
		by `D47data.refresh_sessions()`as equal to `self.d13C_STANDARDIZATION_METHOD`, but
		may be redefined abitrarily at a later stage.
		'''
		for s in self.sessions:
			if self.sessions[s]['d13C_standardization_method'] in ['1pt', '2pt']:
				XY = [(r['d13C_VPDB'], self.Nominal_d13C_VPDB[r['Sample']]) for r in self.sessions[s]['data'] if r['Sample'] in self.Nominal_d13C_VPDB]
				X,Y = zip(*XY)
				if self.sessions[s]['d13C_standardization_method'] == '1pt':
					offset = np.mean(Y) - np.mean(X)
					for r in self.sessions[s]['data']:
						r['d13C_VPDB'] += offset
				elif self.sessions[s]['d13C_standardization_method'] == '2pt':
					a,b = np.polyfit(X,Y,1)
					for r in self.sessions[s]['data']:
						r['d13C_VPDB'] = a * r['d13C_VPDB'] + b

	def standardize_d18O(self):
		'''
		Perform δ18O standadization within each session `s` according to
		`self.ALPHA_18O_ACID_REACTION` and `self.sessions[s]['d18O_standardization_method']`,
		which is defined by default by `D47data.refresh_sessions()`as equal to
		`self.d18O_STANDARDIZATION_METHOD`, but may be redefined abitrarily at a later stage.
		'''
		for s in self.sessions:
			if self.sessions[s]['d18O_standardization_method'] in ['1pt', '2pt']:
				XY = [(r['d18O_VSMOW'], self.Nominal_d18O_VPDB[r['Sample']]) for r in self.sessions[s]['data'] if r['Sample'] in self.Nominal_d18O_VPDB]
				X,Y = zip(*XY)
				Y = [(1000+y) * self.R18_VPDB * self.ALPHA_18O_ACID_REACTION / self.R18_VSMOW - 1000 for y in Y]
				if self.sessions[s]['d18O_standardization_method'] == '1pt':
					offset = np.mean(Y) - np.mean(X)
					for r in self.sessions[s]['data']:
						r['d18O_VSMOW'] += offset
				elif self.sessions[s]['d18O_standardization_method'] == '2pt':
					a,b = np.polyfit(X,Y,1)
					for r in self.sessions[s]['data']:
						r['d18O_VSMOW'] = a * r['d18O_VSMOW'] + b


	def compute_bulk_and_clumping_deltas(self, r):
		'''
		Compute δ13C_VPDB, δ18O_VSMOW, and raw Δ47, Δ48, Δ49 values for a single analysis `r`.
		'''

		# Compute working gas R13, R18, and isobar ratios
		R13_wg = self.R13_VPDB * (1 + r['d13Cwg_VPDB'] / 1000)
		R18_wg = self.R18_VSMOW * (1 + r['d18Owg_VSMOW'] / 1000)
		R45_wg, R46_wg, R47_wg, R48_wg, R49_wg = self.compute_isobar_ratios(R13_wg, R18_wg)

		# Compute analyte isobar ratios
		R45 = (1 + r['d45'] / 1000) * R45_wg
		R46 = (1 + r['d46'] / 1000) * R46_wg
		R47 = (1 + r['d47'] / 1000) * R47_wg
		R48 = (1 + r['d48'] / 1000) * R48_wg
		R49 = (1 + r['d49'] / 1000) * R49_wg

		r['d13C_VPDB'], r['d18O_VSMOW'] = self.compute_bulk_delta(R45, R46, D17O = r['D17O'])
		R13 = (1 + r['d13C_VPDB'] / 1000) * self.R13_VPDB
		R18 = (1 + r['d18O_VSMOW'] / 1000) * self.R18_VSMOW

		# Compute stochastic isobar ratios of the analyte
		R45stoch, R46stoch, R47stoch, R48stoch, R49stoch = self.compute_isobar_ratios(
			R13, R18, D17O = r['D17O']
		)

		# Check that R45/R45stoch and R46/R46stoch are undistinguishable from 1,
		# and raise a warning if the corresponding anomalies exceed 0.02 ppm.
		if (R45 / R45stoch - 1) > 5e-8:
			self.vmsg(f'This is unexpected: R45/R45stoch - 1 = {1e6 * (R45 / R45stoch - 1):.3f} ppm')
		if (R46 / R46stoch - 1) > 5e-8:
			self.vmsg(f'This is unexpected: R46/R46stoch - 1 = {1e6 * (R46 / R46stoch - 1):.3f} ppm')

		# Compute raw clumped isotope anomalies
		r['D47raw'] = 1000 * (R47 / R47stoch - 1)
		r['D48raw'] = 1000 * (R48 / R48stoch - 1)
		r['D49raw'] = 1000 * (R49 / R49stoch - 1)


	def compute_isobar_ratios(self, R13, R18, D17O=0, D47=0, D48=0, D49=0):
		'''
		Compute isobar ratios for a sample with isotopic ratios `R13` and `R18`,
		optionally accounting for non-zero values of Δ17O (`D17O`) and clumped isotope
		anomalies (`D47`, `D48`, `D49`), all expressed in permil.
		'''

		# Compute R17
		R17 = self.R17_VSMOW * np.exp(D17O / 1000) * (R18 / self.R18_VSMOW) ** self.LAMBDA_17

		# Compute isotope concentrations
		C12 = (1 + R13) ** -1
		C13 = C12 * R13
		C16 = (1 + R17 + R18) ** -1
		C17 = C16 * R17
		C18 = C16 * R18

		# Compute stochastic isotopologue concentrations
		C626 = C16 * C12 * C16
		C627 = C16 * C12 * C17 * 2
		C628 = C16 * C12 * C18 * 2
		C636 = C16 * C13 * C16
		C637 = C16 * C13 * C17 * 2
		C638 = C16 * C13 * C18 * 2
		C727 = C17 * C12 * C17
		C728 = C17 * C12 * C18 * 2
		C737 = C17 * C13 * C17
		C738 = C17 * C13 * C18 * 2
		C828 = C18 * C12 * C18
		C838 = C18 * C13 * C18

		# Compute stochastic isobar ratios
		R45 = (C636 + C627) / C626
		R46 = (C628 + C637 + C727) / C626
		R47 = (C638 + C728 + C737) / C626
		R48 = (C738 + C828) / C626
		R49 = C838 / C626

		# Account for stochastic anomalies
		R47 *= 1 + D47 / 1000
		R48 *= 1 + D48 / 1000
		R49 *= 1 + D49 / 1000

		# Return isobar ratios
		return R45, R46, R47, R48, R49


	def split_samples(self, samples_to_split = 'all', grouping = 'by_session'):
		'''
		Split unknown samples by UID (treat all analyses as different samples)
		or by session (treat analyses of a given sample in different sessions as
		different samples).

		**Parameters**

		+ `samples_to_split`: a list of samples to split, e.g., `['IAEA-C1', 'IAEA-C2']`
		+ `grouping`: `by_uid` | `by_session`
		'''
		if samples_to_split == 'all':
			samples_to_split = [s for s in self.unknowns]
		gkeys = {'by_uid':'UID', 'by_session':'Session'}
		self.grouping = grouping.lower()
		if self.grouping in gkeys:
			gkey = gkeys[self.grouping]
		for r in self:
			if r['Sample'] in samples_to_split:
				r['Sample_original'] = r['Sample']
				r['Sample'] = f"{r['Sample']}__{r[gkey]}"
			elif r['Sample'] in self.unknowns:
				r['Sample_original'] = r['Sample']
		self.refresh_samples()


	def unsplit_samples(self, tables = False):
		'''
		Reverse the effects of `D47data.split_samples()`.

		This should only be used after `D4xdata.standardize()` with `method='pooled'`.

		After `D4xdata.standardize()` with `method='indep_sessions'`, one should
		probably use `D4xdata.combine_samples()` instead to reverse the effects of
		`D47data.split_samples()` with `grouping='by_uid'`, or `w_avg()` to reverse the
		effects of `D47data.split_samples()` with `grouping='by_sessions'` (because in
		that case session-averaged Δ4x values are statistically independent).
		'''
		unknowns_old = sorted({s for s in self.unknowns})
		CM_old = self.standardization.covar[:,:]
		VD_old = self.standardization.params.valuesdict().copy()
		vars_old = self.standardization.var_names

		unknowns_new = sorted({r['Sample_original'] for r in self if 'Sample_original' in r})

		Ns = len(vars_old) - len(unknowns_old)
		vars_new = vars_old[:Ns] + [f'D{self._4x}_{pf(u)}' for u in unknowns_new]
		VD_new = {k: VD_old[k] for k in vars_old[:Ns]}

		W = np.zeros((len(vars_new), len(vars_old)))
		W[:Ns,:Ns] = np.eye(Ns)
		for u in unknowns_new:
			splits = sorted({r['Sample'] for r in self if 'Sample_original' in r and r['Sample_original'] == u})
			if self.grouping == 'by_session':
				weights = [self.samples[s][f'SE_D{self._4x}']**-2 for s in splits]
			elif self.grouping == 'by_uid':
				weights = [1 for s in splits]
			sw = sum(weights)
			weights = [w/sw for w in weights]
			W[vars_new.index(f'D{self._4x}_{pf(u)}'),[vars_old.index(f'D{self._4x}_{pf(s)}') for s in splits]] = weights[:]

		CM_new = W @ CM_old @ W.T
		V = W @ np.array([[VD_old[k]] for k in vars_old])
		VD_new = {k:v[0] for k,v in zip(vars_new, V)}

		self.standardization.covar = CM_new
		self.standardization.params.valuesdict = lambda : VD_new
		self.standardization.var_names = vars_new

		for r in self:
			if r['Sample'] in self.unknowns:
				r['Sample_split'] = r['Sample']
				r['Sample'] = r['Sample_original']

		self.refresh_samples()
		self.consolidate_samples()
		self.repeatabilities()

		if tables:
			self.table_of_analyses()
			self.table_of_samples()

	def assign_timestamps(self):
		'''
		Assign a time field `t` of type `float` to each analysis.

		If `TimeTag` is one of the data fields, `t` is equal within a given session
		to `TimeTag` minus the mean value of `TimeTag` for that session.
		Otherwise, `TimeTag` is by default equal to the index of each analysis
		in the dataset and `t` is defined as above.
		'''
		for session in self.sessions:
			sdata = self.sessions[session]['data']
			try:
				t0 = np.mean([r['TimeTag'] for r in sdata])
				for r in sdata:
					r['t'] = r['TimeTag'] - t0
			except KeyError:
				t0 = (len(sdata)-1)/2
				for t,r in enumerate(sdata):
					r['t'] = t - t0


	def report(self):
		'''
		Prints a report on the standardization fit.
		Only applicable after `D4xdata.standardize(method='pooled')`.
		'''
		report_fit(self.standardization)


	def combine_samples(self, sample_groups):
		'''
		Combine analyses of different samples to compute weighted average Δ4x
		and new error (co)variances corresponding to the groups defined by the `sample_groups`
		dictionary.

		Caution: samples are weighted by number of replicate analyses, which is a
		reasonable default behavior but is not always optimal (e.g., in the case of strongly
		correlated analytical errors for one or more samples).

		Returns a tuplet of:

		+ the list of group names
		+ an array of the corresponding Δ4x values
		+ the corresponding (co)variance matrix

		**Parameters**

		+ `sample_groups`: a dictionary of the form:
		```py
		{'group1': ['sample_1', 'sample_2'],
		 'group2': ['sample_3', 'sample_4', 'sample_5']}
		```
		'''

		samples = [s for k in sorted(sample_groups.keys()) for s in sorted(sample_groups[k])]
		groups = sorted(sample_groups.keys())
		group_total_weights = {k: sum([self.samples[s]['N'] for s in sample_groups[k]]) for k in groups}
		D4x_old = np.array([[self.samples[x][f'D{self._4x}']] for x in samples])
		CM_old = np.array([[self.sample_D4x_covar(x,y) for x in samples] for y in samples])
		W = np.array([
			[self.samples[i]['N']/group_total_weights[j] if i in sample_groups[j] else 0 for i in samples]
			for j in groups])
		D4x_new = W @ D4x_old
		CM_new = W @ CM_old @ W.T

		return groups, D4x_new[:,0], CM_new


	@make_verbal
	def standardize(self,
		method = 'pooled',
		weighted_sessions = [],
		consolidate = True,
		consolidate_tables = False,
		consolidate_plots = False,
		mcmc_sample_kw = {},
		constraints = {},
		):
		'''
		Compute absolute Δ4x values for all replicate analyses and for sample averages.

		**Parameters**

		+ `method`:
			- `'pooled'`: processes all sessions in a single step, assuming that all samples
				(anchors and unknowns alike) are homogeneous, i.e. that their true Δ4x values do
				not change between sessions ([Daëron, 2021](https://doi.org/10.1029/2020GC009592)).
			- `'bayes'`: use Bayesian approach to standardization, which accounts for anchors with
				uncertain nominal Δ4x values. See "Bayesian requirements" below.
			- `'indep_sessions'`: processes each session independently, based only on anchor analyses.

		> [!CAUTION]
		> `method = 'indep_sessions'` will eventually be deprecated.

		+ `weighted_sessions` (not implemented for `method = 'bayes'`):
			grouping of sessions (e.g., `[['S1', 'S2'], ['S3', 'S4', 'S5']]`) assumed
			to share the same pooled reproducibility of Δ4x measurements. This is intended for cases where
			different datasets with different intrumental performance levels (e.g., from different mass
			spectrometers) are to be combined.
		+ `consolidate`: Whether to collect information about samples, sessions and repeatabilities after standardization.
			Not to be changed unless you know exactly what you're doing.
		+ `consolidate_tables`: Whether to ouput tables during `D4xdata.consolidate()`.
		+ `consolidate_plots`: Whether to ouput session plots during `D4xdata.consolidate()`.
		+ `weak_anchors` (only applicable when `method = 'bayes'`): dict of `{sample: (mu, sigma)}` items, with `(mu, sigma)`
			being the mean and 1-σ uncertainty of the nominal Δ4x value of sample `sample`.
			Example: `weak_anchors = {'ETH-4': (0.4511, 0.0011), 'TAC-1': (0.7, 0.02)}`
			(NB: these are made-up values for `TAC-1` for now).
		+ `mcmc_sample_kw` (only applicable when `method = 'bayes'`): dict of parameters to be passed on to `pymc.sample()`,
			e.g., `{'draws': 2000, 'random_seed': 1234}`.
		+ `constraints`: specify additional mathematical constraints linking the standardization parameters
			(session parameters and/or unknown samples' Δ4x values). Currently, the format for these constraints
			depends on the `method` parameter:
			- `method = 'pooled'`: a dict of `{param: expr}` items, where `param` is the name of the parameter
				(normalized using `pf()`) and `expr` is a mathematical expression used to compute `param` exactly.
				Internally impremented by the `lmfit` package.
			- `method = 'bayes'`: also a dict of `{param: expr}` items, but using a different, bracket-based notation
				for fit parameters (and without `pf()` normalization). Internally impremented by the `sympy` package.

		Below are examples for the two cases above:

		```py
		# POOLED METHOD
		D4xdata.standardize(
		  method = 'pooled',
		  constraints = {
		    # ensure that Session_01 and Session_02 share the same WG D4x value:
		    'c_Session_02': 'c_Session_01 / a_Session_01 * a_Session_02',
		    # Force the correct scaling between 25 °C equilibrated and 1000 °C heated gases:
		    'D4x_EG_25C' : 'D4x_HG_1000C + 0.893',
		  }
		)

		# BAYES METHOD
		D4xdata.standardize(
		  method = 'bayes',
		  constraints = {
		    # ensure that Session_01 and Session_02 share the same WG D4x value:
		    "c['Session_02']": "c['Session_01'] / a['Session_01'] * a['Session_02']",
		    # Force the correct scaling between 25 °C equilibrated and 1000 °C heated gases:
		    "D47['EG_25C']" : "D47['HG_1000C'] + 0.893",
		  }
		)
		}
		```
		'''

		self.assign_timestamps()

		if method == 'pooled':
			self._pooled_standardization(
				constraints = constraints,
				weighted_sessions = weighted_sessions,
				consolidate = consolidate,
				consolidate_tables = consolidate_tables,
				consolidate_plots = consolidate_plots,
			)

		elif method == 'bayes':
			self._bayesian_standardization(
				constraints = constraints,
				mcmc_sample_kw = mcmc_sample_kw,
				consolidate = consolidate,
				consolidate_tables = consolidate_tables,
				consolidate_plots = consolidate_plots,
			)

	def _pooled_standardization(
		self,
		constraints = {},
		weighted_sessions = [],
		consolidate = True,
		consolidate_tables = False,
		consolidate_plots = False,
	):

		if weighted_sessions:
			for session_group in weighted_sessions:
				if self._4x == '47':
					X = D47data([r for r in self if r['Session'] in session_group])
				elif self._4x == '48':
					X = D48data([r for r in self if r['Session'] in session_group])
				X.RMs = self.RMs.copy()
				X.refresh()
				result = X.standardize(method = 'pooled', weighted_sessions = [], consolidate = False)
				w = np.sqrt(result.redchi)
				self.msg(f'Session group {session_group} MRSWD = {w:.4f}')
				for r in X:
					r[f'wD{self._4x}raw'] *= w
		else:
			self.msg(f'All D{self._4x}raw weights set to 1 ‰')
			for r in self:
				r[f'wD{self._4x}raw'] = 1.

		params = Parameters()
		for k,session in enumerate(self.sessions):
			self.msg(f"Session {session}: scrambling_drift is {self.sessions[session]['scrambling_drift']}.")
			self.msg(f"Session {session}: slope_drift is {self.sessions[session]['slope_drift']}.")
			self.msg(f"Session {session}: wg_drift is {self.sessions[session]['wg_drift']}.")
			s = pf(session)
			params.add(f'a_{s}', value = 0.9)
			params.add(f'b_{s}', value = 0.)
			params.add(f'c_{s}', value = -0.9)
			params.add(f'a2_{s}', value = 0.,
# 					vary = self.sessions[session]['scrambling_drift'],
				)
			params.add(f'b2_{s}', value = 0.,
# 					vary = self.sessions[session]['slope_drift'],
				)
			params.add(f'c2_{s}', value = 0.,
# 					vary = self.sessions[session]['wg_drift'],
				)
			if not self.sessions[session]['scrambling_drift']:
				params[f'a2_{s}'].expr = '0'
			if not self.sessions[session]['slope_drift']:
				params[f'b2_{s}'].expr = '0'
			if not self.sessions[session]['wg_drift']:
				params[f'c2_{s}'].expr = '0'

		for sample in self.unknowns:
			params.add(f'D{self._4x}_{pf(sample)}', value = 0.5)

		for k in constraints:
			params[k].expr = constraints[k]

		def residuals(p):
			R = []
			for r in self:
				session = pf(r['Session'])
				sample = pf(r['Sample'])
				if r['Sample'] in self.bare_RMs:
					R += [ (
						r[f'D{self._4x}raw'] - (
							p[f'a_{session}'] * self.bare_RMs[r['Sample']]
							+ p[f'b_{session}'] * r[f'd{self._4x}']
							+	p[f'c_{session}']
							+ r['t'] * (
								p[f'a2_{session}'] * self.bare_RMs[r['Sample']]
								+ p[f'b2_{session}'] * r[f'd{self._4x}']
								+	p[f'c2_{session}']
								)
							)
						) / r[f'wD{self._4x}raw'] ]
				else:
					R += [ (
						r[f'D{self._4x}raw'] - (
							p[f'a_{session}'] * p[f'D{self._4x}_{sample}']
							+ p[f'b_{session}'] * r[f'd{self._4x}']
							+	p[f'c_{session}']
							+ r['t'] * (
								p[f'a2_{session}'] * p[f'D{self._4x}_{sample}']
								+ p[f'b2_{session}'] * r[f'd{self._4x}']
								+	p[f'c2_{session}']
								)
							)
						) / r[f'wD{self._4x}raw'] ]
			return R

		M = Minimizer(residuals, params)
		result = M.least_squares()
		result.J = result.jac                                # best fit Jacobian
		result.Q, _ = np.linalg.qr(result.J)                 # QR decomposition of J
		result.h = np.einsum('ij,ij->i', result.Q, result.Q) # leverage h matrix, aligned row-for-row with result.residual

		# sanity check
		assert np.isclose(result.h.sum(), len(result.var_names))

		self.standardization['pooled'] = dict(method = 'pooled')
		S = self.standardization['pooled']
		self.standardization['latest'] = S

		S['lmfit'] = result
		S['Nf'] = result.nfree
		S['t95'] = tstudent.ppf(1 - 0.05/2, result.nfree)

		new_names, new_covar, new_se = _fullcovar(result)[:3]

		uparams = uncertainties.correlated_values(
			[result.params.valuesdict()[k] for k in new_names],
			new_covar,
		)
		S['uparams'] = {k: v for k,v in zip(new_names, uparams)}

		for r in self:
			s = pf(r["Session"])
			a = result.params.valuesdict()[f'a_{s}']
			b = result.params.valuesdict()[f'b_{s}']
			c = result.params.valuesdict()[f'c_{s}']
			a2 = result.params.valuesdict()[f'a2_{s}']
			b2 = result.params.valuesdict()[f'b2_{s}']
			c2 = result.params.valuesdict()[f'c2_{s}']
			r[f'D{self._4x}'] = (r[f'D{self._4x}raw'] - c - b * r[f'd{self._4x}'] - c2 * r['t'] - b2 * r['t'] * r[f'd{self._4x}']) / (a + a2 * r['t'])

		S['samples'] = {}
		_D4x_ = f'D{self._4x}'
		for sample in self.samples:
			S['samples'][sample] = {}
			if sample in self.bare_RMs:
				with warnings.catch_warnings():
					warnings.filterwarnings("ignore", message="Using UFloat objects with std_dev==0")
					S['samples'][sample][_D4x_] = uncertainties.ufloat(self.bare_RMs[sample], 0.)
			else:
				S['samples'][sample][_D4x_] = S['uparams'][f'{_D4x_}_{pf(sample)}']
			S['samples'][sample][f'95CL_{_D4x_}'] = S['samples'][sample][_D4x_].s * S['t95']

		S['sessions'] = {}
		for session in self.sessions:
			S['sessions'][session] = {'Np': 3}
			for k in ['scrambling', 'slope', 'wg']:
				if self.sessions[session][f'{k}_drift']:
					S['sessions'][session]['Np'] += 1

			for k in ['a', 'b', 'c', 'a2', 'b2', 'c2']:
				S['sessions'][session][k] = S['uparams'][f'{k}_{pf(session)}']

			S['sessions'][session]['CM'] = np.array(uncertainties.covariance_matrix([
				S['sessions'][session][k]
				for k in ['a', 'b', 'c', 'a2', 'b2', 'c2']
			]))


		if consolidate:
			self.consolidate(target = 'pooled', tables = consolidate_tables, plots = consolidate_plots)

		return result

	def _bayesian_standardization(
		self,
		constraints = {},
		sigma_session_groups = None,
		consolidate = True,
		consolidate_tables = False,
		consolidate_plots = False,
		mcmc_sample_kw = {},
		default_mcmc_sample_kw = {
			'draws': 2000,
			'tune': 1000,
			'random_seed': None,
			'target_accept': 0.95,
		},
	):
		'''
		Compute absolute Δ4x values as when calling `standardize()`, but using bayesian methods
		accounting for uncertainties in the nominal Δ4x values of "loose" anchors.

		**Parameters**

		+ `constraints`: a dict specifying exact algebraic relationships between elements of
		  `a`, `b`, `c`, or `D4x`, keyed by the constrained element and valued by an expression
		  string defining it, e.g.:

		    constraints = {
		        "c['Session_02']": "c['Session_01'] / a['Session_01'] * a['Session_02']",
		        "D47['Sample_02']": "D47['Sample_01'] + 0.987",
		    }

		Each key/value may reference elements of `a`, `b`, `c`, or `D{4x}` by session/sample
		label (e.g. `a['Session_01']`) or by integer position (e.g. `a[0]`).

		+ `sigma_session_groups`: a list of lists of session names specifying which sessions
		  share a common value of `sigma` (the corrected-Δ4x-space analytical noise), e.g.:

		    sigma_session_groups = [
		        ['Session_01', 'Session_02'],
		        ['Session_03', 'Session_04', 'Session_05'],
		    ]

		Sessions not listed in any group are silently collected into one additional group.
		If `None` (default), a single group containing all sessions is used, i.e. `sigma`
		is a single value shared by every session.
		'''

		# lazy imports:
		import re
		import pymc as pm
		import arviz as az
		import sympy as sp
		import pytensor.tensor as pt
		from sympy.parsing.sympy_parser import parse_expr, standard_transformations

		_d4x_ = f'd{self._4x}' # 'd47' or 'd48' or 'd49'
		_D4x_ = f'D{self._4x}' # 'D47' or 'D48' or 'D49'

		# arrays of δ4x and Δ4x values
		d4x = np.array([_[_d4x_] for _ in self])
		D4x_raw = np.array([_[f'{_D4x_}raw'] for _ in self])

		# within-session timetag, used to model drifts of a/b/c via a2/b2/c2
		t = np.array([_['t'] for _ in self])

		# unknowns = samples not in fixed nor loose anchors
		unknowns = {
			s: (0.5, 2.)
			for s in self.samples
			if s not in self.RMs
		}

		# list of sessions:
		sessions = [s for s in self.sessions]
		n_sessions = len(sessions)
		# bidirectional session search:
		session_search = {s:k for k,s in enumerate(sessions)} | {k:s for k,s in enumerate(sessions)}
		# session index for all analyses:
		session_idx = np.array([session_search[_['Session']] for _ in self])

		# list of samples:
		samples = [s for s in self.samples]
		n_samples = len(samples)
		# bidirectional sample search:
		sample_search = {s:k for k,s in enumerate(samples)} | {k:s for k,s in enumerate(samples)}
		# sample index for all analyses:
		sample_idx = np.array([sample_search[_['Sample']] for _ in self])

		#### TRANSLATE PER-SESSION DRIFT FLAGS INTO CONSTRAINTS ####

		# Sessions with scrambling_drift, slope_drift and/or wg_drift set to False get
		# their corresponding a2/b2/c2 element constrained to '0.0' via the same mechanism
		# as any other constraint
		drift_constraints = {}
		for session in sessions:
			if not self.sessions[session]['scrambling_drift']:
				drift_constraints[f"a2['{session}']"] = '0.0'
			if not self.sessions[session]['slope_drift']:
				drift_constraints[f"b2['{session}']"] = '0.0'
			if not self.sessions[session]['wg_drift']:
				drift_constraints[f"c2['{session}']"] = '0.0'

		# Merge drift-derived constraints with user-specified ones
		constraints = drift_constraints | constraints

		#### HELPERS FOR SESSION SIGMA GROUPS ####

		if sigma_session_groups is None:
			# default: one single group containing every session
			# (equivalent to the original single-scalar-sigma behavior)
			sigma_session_groups = [sessions]
		else:
			# copy the input so we never mutate the caller's list/sublists in place
			sigma_session_groups = [list(group) for group in sigma_session_groups]

		# validate group contents:
		# every listed session must be real, and no session may be assigned to more than one group
		sessions_already_grouped = set()
		for group in sigma_session_groups:
			for session in group:
				if session not in sessions:
					raise ValueError(f"sigma_session_groups: unknown session '{session}'")
				if session in sessions_already_grouped:
					raise ValueError(f"sigma_session_groups: session '{session}' appears in more than one group")
				sessions_already_grouped.add(session)

		# silently collect any session absent from all groups into one extra group
		ungrouped_sessions = [session for session in sessions if session not in sessions_already_grouped]
		if ungrouped_sessions:
			sigma_session_groups.append(ungrouped_sessions)

		# total number of distinct sigma values to estimate
		n_sigma_groups = len(sigma_session_groups)

		# map each session name to its sigma-group index
		session_to_sigma_group = {}
		for group_index, group in enumerate(sigma_session_groups):
			for session in group:
				session_to_sigma_group[session] = group_index

		# array mapping each *position* in `sessions` (0..n_sessions-1) to its sigma-group index;
		# used to broadcast the per-group free sigma values onto the full `sessions`-dims vector
		sigma_group_of_session = np.array([session_to_sigma_group[session] for session in sessions])

		#### HELPERS FOR PARSING CONSTRAINTS ####

		# regex pattern, matches a[0], a['Session_01'], a["Session_01"]...
		_INDEX_PATTERN = re.compile(r"(\w+)\[(?:'([^']+)'|\"([^\"]+)\"|(\d+))\]")

		# functions users are allowed to reference in constraint expressions
		_ALLOWED_FUNCS = {
			# 'sqrt': pt.sqrt,
			# 'exp': pt.exp,
			# 'log': pt.log,
			# 'abs': pt.abs,
		}

		# base namespace: sympy's own names (Symbol, Integer, Rational, etc.), with
		# Python builtins stripped out so eval'd expressions can't reach __import__ etc.
		_SYMPY_NAMESPACE = {}
		exec("from sympy import *", _SYMPY_NAMESPACE)
		_SYMPY_NAMESPACE.pop("__builtins__", None)

		# which coordinate each constrainable variable is indexed along
		variable_dim = {
			'a':  'sessions',
			'b':  'sessions',
			'c':  'sessions',
			'a2': 'sessions',
			'b2': 'sessions',
			'c2': 'sessions',
			_D4x_: 'samples',
		}

		# lookup tables for each coordinate
		# {coord -> {label -> index}}
		label_to_pos = {
			'sessions': {s: k for k, s in enumerate(sessions)},
			'samples': {s: k for k, s in enumerate(samples)},
		}

		def _resolve_ref(
			base_name,
			label_or_index,
		):
			'''Translate and validate a parsed base_name + label/index into (base_name, index)'''

			if base_name not in variable_dim:
				raise ValueError(f"Constraints: unknown variable '{base_name}'. Must be one of {sorted(variable_dim)}.")

			dim = variable_dim[base_name]

			if label_or_index.isdigit():
				pos = int(label_or_index)
				if not (0 <= pos < len(label_to_pos[dim])):
					raise ValueError(f"Constraints: index {pos} out of range for '{base_name}' (dim '{dim}').")
			else:
				if label_or_index not in label_to_pos[dim]:
					raise ValueError(f"Constraints: label '{label_or_index}' not found in dim '{dim}' (variable '{base_name}').")
				pos = label_to_pos[dim][label_or_index]
			return base_name, pos

		def _parse_single_ref(ref_str):
			'''Parse a standalone reference string, e.g. "c['Session_02']", used for constraint keys.'''
			match = _INDEX_PATTERN.fullmatch(ref_str.strip())
			if not match:
				raise ValueError(f"Constraints: invalid reference syntax '{ref_str}' (expected e.g. \"a['Session_01']\" or \"b[1]\")")

			base_name = match.group(1)

			# pick out whichever of the three alternative index-formats actually matched,
			# since the regex has three mutually exclusive capture groups for the bracket content.
			label_or_index = match.group(2) or match.group(3) or match.group(4)

			return _resolve_ref(base_name, label_or_index)

		def _safe_parse(expr_str):
			'''Parse an algebraic expression string into sympy, without exposing Python builtins.'''
			return parse_expr(
				expr_str,
				transformations = standard_transformations,
				global_dict = _SYMPY_NAMESPACE,
				local_dict = {},
			)

		def _preprocess_expression(expr_str):
			'''Replace every 'name[...]' reference in expr_str with a safe placeholder symbol.
			Returns (rewritten_string, {placeholder: (base_name, position)}).'''

			ref_map = {}

			def _replace(match):

				base_name = match.group(1)
				label_or_index = match.group(2) or match.group(3) or match.group(4) # see explanation above

				# translate and validate parsed base_name + label/index into (base_name, index):
				base_name, pos = _resolve_ref(base_name, label_or_index)

				placeholder = f'{base_name}__{pos}'

				# update ref_map
				ref_map[placeholder] = (base_name, pos)

				# update ref_map
				return placeholder

			rewritten = _INDEX_PATTERN.sub(_replace, expr_str)

			return rewritten, ref_map

		#### PARSE ALL CONSTRAINTS ####

		parsed_constraints = {}
		# Target format = {
		#     (base_name, pos) -> {
		#         'expr': sympy expr,
		#         'symbols': [...],
		#         'ref_map': {...}
		#     }
		# }

		for target_str, expr_str in constraints.items():

			target = _parse_single_ref(target_str) # -> (base_name, pos)

			# Raise error if target is a fixed anchor
			if target[0] == _D4x_ and samples[target[1]] in self.fixed_RMs:
				raise ValueError(
					f"Constraints: '{target_str}' cannot be constrained because '{samples[target[1]]}' is a fixed anchor."
				)

			# Raise error if target has already been constrained
			if target in parsed_constraints:
				raise ValueError(f"Constraints: '{target_str}' cannot be the target of more than one constraint.")

			# Preprocess the constraint expression
			rewritten_expr, ref_map = _preprocess_expression(expr_str)

			# Parse the preprocessed expression
			try:
				sympy_expr = _safe_parse(rewritten_expr)
			except Exception as e:
				raise ValueError(f"Constraints: could not parse expression '{expr_str}' for target '{target_str}' ({e}).")

			# Check that no symbol remains undefined
			free_symbol_names = sorted(str(s) for s in sympy_expr.free_symbols)
			unresolved = [n for n in free_symbol_names if n not in ref_map]
			if unresolved:
				raise ValueError(f"Constraints: unrecognized term(s) {unresolved} in expression '{expr_str}'.")

			# Update parsed_constraints
			parsed_constraints[target] = {
				'expr': sympy_expr,
				'symbols': free_symbol_names,
				'ref_map': ref_map,
			}

		# Make a note of all constrained keys
		constrained_keys = set(parsed_constraints)

		# lookup table for which positions of a given variable name are constrained:
		def _constrained_positions(base_name):
			'''Positions of `base_name` that are the target of a constraint.'''
			return {pos for (b, pos) in constrained_keys if b == base_name}

		#### TOPOLOGICALLY ORDER CONSTRAINTS (DEPENDENCIES RESOLVED BEFORE DEPENDENTS) ####

		resolved_order = []
		visiting = set()

		def _visit(key, chain):
			'''Depth-first traversal building a dependency-respecting resolution order; raises on cycles.'''
			if key in resolved_order:
				return
			if key in visiting:
				cycle = ' -> '.join(f'{b}[{samples[p] if b == _D4x_ else sessions[p]}]' for b, p in chain + [key])
				raise ValueError(f'constraints: circular dependency detected ({cycle})')
			visiting.add(key)
			for symbol_name in parsed_constraints[key]['symbols']:
				dep = parsed_constraints[key]['ref_map'][symbol_name]
				if dep in constrained_keys:
					_visit(dep, chain + [key])
			visiting.discard(key)
			resolved_order.append(key)

		for key in constrained_keys:
			_visit(key, [])
		# at this point, resolved_order should be populated in the correct order

		#### GENERIC BUILDER FOR PARTIALLY-FREE VECTORS (USED FOR a, b, c, a2, b2, c2) ####

		def _build_vector(base_name, dist_fn, size, dims, **dist_kwargs):
			'''
			Build `size` slots for `base_name`. Positions not targeted by a constraint are
			drawn from dist_fn; constrained positions are left as None, to be filled in
			later once their defining expression can be evaluated.
			Returns (free_rv_or_None, slots).
			'''
			constrained_positions = _constrained_positions(base_name)
			free_positions = [i for i in range(size) if i not in constrained_positions]

			# subset any per-position kwarg (e.g. `mu`) down to the free positions only
			free_kwargs = {}
			for key, value in dist_kwargs.items():
				if hasattr(value, '__len__') and len(value) == size:
					free_kwargs[key] = [value[i] for i in free_positions]
				else:
					free_kwargs[key] = value

			free_name = base_name if not constrained_positions else f'{base_name}_free'
			free_dims = dims if not constrained_positions else None

			slots = [None] * size
			free_rv = None
			if free_positions:
				free_rv = dist_fn(free_name, shape = len(free_positions), dims = free_dims, **free_kwargs)
				for pos_in_free, pos in enumerate(free_positions):
					slots[pos] = free_rv[pos_in_free]

			return free_rv, slots

		def _finalize_vector(base_name, free_rv, slots, dims):
			'''Register the completed vector: a Deterministic if any position was constraint-derived,
			otherwise the free RV itself (unchanged from the non-constrained code path).'''
			if _constrained_positions(base_name):
				return pm.Deterministic(base_name, pt.stack(slots), dims = dims)
			return free_rv

		#### MODEL ####

		with pm.Model(coords = {'sessions': sessions, 'samples': samples}) as model:

			a_free, a_slots = _build_vector(
				base_name = 'a',
				dist_fn = pm.Uniform,
				size = n_sessions,
				dims = 'sessions',
				lower = 0.1,
				upper = 1.5,
			)

			b_free, b_slots = _build_vector(
				base_name = 'b',
				dist_fn = pm.Normal,
				size = n_sessions,
				dims = 'sessions',
				mu = [0. for session in self.sessions],
				sigma = 0.1,
			)

			c_free, c_slots = _build_vector(
				base_name = 'c',
				dist_fn = pm.Normal,
				size = n_sessions,
				dims = 'sessions',
				mu = [0.9 for session in self.sessions],
				sigma = 2,
			)

			# a2/b2/c2: per-session drift rates (multiplying `t`), analogous to a/b/c but
			# centered on zero, since "no drift" is the default expectation.
			a2_free, a2_slots = _build_vector(
				base_name = 'a2',
				dist_fn = pm.Normal,
				size = n_sessions,
				dims = 'sessions',
				mu = [0. for session in self.sessions],
				sigma = 0.01,
			)

			b2_free, b2_slots = _build_vector(
				base_name = 'b2',
				dist_fn = pm.Normal,
				size = n_sessions,
				dims = 'sessions',
				mu = [0. for session in self.sessions],
				sigma = 0.05,
			)

			c2_free, c2_slots = _build_vector(
				base_name = 'c2',
				dist_fn = pm.Normal,
				size = n_sessions,
				dims = 'sessions',
				mu = [0. for session in self.sessions],
				sigma = 0.5,
			)

			# one free Δ4x sigma value per group defined by sigma_session_groups,
			# shared by every session within its group
			sigma_group = pm.HalfNormal('sigma_group', sigma = 0.2, shape = n_sigma_groups)

			# broadcast each group's sigma value onto every session belonging to that group,
			# giving a full-length vector indexed like `a`, `b`, `c` (dims = 'sessions');
			# sigma_group_of_session[k] is the group index of sessions[k]
			sigma = pm.Deterministic('sigma', sigma_group[sigma_group_of_session], dims = 'sessions')

			# D4x: constants for fixed anchors, free Normals for the rest,
			# except positions targeted by a constraint, left as None for now
			D4x_slots = [None] * n_samples

			for i, sample in enumerate(samples):
				if (_D4x_, i) in parsed_constraints:
					continue   # filled in during constraint resolution below
				s = pf(sample)
				if sample in self.fixed_RMs:
					D4x_slots[i] = pt.constant(self.fixed_RMs[sample], name = f'D4x_{s}')
				else:
					mu, sig = (self.loose_RMs | unknowns)[sample]
					D4x_slots[i] = pm.Normal(f'D4x_{s}', mu = mu, sigma = sig)

			# Resolve constrained positions in the correct dependency order
			# (a2/b2/c2 slots added so constraints can reference/target them too)
			slots_by_var = {
				'a' :  a_slots, 'b' :  b_slots, 'c' :  c_slots,
				'a2': a2_slots, 'b2': b2_slots, 'c2': c2_slots,
				_D4x_: D4x_slots,
			}

			for base_name, pos in resolved_order:
				info = parsed_constraints[(base_name, pos)]
				symbol_objs = [sp.Symbol(n) for n in info['symbols']]
				numeric_func = sp.lambdify(symbol_objs, info['expr'], modules = [_ALLOWED_FUNCS, 'numpy'])

				arg_values = []
				for symbol_name in info['symbols']:
					dep_base, dep_pos = info['ref_map'][symbol_name]
					value = slots_by_var[dep_base][dep_pos]
					if value is None:
						raise ValueError(
							f'constraints: dependency {dep_base}[{dep_pos}] (needed for {base_name}[{pos}]) was not yet resolved (internal ordering error)'
						)
					arg_values.append(value)

				slots_by_var[base_name][pos] = numeric_func(*arg_values)

			# Finalize a, b, c, a2, b2, c2 (Deterministic if constrained, free RV otherwise)
			a = _finalize_vector('a', a_free, a_slots, 'sessions')
			b = _finalize_vector('b', b_free, b_slots, 'sessions')
			c = _finalize_vector('c', c_free, c_slots, 'sessions')
			a2 = _finalize_vector('a2', a2_free, a2_slots, 'sessions')
			b2 = _finalize_vector('b2', b2_free, b2_slots, 'sessions')
			c2 = _finalize_vector('c2', c2_free, c2_slots, 'sessions')

			# Finalize D4x (Deterministic if constrained, free RV otherwise)
			D4x = pm.Deterministic(_D4x_, pt.stack(D4x_slots), dims = 'samples')
			D4x_idx = {v: k for k, v in enumerate(samples)}

			# Convenience variable (WG composition, computed from a,c)
			# NB: calculated as -c/a (drift-free), corresponding to t=0
			D4x_wg = pm.Deterministic(f'D{self._4x}_wg', -c/a)

			# Build predicted raw D4x values (observations)
			mu = (
				(a[session_idx] + a2[session_idx] * t) * D4x[sample_idx]
				+ (b[session_idx] + b2[session_idx] * t) * d4x
				+ (c[session_idx] + c2[session_idx] * t)
			)

			# sigma is per-session (grouped) but indexed by session, just like a, b, c,
			# and rescaled by |a + a2*t| to convert from corrected to raw Δ4x noise
			# (abs() guards against the effective slope crossing zero when a2 is free).
			pm.Normal(
				'D4xraw',
				mu = mu,
				sigma = sigma[session_idx] * pt.abs(a[session_idx] + a2[session_idx] * t),
				observed = D4x_raw,
			)

			idata = pm.sample(
				**(default_mcmc_sample_kw | mcmc_sample_kw)
			)

		_posterior = idata.posterior

		self.standardization['bayes'] = dict(method = 'bayes')
		self.standardization['latest'] = self.standardization['bayes']

		self.standardization['bayes']['idata'] = idata
		if len(sigma_session_groups) == 1:
			pdf = _posterior['sigma'][:,:,0].values.reshape(-1)
			self.standardization['bayes']['sigma'] = uncertainties.ufloat(
				pdf.mean(),
				pdf.std(ddof = 1)
			)

		with warnings.catch_warnings():
			warnings.filterwarnings(
				'ignore',
				category = RuntimeWarning,
				message = 'invalid value encountered in scalar divide',
			)

			self.standardization['bayes']['summary'] = az.summary(
				idata,
				var_names = ['sigma', 'a', 'b', 'c', 'a2', 'b2', 'c2', f'D{self._4x}_wg', f'D{self._4x}'],
				round_to = 9,
				ci_kind = 'eti',
				ci_prob=0.95,
				# coords={'samples': [s for s in samples if s not in fixed_anchors]}
			)

		# populate self.standardization['bayes']['sessions']
		self.standardization['bayes']['sessions'] = {}
		S = self.standardization['bayes']['sessions']

		for session in self.sessions:

			params = ['a', 'b', 'c', 'a2', 'b2', 'c2', 'sigma']
			session_index = session_search[session]

			draws = np.array([_posterior[p][:,:,session_index].values.reshape(-1) for p in params])
			uparams = uncertainties.correlated_values(
				draws.mean(1),
				np.cov(draws),
			)

			S[session] = {}

			for p_index, p in enumerate(params):
				S[session][p] = uparams[p_index]
				S[session][f'pdf_{p}'] = draws[p_index,:]
				S[session][f'95CL_{p}'] = float(
					np.quantile(np.abs(draws[p_index,:] - uparams[p_index].n), 0.95)
				)

			S[session]['Np'] = 3 + sum([
				self.sessions[session][_]
				for _ in [
					'scrambling_drift',
					'slope_drift',
					'wg_drift',
				]
			])

		# populate self.standardization['bayes']['samples']
		D4x = f'D{self._4x}'
		self.standardization['bayes']['samples'] = {}
		S = self.standardization['bayes']['samples']

		for s in self.fixed_anchors:
			S[s] = {}

		for s in (self.loose_anchors | self.unknowns):
			S[s] = {f'pdf_{D4x}': _posterior[D4x].sel(samples = s).values.reshape(-1)}

		# Caution: we are redefining samples here. Why? Do we really need to?
		samples = [s for s in S]
		sample_index = [D4x_idx[s] for s in samples]

		draws = np.array([_posterior[D4x][:,:,i].values.reshape(-1) for i in sample_index])
		uD4x = uncertainties.correlated_values(
			draws.mean(1),
			np.cov(draws),
		)
		for k,s in enumerate(samples):
			S[s][D4x] = uD4x[k]
			if sample in self.fixed_RMs:
				S[s][f'95CL_{D4x}'] = 0.
			else:
				S[s][f'95CL_{D4x}'] = float(
					np.quantile(np.abs(draws[k,:] - draws[k,:].mean()), 0.95)
				)

		for r in self:
			s = r["Session"]
			a = self.standardization['bayes']['sessions'][s]['a'].n
			b = self.standardization['bayes']['sessions'][s]['b'].n
			c = self.standardization['bayes']['sessions'][s]['c'].n
			a2 = self.standardization['bayes']['sessions'][s]['a2'].n
			b2 = self.standardization['bayes']['sessions'][s]['b2'].n
			c2 = self.standardization['bayes']['sessions'][s]['c2'].n
			r[D4x] = (r[f'{D4x}raw'] - c - b * r[f'd{self._4x}'] - c2 * r['t'] - b2 * r['t'] * r[f'd{self._4x}']) / (a + a2 * r['t'])

	@make_verbal
	def table_of_least_squares_vs_bayesian_results(
		self,
		dir = 'output',
		filename = None,
		save_to_file = True,
		print_out = True,
		output = None,
	):
		D4x = f'D{self._4x}'
		out = [['Sample','N', f'D{self._4x} (LS)','SE','95% CL', f'D{self._4x} (Bayes)','SE','95% CL', 'difference']]
		pooled = self.standardization['pooled']['samples']
		bayes = self.standardization['bayes']['samples']
		for sample in self.standardization['bayes']['samples']:
			out += [[
				f"{sample}",
				f"{self.samples[sample]['N']}",
				f"{pooled[sample][D4x].n:.4f}",
				f"{pooled[sample][D4x].s:.4f}" if sample in self.unknowns else '',
				f"± {pooled[sample][f'95CL_{D4x}']:.4f}" if sample in self.unknowns else '',
				f"{bayes[sample][D4x].n:.4f}",
				f"{bayes[sample][D4x].s:.4f}" if sample not in self.fixed_RMs else '',
				f"± {bayes[sample][f'95CL_{D4x}']:.4f}" if sample not in self.fixed_RMs else '',
				f"{round(bayes[sample][D4x].n - pooled[sample][D4x].n, 4) + 0.0:.4f}",
			]]
		if save_to_file:
			if not os.path.exists(dir):
				os.makedirs(dir)
			if filename is None:
				filename = f'D{self._4x}_ls_vs_bayes.csv'
			with open(f'{dir}/{filename}', 'w') as fid:
				fid.write(make_csv(out))
		if print_out:
			self.msg('\n'+pretty_table(out))
		if output == 'raw':
			return out
		elif output == 'pretty':
			return pretty_table(out)



	def plot_least_squares_vs_bayesian_results(
		self,
		target = 'samples',
		figsize = None,
		columns = 1,
		left_margin = 0.2,
		right_margin = 0.2,
		top_margin = 0.5,
		bottom_margin = 0.7,
		bayes_on_top = True,
		cell_width = 6,
		cell_height = 0.5,
		dir = 'output',
		savefig = True,
		filename = '',
		dpi = 100,
		weak_anchor_color = (1, 0.25, 0),
		ls_color = (0.5, 0.5, 0.5),
		bayes_color = (0, 0.75, 1),
	):

		from scipy.stats import gaussian_kde
		from coloraide import Color

		samples = [s for s in self.samples if s not in self.fixed_RMs]
		N = len(samples)
		lines  = N//columns
		if target == 'samples':
			if figsize is None:
				figsize = (
					columns * cell_width + (columns + 1) * (left_margin + right_margin)/2,
					lines * cell_height + (lines + 1) * (top_margin + bottom_margin)/2,
				)
			fig = ppl.figure(figsize = figsize)
			ppl.subplots_adjust(
				left_margin/figsize[0],
				bottom_margin/figsize[1],
				1 - right_margin/figsize[0],
				1 - top_margin/figsize[1],
				(left_margin + right_margin)/2/cell_width,
				(top_margin + bottom_margin)/2/cell_height,
			)
			axs = [ppl.subplot(lines, columns, _+1) for _ in range(N)]
			xmin, xmax = 1000, -1000
			for sample in samples:
				x = self.standardization['bayes']['samples'][sample][f'pdf_D{self._4x}']
				x0 = x.min() - 0.02
				x1 = x.max() + 0.02
				xmin = min(xmin, x0)
				xmax = max(xmax, x1)

			for sample, ax in zip(samples, axs):
				ppl.sca(ax)
				ppl.yticks([])
				ppl.title(sample, weight = 'bold', size = 10)
				xi = np.linspace(xmin, xmax, 1001)
				x = self.standardization['bayes']['samples'][sample][f'pdf_D{self._4x}']
				yi = gaussian_kde(x).evaluate(xi)
				_color = Color('srgb', bayes_color)
				kw = dict(
					ec = Color(_color).mix('white', 0.2, space = 'srgb'),
					fc = Color(_color).set('alpha', 0.3),
					lw = 1,
					zorder = 5,
				)
				ax.fill_between(
					xi,
					yi*(1 if bayes_on_top else 0),
					yi*(0 if bayes_on_top else -1),
					**kw,
				)
				xm = self.standardization['bayes']['samples'][sample][f'D{self._4x}'].n
				ha = 'left' if (xm - xmin) > (xmax - xm) else 'right'
				x = xmin if (xm - xmin) > (xmax - xm) else xmax
				ax.text(
					x,
					0,
					'  Bayesian posterior\n' if ha == 'left' else 'Bayesian posterior  \n',
					ha = ha,
					va = 'center',
					color = _color,
					linespacing = 1.9,
				)

				if sample in self.loose_RMs:
					_color = Color('srgb', weak_anchor_color)
					kw = dict(
						ec = Color(_color).mix('white', 0.3, space = 'srgb'),
						fc = Color(_color).set('alpha', 0.15),
						lw = 1,
						zorder = 2,
					)
					yi = yi.max() * np.exp(-0.5 * ((xi-self.loose_RMs[sample][0])/self.loose_RMs[sample][1])**2)
					ax.fill_between(
						xi,
						yi*(1 if bayes_on_top else 0),
						yi*(0 if bayes_on_top else -1),
						**kw,
					)
					xm = self.loose_RMs[sample][0]
					ha = 'left' if (xm - xmin) > (xmax - xm) else 'right'
					x = xmin if (xm - xmin) > (xmax - xm) else xmax
					ax.text(
						x,
						0,
						' '*33 + '(and prior)\n' if ha == 'left' else '(and prior)' + ' '*33 + '\n',
						ha = ha,
						va = 'center',
						color = _color,
						linespacing = 1.9,
					)
				if sample in self.unknowns:
					_color = Color('srgb', ls_color)
					kw = dict(
						ec = Color(_color).mix('white', 0.5, space = 'srgb'),
						fc = Color(_color).set('alpha', 0.2),
						lw = 1,
						zorder = 4,
					)
					mu = self.samples[sample][f'D{self._4x}']
					sigma = self.samples[sample][f'SE_D{self._4x}']
					yi = yi.max() * np.exp(-0.5 * ((xi - mu)/sigma)**2)
					ax.fill_between(
						xi,
						yi*(0 if bayes_on_top else 1),
						yi*(-1 if bayes_on_top else 0),
						**kw,
					)
					xm = self.standardization['pooled']['samples'][sample][f'D{self._4x}'].n
					ha = 'left' if (xm - xmin) > (xmax - xm) else 'right'
					x = xmin if (xm - xmin) > (xmax - xm) else xmax
					ax.text(
						x,
						0,
						'\n  Least squares' if ha == 'left' else '\nLeast squares  ',
						ha = ha,
						va = 'center',
						color = _color,
						linespacing = 1.9,
					)

			for ax in axs:
				ppl.sca(ax)
				ppl.grid(alpha = 0.2)
				ppl.axis([xmin, xmax, None, None])

			ppl.xlabel(f'Δ{self._4x} [‰]')

			if savefig:
				if not os.path.exists(dir):
					os.makedirs(dir)
				if filename is None:
					return fig
				elif filename == '':
					filename = f'D{self._4x}_ls_vs_bayes.pdf'
				ppl.savefig(f'{dir}/{filename}', dpi = dpi)
				ppl.close(fig)
			else:
				return fig

	def plot_single_bayesian_session(self,
		session,
		kw_plot_fixed_anchors = dict(ls='None', marker='x', mec=(.75, 0, 0), mew = .75, ms = 4),
		kw_plot_weak_anchors = dict(ls='None', marker='x', mec=(1, 0, .25), mew = .75, ms = 4),
		kw_plot_unknowns = dict(ls='None', marker='x', mec=(0, 0, .75), mew = .75, ms = 4),
		kw_plot_fixed_anchor_avg = dict(ls='-', marker='None', color=(.75, 0, 0), lw = 2, alpha = 1/3),
		kw_fill_weak_anchor_avg = dict(color=(1, 0, .25), lw = 0, alpha = 1/3),
		kw_fill_unknown_avg = dict(color=(0, 0, .75), lw = 0, alpha = 0.2),
		kw_contour_error = dict(colors = [[0, 0, 0]], alpha = .5, linewidths = 0.75),
		xylimits = 'free', # | 'constant'
		x_label = None,
		y_label = None,
		error_contour_interval = 'auto',
		fig = 'new',
		):
		'''
		Generate plot for a single session after Bayesian standardization
		'''
		if x_label is None:
			x_label = f'δ$_{{{self._4x}}}$ (‰)'
		if y_label is None:
			y_label = f'Δ$_{{{self._4x}}}$ (‰)'

		out = _SessionPlot()

		fixed_anchors = [a for a in self.fixed_anchors if [r for r in self.sessions[session]['data'] if r['Sample'] == a]]
		weak_anchors   = [a for a in self.loose_anchors if [r for r in self.sessions[session]['data'] if r['Sample'] == a]]
		unknowns       = [u for u in self.unknowns if [r for r in self.sessions[session]['data'] if r['Sample'] == u]]

		fixed_anchors_d = [r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in fixed_anchors]
		fixed_anchors_D = [r[f'D{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in fixed_anchors]
		weak_anchors_d   = [r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in weak_anchors]
		weak_anchors_D   = [r[f'D{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in weak_anchors]
		unknowns_d       = [r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in unknowns]
		unknowns_D       = [r[f'D{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in unknowns]

		fixed_anchor_avg = (
			np.array([
				[
					np.min([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) - 0.5,
					np.max([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) + 0.5,
				]
				for sample in fixed_anchors
			]).T,
			np.array([
				self.bare_RMs[sample] + np.array([0, 0])
				for sample in fixed_anchors
			]).T
		)

		weak_anchor_avg = (
			np.array([
				[
					np.min([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) - 0.5,
					np.min([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) - 0.5,
					np.max([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) + 0.5,
					np.max([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) + 0.5,
					np.min([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) - 0.5,
				]
				for sample in weak_anchors
			]).T,
			np.array([
				self.standardization['bayes']['samples'][sample][f'D{self._4x}'].n
				+ np.array([-1, 1, 1, -1, -1])/2
				* self.standardization['bayes']['samples'][sample][f'95CL_D{self._4x}']
				for sample in weak_anchors
			]).T
		)

		unknown_avg = (
			np.array([
				[
					np.min([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) - 0.5,
					np.min([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) - 0.5,
					np.max([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) + 0.5,
					np.max([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) + 0.5,
					np.min([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) - 0.5,
				]
				for sample in unknowns
			]).T,
			np.array([
				self.standardization['bayes']['samples'][sample][f'D{self._4x}'].n
				+ np.array([-1, 1, 1, -1, -1])
				* self.standardization['bayes']['samples'][sample][f'95CL_D{self._4x}']
				for sample in unknowns
			]).T
		)

		if fig == 'new':
			out.fig = ppl.figure(figsize = (6,6))
			ppl.subplots_adjust(.1,.1,.9,.9)

		out.fixed_anchor_analyses, = ppl.plot(
			fixed_anchors_d,
			fixed_anchors_D,
			**kw_plot_fixed_anchors)
		out.weak_anchor_analyses, = ppl.plot(
			weak_anchors_d,
			weak_anchors_D,
			**kw_plot_weak_anchors)
		out.unknown_analyses, = ppl.plot(
			unknowns_d,
			unknowns_D,
			**kw_plot_unknowns)
		out.fixed_anchor_avg = ppl.plot(
			*fixed_anchor_avg,
			**kw_plot_fixed_anchor_avg)
		out.weak_anchor_avg = ppl.fill(
			*weak_anchor_avg,
			**kw_fill_weak_anchor_avg)
		out.unknown_avg = ppl.fill(
			*unknown_avg,
			**kw_fill_unknown_avg)

		if xylimits == 'constant':
			x = [r[f'd{self._4x}'] for r in self]
			y = [r[f'D{self._4x}'] for r in self]
			x1, x2, y1, y2 = np.min(x), np.max(x), np.min(y), np.max(y)
			w, h = x2-x1, y2-y1
			x1 -= w/20
			x2 += w/20
			y1 -= h/20
			y2 += h/20
			ppl.axis([x1, x2, y1, y2])
		elif xylimits == 'free':
			x1, x2, y1, y2 = ppl.axis()
		else:
			x1, x2, y1, y2 = ppl.axis(xylimits)

		if error_contour_interval != 'none':
			xi, yi = np.linspace(x1, x2), np.linspace(y1, y2)
			XI,YI = np.meshgrid(xi, yi)
			# SI = np.array([[self.standardization_error(session, x, y) for x in xi] for y in yi])
			_a = self.standardization['bayes']['sessions'][session]['pdf_a']
			_b = self.standardization['bayes']['sessions'][session]['pdf_b']
			_c = self.standardization['bayes']['sessions'][session]['pdf_c']
			_ZI = _a.mean()*YI + _b.mean()*XI + _c.mean()
			_YI = (_ZI[:,:,None] - _b*XI[:,:,None] - _c) / _a
			SI = (_YI).std(axis = -1, ddof = 1)
			if error_contour_interval == 'auto':
				rng = np.max(SI) - np.min(SI)
				if rng <= 0.01:
					cinterval = 0.001
				elif rng <= 0.03:
					cinterval = 0.004
				elif rng <= 0.1:
					cinterval = 0.01
				elif rng <= 0.3:
					cinterval = 0.03
				elif rng <= 1.:
					cinterval = 0.1
				else:
					cinterval = 0.5
			else:
				cinterval = error_contour_interval

			cval = np.arange(np.ceil(SI.min() / .001) * .001, np.ceil(SI.max() / .001 + 1) * .001, cinterval)
			out.contour = ppl.contour(XI, YI, SI, cval, **kw_contour_error)
			out.clabel = ppl.clabel(out.contour)
			contour = (XI, YI, SI, cval, cinterval)

		if fig == None:
			return {
			'anchors':anchors,
			'unknowns':unknowns,
			'anchors_d':anchors_d,
			'anchors_D':anchors_D,
			'unknowns_d':unknowns_d,
			'unknowns_D':unknowns_D,
			'anchor_avg':anchor_avg,
			'unknown_avg':unknown_avg,
			'contour':contour,
			}

		ppl.xlabel(x_label)
		ppl.ylabel(y_label)
		ppl.title(session, weight = 'bold')
		ppl.grid(alpha = .2)
		out.ax = ppl.gca()

		return out

	def plot_bayesian_sessions(self, dir = 'output', figsize = (8,8), filetype = 'pdf', dpi = 100):
		'''
		Generate Bayesian session plots and save them to disk.

		**Parameters**

		+ `dir`: the directory in which to save the plots
		+ `figsize`: the width and height (in inches) of each plot
		+ `filetype`: 'pdf' or 'png'
		+ `dpi`: resolution for PNG output
		'''
		if not os.path.exists(dir):
			os.makedirs(dir)

		for session in self.sessions:
			sp = self.plot_single_bayesian_session(session, xylimits = 'constant')
			ppl.savefig(f'{dir}/D{self._4x}_plot_{session}.{filetype}', **({'dpi': dpi} if filetype.lower() == 'png' else {}))
			ppl.close(sp.fig)

	def standardization_error(self, session, d4x, D4x, t = 0, target = 'latest'):
		'''
		Compute standardization error for a given session and
		(δ47, Δ47) composition.
		'''

		target = self._resolve_target(target)
		stdz = self.standardization[target]

		a = stdz['sessions'][session]['a']
		b = stdz['sessions'][session]['b']
		c = stdz['sessions'][session]['c']
		a2 = stdz['sessions'][session]['a2']
		b2 = stdz['sessions'][session]['b2']
		c2 = stdz['sessions'][session]['c2']

		x, y = D4x, d4x
		D4x_raw = (a * x + b * y + c + a2 * x * t + b2 * y * t + c2 * t).n
		sx = ((D4x_raw - b*y - b2*y*t - c - c2*t) / (a + a2*t)).s

		return sx

	@make_verbal
	def summary(self,
		dir = 'output',
		filename = None,
		save_to_file = True,
		print_out = True,
		):
		'''
		Print out an/or save to disk a summary of the standardization results.

		**Parameters**

		+ `dir`: the directory in which to save the table
		+ `filename`: the name to the csv file to write to
		+ `save_to_file`: whether to save the table to disk
		+ `print_out`: whether to print out the table
		'''

		out = []
		out += [['N samples (anchors + unknowns)', f"{len(self.samples)} ({len(self.anchors)} + {len(self.unknowns)})"]]
		out += [['N analyses (anchors + unknowns)', f"{len(self)} ({len([r for r in self if r['Sample'] in self.anchors])} + {len([r for r in self if r['Sample'] in self.unknowns])})"]]
		out += [['Repeatability of δ13C_VPDB', f"{1000 * self.repeatability['r_d13C_VPDB']:.1f} ppm"]]
		out += [['Repeatability of δ18O_VSMOW', f"{1000 * self.repeatability['r_d18O_VSMOW']:.1f} ppm"]]
		out += [[f'Repeatability of Δ{self._4x} (anchors)', f"{1000 * self.repeatability[f'r_D{self._4x}a']:.1f} ppm"]]
		out += [[f'Repeatability of Δ{self._4x} (unknowns)', f"{1000 * self.repeatability[f'r_D{self._4x}u']:.1f} ppm"]]
		out += [[f'Repeatability of Δ{self._4x} (all)', f"{1000 * self.repeatability[f'r_D{self._4x}']:.1f} ppm"]]
		out += [['Model degrees of freedom', f"{self.Nf}"]]
		out += [['Student\'s 95% t-factor', f"{self.standardization['latest']['t95']:.2f}"]]
		out += [['Standardization method', self.standardization_method]]

		if save_to_file:
			if not os.path.exists(dir):
				os.makedirs(dir)
			if filename is None:
				filename = f'D{self._4x}_summary.csv'
			with open(f'{dir}/{filename}', 'w') as fid:
				fid.write(make_csv(out))
		if print_out:
			self.msg('\n' + pretty_table(out, header = 0))


	@make_verbal
	def table_of_sessions(self,
		dir = 'output',
		filename = None,
		save_to_file = True,
		print_out = True,
		output = None,
		target = 'latest',
		):
		'''
		Print out an/or save to disk a table of sessions.

		**Parameters**

		+ `dir`: the directory in which to save the table
		+ `filename`: the name to the csv file to write to
		+ `save_to_file`: whether to save the table to disk
		+ `print_out`: whether to print out the table
		+ `output`: if set to `'pretty'`: return a pretty text table (see `pretty_table()`);
		    if set to `'raw'`: return a list of list of strings
		    (e.g., `[['header1', 'header2'], ['0.1', '0.2']]`)
		'''

		target = self._resolve_target(target)

		include_a2 = any([self.sessions[session]['scrambling_drift'] for session in self.sessions])
		include_b2 = any([self.sessions[session]['slope_drift'] for session in self.sessions])
		include_c2 = any([self.sessions[session]['wg_drift'] for session in self.sessions])

		out = [['Session','Na','Nu','d13Cwg_VPDB','d18Owg_VSMOW','r_d13C','r_d18O',f'r_D{self._4x}','a ± SE','1e3 x b ± SE','c ± SE']]
		if include_a2:
			out[-1] += ['a2 ± SE']
		if include_b2:
			out[-1] += ['b2 ± SE']
		if include_c2:
			out[-1] += ['c2 ± SE']
		for session in self.sessions:
			S = self.standardization[target]['sessions'][session]
			out += [[
				session,
				f"{self.sessions[session]['Na']}",
				f"{self.sessions[session]['Nu']}",
				f"{self.sessions[session]['d13Cwg_VPDB']:.3f}",
				f"{self.sessions[session]['d18Owg_VSMOW']:.3f}",
				f"{self.sessions[session]['r_d13C_VPDB']:.4f}",
				f"{self.sessions[session]['r_d18O_VSMOW']:.4f}",
				f"{self.sessions[session][f'r_D{self._4x}']:.4f}",
				f"{S['a'].n:.3f} ± {S['a'].s:.3f}",
				f"{1e3*S['b'].n:.3f} ± {1e3*S['b'].s:.3f}",
				f"{S['c'].n:.3f} ± {S['c'].s:.3f}",
				]]
			if include_a2:
				if self.sessions[session]['scrambling_drift']:
					out[-1] += [f"{S['a2'].n:.1e} ± {S['a2'].s:.1e}"]
				else:
					out[-1] += ['']
			if include_b2:
				if self.sessions[session]['slope_drift']:
					out[-1] += [f"{S['b2'].n:.1e} ± {S['b2'].s:.1e}"]
				else:
					out[-1] += ['']
			if include_c2:
				if self.sessions[session]['wg_drift']:
					out[-1] += [f"{S['c2'].n:.1e} ± {S['c2'].s:.1e}"]
				else:
					out[-1] += ['']

		if save_to_file:
			if not os.path.exists(dir):
				os.makedirs(dir)
			if filename is None:
				filename = f'D{self._4x}_sessions_{target}.csv'
			with open(f'{dir}/{filename}', 'w') as fid:
				fid.write(make_csv(out))
		if print_out:
			self.msg('\n' + pretty_table(out))
		if output == 'raw':
			return out
		elif output == 'pretty':
			return pretty_table(out)


	@make_verbal
	def table_of_analyses(
		self,
		dir = 'output',
		filename = None,
		save_to_file = True,
		print_out = True,
		output = None,
		):
		'''
		Print out an/or save to disk a table of analyses.

		**Parameters**

		+ `dir`: the directory in which to save the table
		+ `filename`: the name to the csv file to write to
		+ `save_to_file`: whether to save the table to disk
		+ `print_out`: whether to print out the table
		+ `output`: if set to `'pretty'`: return a pretty text table (see `pretty_table()`);
		    if set to `'raw'`: return a list of list of strings
		    (e.g., `[['header1', 'header2'], ['0.1', '0.2']]`)
		'''

		out = [['UID','Session','Sample']]
		extra_fields = [f for f in [('SampleMass','.2f'),('ColdFingerPressure','.1f'),('AcidReactionYield','.3f')] if f[0] in {k for r in self for k in r}]
		for f in extra_fields:
			out[-1] += [f[0]]
		out[-1] += ['d13Cwg_VPDB','d18Owg_VSMOW','d45','d46','d47','d48','d49','d13C_VPDB','d18O_VSMOW','D47raw','D48raw','D49raw',f'D{self._4x}']
		for r in self:
			out += [[f"{r['UID']}",f"{r['Session']}",f"{r['Sample']}"]]
			for f in extra_fields:
				out[-1] += [f"{r[f[0]]:{f[1]}}"]
			out[-1] += [
				f"{r['d13Cwg_VPDB']:.3f}",
				f"{r['d18Owg_VSMOW']:.3f}",
				f"{r['d45']:.6f}",
				f"{r['d46']:.6f}",
				f"{r['d47']:.6f}",
				f"{r['d48']:.6f}",
				f"{r['d49']:.6f}",
				f"{r['d13C_VPDB']:.6f}",
				f"{r['d18O_VSMOW']:.6f}",
				f"{r['D47raw']:.6f}",
				f"{r['D48raw']:.6f}",
				f"{r['D49raw']:.6f}",
				f"{r[f'D{self._4x}']:.6f}"
				]
		if save_to_file:
			if not os.path.exists(dir):
				os.makedirs(dir)
			if filename is None:
				filename = f'D{self._4x}_analyses.csv'
			with open(f'{dir}/{filename}', 'w') as fid:
				fid.write(make_csv(out))
		if print_out:
			self.msg('\n' + pretty_table(out))
		return out

	@make_verbal
	def covar_table(
		self,
		correl = False,
		dir = 'output',
		filename = None,
		save_to_file = True,
		print_out = True,
		output = None,
		):
		'''
		Print out, save to disk and/or return the variance-covariance matrix of D4x
		for all unknown samples.

		**Parameters**

		+ `dir`: the directory in which to save the csv
		+ `filename`: the name of the csv file to write to
		+ `save_to_file`: whether to save the csv
		+ `print_out`: whether to print out the matrix
		+ `output`: if set to `'pretty'`: return a pretty text matrix (see `pretty_table()`);
		    if set to `'raw'`: return a list of list of strings
		    (e.g., `[['header1', 'header2'], ['0.1', '0.2']]`)
		'''
		samples = sorted([u for u in self.unknowns])
		out = [[''] + samples]
		for s1 in samples:
			out.append([s1])
			for s2 in samples:
				if correl:
					out[-1].append(f'{self.sample_D4x_correl(s1, s2):.6f}')
				else:
					out[-1].append(f'{self.sample_D4x_covar(s1, s2):.8e}')

		if save_to_file:
			if not os.path.exists(dir):
				os.makedirs(dir)
			if filename is None:
				if correl:
					filename = f'D{self._4x}_correl.csv'
				else:
					filename = f'D{self._4x}_covar.csv'
			with open(f'{dir}/{filename}', 'w') as fid:
				fid.write(make_csv(out))
		if print_out:
			self.msg('\n'+pretty_table(out))
		if output == 'raw':
			return out
		elif output == 'pretty':
			return pretty_table(out)

	@make_verbal
	def table_of_samples(
		self,
		dir = 'output',
		filename = None,
		save_to_file = True,
		print_out = True,
		output = None,
		target = 'latest',
		):
		'''
		Print out, save to disk and/or return a table of samples.

		**Parameters**

		+ `dir`: the directory in which to save the csv
		+ `filename`: the name of the csv file to write to
		+ `save_to_file`: whether to save the csv
		+ `print_out`: whether to print out the table
		+ `output`: if set to `'pretty'`: return a pretty text table (see `pretty_table()`);
		    if set to `'raw'`: return a list of list of strings
		    (e.g., `[['header1', 'header2'], ['0.1', '0.2']]`)
		'''

		target = self._resolve_target(target)

		out = [[
			'Sample',
			'N',
			'd13C_VPDB',
			'd18O_VSMOW',
			f'D{self._4x}',
			'SE',
			'95% CL',
			'SD',
			# 'p_Levene',
		]]

		S = self.standardization[target]['samples']
		match target:
			case 'bayes':
				for sample in self.fixed_anchors:
					out += [[
						f"{sample}",
						f"{self.samples[sample]['N']}",
						f"{self.samples[sample]['d13C_VPDB']:.2f}",
						f"{self.samples[sample]['d18O_VSMOW']:.2f}",
						f"{S[sample][f'D{self._4x}'].n:.4f}",'','',
						f"{self.samples[sample][f'SD_D{self._4x}']:.4f}" if self.samples[sample]['N'] > 1 else '', ''
						]]
				for sample in (self.unknowns | self.loose_anchors):
					out += [[
						f"{sample}",
						f"{self.samples[sample]['N']}",
						f"{self.samples[sample]['d13C_VPDB']:.2f}",
						f"{self.samples[sample]['d18O_VSMOW']:.2f}",
						f"{S[sample][f'D{self._4x}'].n:.4f}",
						f"{S[sample][f'D{self._4x}'].s:.4f}",
						f"± {S[sample][f'95CL_D{self._4x}']:.4f}",
						f"{self.samples[sample][f'SD_D{self._4x}']:.4f}" if self.samples[sample]['N'] > 1 else '',
						# f"{self.samples[sample]['p_Levene']:.3f}" if self.samples[sample]['N'] > 2 else ''
						]]
			case 'pooled':
				for sample in self.anchors:
					out += [[
						f"{sample}",
						f"{self.samples[sample]['N']}",
						f"{self.samples[sample]['d13C_VPDB']:.2f}",
						f"{self.samples[sample]['d18O_VSMOW']:.2f}",
						f"{S[sample][f'D{self._4x}'].n:.4f}",'','',
						f"{self.samples[sample][f'SD_D{self._4x}']:.4f}" if self.samples[sample]['N'] > 1 else '', ''
						]]
				for sample in self.unknowns:
					out += [[
						f"{sample}",
						f"{self.samples[sample]['N']}",
						f"{self.samples[sample]['d13C_VPDB']:.2f}",
						f"{self.samples[sample]['d18O_VSMOW']:.2f}",
						f"{S[sample][f'D{self._4x}'].n:.4f}",
						f"{S[sample][f'D{self._4x}'].s:.4f}",
						f"± {S[sample][f'95CL_D{self._4x}']:.4f}",
						f"{self.samples[sample][f'SD_D{self._4x}']:.4f}" if self.samples[sample]['N'] > 1 else '',
						]]

		if save_to_file:
			if not os.path.exists(dir):
				os.makedirs(dir)
			if filename is None:
				filename = f'D{self._4x}_samples.csv'
			with open(f'{dir}/{filename}', 'w') as fid:
				fid.write(make_csv(out))
		if print_out:
			self.msg('\n'+pretty_table(out))
		if output == 'raw':
			return out
		elif output == 'pretty':
			return pretty_table(out)

	def _resolve_target(self, target):
		if target == 'latest':
			return self.standardization['latest']['method']
		return target

	def plot_sessions(
		self,
		target = 'latest',
		dir = 'output',
		figsize = (8,8),
		filetype = 'pdf',
		dpi = 100,
	):
		'''
		Generate session plots and save them to disk.

		**Parameters**

		+ `target`: which standardization method to target. May be:
		  - `'pooled'`: least-squares pooled regression method
		  - `'bayes': Bayesian regression method
		  - `None`: automatically target `'pooled'` or `'bayes'` depending on what is available,
		    and raise an error if both are available.
		+ `dir`: the directory in which to save the plots
		+ `figsize`: the width and height (in inches) of each plot
		+ `filetype`: 'pdf' or 'png'
		+ `dpi`: resolution for PNG output
		'''

		target = self._resolve_target(target)

		if not os.path.exists(dir):
			os.makedirs(dir)

		for session in self.sessions:
			match target:
				case 'pooled':
					sp = self.plot_single_pooled_session(session, xylimits = 'constant')
				case 'bayes':
					sp = self.plot_single_bayesian_session(session, xylimits = 'constant')
			ppl.savefig(f'{dir}/D{self._4x}_plot_{session}.{filetype}', **({'dpi': dpi} if filetype.lower() == 'png' else {}))
			ppl.close(sp.fig)

	@make_verbal
	def consolidate_samples(self, target = None):
		'''
		Compile various statistics for each sample.

		For each anchor sample:

		+ `D47` or `D48`: the nominal Δ4x value for this anchor, specified by `self.RMs`
		+ `SE_D47` or `SE_D48`: set to zero by definition

		For each unknown sample:

		+ `D47` or `D48`: the standardized Δ4x value for this unknown
		+ `SE_D47` or `SE_D48`: the standard error of Δ4x for this unknown

		For each anchor and unknown:

		+ `N`: the total number of analyses of this sample
		+ `SD_D47` or `SD_D48`: the “sample” (in the statistical sense) standard deviation for this sample
		+ `d13C_VPDB`: the average δ13C_VPDB value for this sample
		+ `d18O_VSMOW`: the average δ18O_VSMOW value for this sample (as CO2)
		'''

		target = self._resolve_target(target)
		_D4x_ = f'D{self._4x}'

		for sample in self.samples:

			for k in self.samples[sample]:
				if k != 'data':
					print(f'Popping {k} from {sample}.')
					self.samples[sample].pop(k)

			self.samples[sample]['N'] = len(self.samples[sample]['data'])

			if self.samples[sample]['N'] > 1:
				self.samples[sample][f'SD_D{self._4x}'] = stdev([r[f'D{self._4x}'] for r in self.samples[sample]['data']])

			self.samples[sample]['d13C_VPDB'] = np.mean([r['d13C_VPDB'] for r in self.samples[sample]['data']])
			self.samples[sample]['d18O_VSMOW'] = np.mean([r['d18O_VSMOW'] for r in self.samples[sample]['data']])

			self.samples[sample][_D4x_] = self.standardization[target]['samples'][sample][_D4x_].n
			self.samples[sample][f'SE_{_D4x_}'] = self.standardization[target]['samples'][sample][_D4x_].s
			self.samples[sample][f'95CL_{_D4x_}'] = self.standardization[target]['samples'][sample][f'95CL_{_D4x_}']

		for r in self:
			r[f'{_D4x_}_residual'] = r[f'D{self._4x}'] - self.standardization[target]['samples'][r['Sample']][f'D{self._4x}'].n

	def consolidate_sessions(self, target = 'latest'):
		'''
		Compute various statistics for each session.

		+ `Na`: Number of anchor analyses in the session
		+ `Nu`: Number of unknown analyses in the session
		+ `r_d13C_VPDB`: δ13C_VPDB repeatability of analyses within the session
		+ `r_d18O_VSMOW`: δ18O_VSMOW repeatability of analyses within the session
		+ `r_D47` or `r_D48`: Δ4x repeatability of analyses within the session
		+ `a`: scrambling factor
		+ `b`: compositional slope
		+ `c`: WG offset
		+ `SE_a`: Model stadard erorr of `a`
		+ `SE_b`: Model stadard erorr of `b`
		+ `SE_c`: Model stadard erorr of `c`
		+ `scrambling_drift` (boolean): whether to allow a temporal drift in the scrambling factor (`a`)
		+ `slope_drift` (boolean): whether to allow a temporal drift in the compositional slope (`b`)
		+ `wg_drift` (boolean): whether to allow a temporal drift in the WG offset (`c`)
		+ `a2`: scrambling factor drift
		+ `b2`: compositional slope drift
		+ `c2`: WG offset drift
		+ `Np`: Number of standardization parameters to fit
		+ `CM`: model covariance matrix for (`a`, `b`, `c`, `a2`, `b2`, `c2`)
		+ `d13Cwg_VPDB`: δ13C_VPDB of WG
		+ `d18Owg_VSMOW`: δ18O_VSMOW of WG
		'''

		for session in self.sessions:

			self.sessions[session] = {
				k: self.sessions[session][k]
				for k in ['data', 'N', 'scrambling_drift', 'slope_drift', 'wg_drift']
			}

			if 'd13Cwg_VPDB' not in self.sessions[session]:
				self.sessions[session]['d13Cwg_VPDB'] = self.sessions[session]['data'][0]['d13Cwg_VPDB']
			if 'd18Owg_VSMOW' not in self.sessions[session]:
				self.sessions[session]['d18Owg_VSMOW'] = self.sessions[session]['data'][0]['d18Owg_VSMOW']
			self.sessions[session]['Na'] = len([r for r in self.sessions[session]['data'] if r['Sample'] in self.anchors])
			self.sessions[session]['Nu'] = len([r for r in self.sessions[session]['data'] if r['Sample'] in self.unknowns])

			self.msg(f'Computing repeatabilities for session {session}')
			self.sessions[session]['r_d13C_VPDB'] = self.compute_r('d13C_VPDB', samples = 'anchors', sessions = [session])
			self.sessions[session]['r_d18O_VSMOW'] = self.compute_r('d18O_VSMOW', samples = 'anchors', sessions = [session])
			self.sessions[session][f'r_D{self._4x}'] = self.compute_r(f'D{self._4x}', sessions = [session])

		# target = self._resolve_target(target)
		# match target:
		# 	case 'bayes':
		# 		raise NotImplementedError
		# 	case 'pooled':
		# 		for session in self.sessions:
		# 			_s_ = self.standardization[target]['sessions'][session]

		# 			_s_['a'] = self.standardization[target]['lmfit'].params.valuesdict()[f'a_{pf(session)}']
		# 			i = self.standardization[target]['var_names'].index(f'a_{pf(session)}')
		# 			_s_['SE_a'] = self.standardization[target]['covar'][i,i]**.5

		# 			_s_['b'] = self.standardization[target]['lmfit'].params.valuesdict()[f'b_{pf(session)}']
		# 			i = self.standardization[target]['var_names'].index(f'b_{pf(session)}')
		# 			_s_['SE_b'] = self.standardization[target]['covar'][i,i]**.5

		# 			_s_['c'] = self.standardization[target]['lmfit'].params.valuesdict()[f'c_{pf(session)}']
		# 			i = self.standardization[target]['var_names'].index(f'c_{pf(session)}')
		# 			_s_['SE_c'] = self.standardization[target]['covar'][i,i]**.5

		# 			_s_['a2'] = self.standardization[target]['lmfit'].params.valuesdict()[f'a2_{pf(session)}']
		# 			if self.sessions[session]['scrambling_drift']:
		# 				i = self.standardization[target]['var_names'].index(f'a2_{pf(session)}')
		# 				_s_['SE_a2'] = self.standardization[target]['covar'][i,i]**.5
		# 			else:
		# 				_s_['SE_a2'] = 0.

		# 			_s_['b2'] = self.standardization[target]['lmfit'].params.valuesdict()[f'b2_{pf(session)}']
		# 			if self.sessions[session]['slope_drift']:
		# 				i = self.standardization[target]['var_names'].index(f'b2_{pf(session)}')
		# 				_s_['SE_b2'] = self.standardization[target]['covar'][i,i]**.5
		# 			else:
		# 				_s_['SE_b2'] = 0.

		# 			_s_['c2'] = self.standardization[target]['lmfit'].params.valuesdict()[f'c2_{pf(session)}']
		# 			if self.sessions[session]['wg_drift']:
		# 				i = self.standardization[target]['var_names'].index(f'c2_{pf(session)}')
		# 				_s_['SE_c2'] = self.standardization[target]['covar'][i,i]**.5
		# 			else:
		# 				_s_['SE_c2'] = 0.

		# 			i = self.standardization[target]['var_names'].index(f'a_{pf(session)}')
		# 			j = self.standardization[target]['var_names'].index(f'b_{pf(session)}')
		# 			k = self.standardization[target]['var_names'].index(f'c_{pf(session)}')
		# 			CM = np.zeros((6,6))
		# 			CM[:3,:3] = self.standardization[target]['covar'][[i,j,k],:][:,[i,j,k]]
		# 			try:
		# 				i2 = self.standardization[target]['var_names'].index(f'a2_{pf(session)}')
		# 				CM[3,[0,1,2,3]] = self.standardization[target]['covar'][i2,[i,j,k,i2]]
		# 				CM[[0,1,2,3],3] = self.standardization[target]['covar'][[i,j,k,i2],i2]
		# 				try:
		# 					j2 = self.standardization[target]['var_names'].index(f'b2_{pf(session)}')
		# 					CM[3,4] = self.standardization[target]['covar'][i2,j2]
		# 					CM[4,3] = self.standardization[target]['covar'][j2,i2]
		# 				except ValueError:
		# 					pass
		# 				try:
		# 					k2 = self.standardization[target]['var_names'].index(f'c2_{pf(session)}')
		# 					CM[3,5] = self.standardization[target]['covar'][i2,k2]
		# 					CM[5,3] = self.standardization[target]['covar'][k2,i2]
		# 				except ValueError:
		# 					pass
		# 			except ValueError:
		# 				pass
		# 			try:
		# 				j2 = self.standardization[target]['var_names'].index(f'b2_{pf(session)}')
		# 				CM[4,[0,1,2,4]] = self.standardization[target]['covar'][j2,[i,j,k,j2]]
		# 				CM[[0,1,2,4],4] = self.standardization[target]['covar'][[i,j,k,j2],j2]
		# 				try:
		# 					k2 = self.standardization[target]['var_names'].index(f'c2_{pf(session)}')
		# 					CM[4,5] = self.standardization[target]['covar'][j2,k2]
		# 					CM[5,4] = self.standardization[target]['covar'][k2,j2]
		# 				except ValueError:
		# 					pass
		# 			except ValueError:
		# 				pass
		# 			try:
		# 				k2 = self.standardization[target]['var_names'].index(f'c2_{pf(session)}')
		# 				CM[5,[0,1,2,5]] = self.standardization[target]['covar'][k2,[i,j,k,k2]]
		# 				CM[[0,1,2,5],5] = self.standardization[target]['covar'][[i,j,k,k2],k2]
		# 			except ValueError:
		# 				pass

		# 			_s_['CM'] = CM

	@make_verbal
	def repeatabilities(self, target = 'latest'):
		'''
		Compute analytical repeatabilities for δ13C_VPDB, δ18O_VSMOW, Δ4x
		(for all samples, for anchors, and for unknowns).
		'''
		self.msg('Computing reproducibilities for all sessions')

		self.repeatability['r_d13C_VPDB'] = self.compute_r('d13C_VPDB', samples = 'anchors')
		self.repeatability['r_d18O_VSMOW'] = self.compute_r('d18O_VSMOW', samples = 'anchors')
		self.repeatability[f'r_D{self._4x}a'] = self.compute_r(f'D{self._4x}', samples = 'anchors')
		self.repeatability[f'r_D{self._4x}u'] = self.compute_r(f'D{self._4x}', samples = 'unknowns')
		self.repeatability[f'r_D{self._4x}'] = self.compute_r(f'D{self._4x}', samples = 'all samples')


	@make_verbal
	def consolidate(
		self,
		target = 'latest',
		tables = True,
		plots = True,
	):
		'''
		Collect information about samples, sessions and repeatabilities.
		'''
		self.consolidate_samples(target = target)
		self.consolidate_sessions(target = target)
		self.repeatabilities(target = target)

		if tables:
			self.summary(target = target)
			self.table_of_sessions(target = target)
			self.table_of_analyses()
			self.table_of_samples(target = target)

		if plots:
			self.plot_sessions(target = target)


	@make_verbal
	def rmswd(self,
		samples = 'all samples',
		sessions = 'all sessions',
		):
		'''
		Compute the χ2, root mean squared weighted deviation
		(i.e. reduced χ2), and corresponding degrees of freedom of the
		Δ4x values for samples in `samples` and sessions in `sessions`.

		Only used in `D4xdata.standardize()` with `method='indep_sessions'`.
		'''
		if samples == 'all samples':
			mysamples = [k for k in self.samples]
		elif samples == 'anchors':
			mysamples = [k for k in self.anchors]
		elif samples == 'unknowns':
			mysamples = [k for k in self.unknowns]
		else:
			mysamples = samples

		if sessions == 'all sessions':
			sessions = [k for k in self.sessions]

		chisq, Nf = 0, 0
		for sample in mysamples :
			G = [ r for r in self if r['Sample'] == sample and r['Session'] in sessions ]
			if len(G) > 1 :
				X, sX = w_avg([r[f'D{self._4x}'] for r in G], [r[f'wD{self._4x}'] for r in G])
				Nf += (len(G) - 1)
				chisq += np.sum([ ((r[f'D{self._4x}']-X)/r[f'wD{self._4x}'])**2 for r in G])
		r = (chisq / Nf)**.5 if Nf > 0 else 0
		self.msg(f'RMSWD of r["D{self._4x}"] is {r:.6f} for {samples}.')
		return {'rmswd': r, 'chisq': chisq, 'Nf': Nf}


	@make_verbal
	def compute_r(self, key, samples = 'all samples', sessions = 'all sessions'):
		'''
		Compute the repeatability of `[r[key] for r in self]`
		'''

		# DESIGN PRINCIPLE:
		# Computing repeatabilities relies on the residuals of each analyses.
		# Thus it only applies to the latest standardization results, and
		# repeatabilities should always be recomputed at the end of standardize().

		if samples == 'all samples':
			mysamples = [k for k in self.samples]
		elif samples == 'anchors':
			mysamples = [k for k in self.anchors]
		elif samples == 'unknowns':
			mysamples = [k for k in self.unknowns]
		else:
			mysamples = samples

		if sessions == 'all sessions':
			sessions = [k for k in self.sessions]


		if key in ['D47', 'D48', 'D49']:

			stdz_method = self.standardization['latest']['method']
			stdz = self.standardization[stdz_method]

			match stdz_method:
				case 'bayes':
					raise NotImplementedError
				case 'pooled':
					# New computation in v3.0, based on leverage hat matrix (https://en.wikipedia.org/wiki/Leverage_(statistics))
					mask = [r['Sample'] in mysamples and r['Session'] in sessions for r in self] # booleans
					chi2 = np.sum(stdz['lmfit'].residual[mask]**2)
					dof = sum(mask) - stdz['lmfit'].h[mask].sum()
					r = np.sqrt(chi2 / dof) if abs(dof) > 1e-3 else 0.

		else: # if key not in ['D47', 'D48', 'D49']
			chisq, Nf = 0, 0
			for sample in mysamples :
				X = [ r[key] for r in self if r['Sample'] == sample and r['Session'] in sessions ]
				if len(X) > 1 :
					Nf += len(X) - 1
					chisq += np.sum([ (x-np.mean(X))**2 for x in X ])
			r = (chisq / Nf)**.5 if Nf > 0 else 0

		self.msg(f'Repeatability of r["{key}"] is {1000*r:.1f} ppm for {samples} in {sessions}.')
		return r

	def sample_average(self, samples, weights = 'equal', normalize = True):
		'''
		Weighted average Δ4x value of a group of samples, accounting for covariance.

		Returns the weighed average Δ4x value and associated SE
		of a group of samples. Weights are equal by default. If `normalize` is
		true, `weights` will be rescaled so that their sum equals 1.

		**Examples**

		```python
		self.sample_average(['X','Y'], [1, 2])
		```

		returns the value and SE of [Δ4x(X) + 2 Δ4x(Y)]/3,
		where Δ4x(X) and Δ4x(Y) are the average Δ4x
		values of samples X and Y, respectively.

		```python
		self.sample_average(['X','Y'], [1, -1], normalize = False)
		```

		returns the value and SE of the difference Δ4x(X) - Δ4x(Y).
		'''
		if weights == 'equal':
			weights = [1/len(samples)] * len(samples)

		if normalize:
			s = sum(weights)
			if s:
				weights = [w/s for w in weights]

		try:
# 			indices = [self.standardization.var_names.index(f'D47_{pf(sample)}') for sample in samples]
# 			C = self.standardization.covar[indices,:][:,indices]
			C = np.array([[self.sample_D4x_covar(x, y) for x in samples] for y in samples])
			X = [self.samples[sample][f'D{self._4x}'] for sample in samples]
			return correlated_sum(X, C, weights)
		except ValueError:
			return (0., 0.)


	def sample_D4x_covar(self, sample1, sample2 = None, target = 'latest'):
		'''
		Covariance between Δ4x values of samples

		Returns the error covariance between the average Δ4x values of two
		samples. If if only `sample_1` is specified, or if `sample_1 == sample_2`),
		returns the Δ4x variance for that sample.
		'''

		target = self._resolve_target(target)

		if sample2 is None:
			sample2 = sample1
		match target:
			case 'pooled':
				i = self.standardization[target]['var_names'].index(f'D{self._4x}_{pf(sample1)}')
				j = self.standardization[target]['var_names'].index(f'D{self._4x}_{pf(sample2)}')
				return self.standardization[target]['covar'][i, j]
			case 'bayes':
				raise NotImplementedError


	def sample_D4x_correl(self, sample1, sample2 = None):
		'''
		Correlation between Δ4x errors of samples

		Returns the error correlation between the average Δ4x values of two samples.
		'''
		if sample2 is None or sample2 == sample1:
			return 1.
		return (
			self.sample_D4x_covar(sample1, sample2)
			/ self.unknowns[sample1][f'SE_D{self._4x}']
			/ self.unknowns[sample2][f'SE_D{self._4x}']
			)

	def plot_single_pooled_session(self,
		session,
		kw_plot_anchors = dict(ls='None', marker='x', mec=(.75, 0, 0), mew = .75, ms = 4),
		kw_plot_unknowns = dict(ls='None', marker='x', mec=(0, 0, .75), mew = .75, ms = 4),
		kw_plot_anchor_avg = dict(ls='-', marker='None', color=(.75, 0, 0), lw = .75),
		kw_plot_unknown_avg = dict(ls='-', marker='None', color=(0, 0, .75), lw = .75),
		kw_contour_error = dict(colors = [[0, 0, 0]], alpha = .5, linewidths = 0.75),
		xylimits = 'free', # | 'constant'
		x_label = None,
		y_label = None,
		error_contour_interval = 'auto',
		fig = 'new',
		):
		'''
		Generate plot for a single session
		'''
		if x_label is None:
			x_label = f'δ$_{{{self._4x}}}$ (‰)'
		if y_label is None:
			y_label = f'Δ$_{{{self._4x}}}$ (‰)'

		out = _SessionPlot()
		anchors = [a for a in self.anchors if [r for r in self.sessions[session]['data'] if r['Sample'] == a]]
		unknowns = [u for u in self.unknowns if [r for r in self.sessions[session]['data'] if r['Sample'] == u]]
		anchors_d = [r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in self.anchors]
		anchors_D = [r[f'D{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in self.anchors]
		unknowns_d = [r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in self.unknowns]
		unknowns_D = [r[f'D{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in self.unknowns]
		anchor_avg = (
			np.array([
				np.array([
					np.min([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) - 1,
					np.max([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) + 1,
				])
				for sample in anchors
			]).T,
			np.array([
				np.array([0, 0]) + self.bare_RMs[sample]
				for sample in anchors
			]).T,
		)
		unknown_avg = (
			np.array([
				np.array([
					np.min([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) - 1,
					np.max([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) + 1,
				])
				for sample in unknowns
			]).T,
			np.array([
				np.array([0, 0]) + self.unknowns[sample][f'D{self._4x}']
				for sample in unknowns
			]).T,
		)

		if fig == 'new':
			out.fig = ppl.figure(figsize = (6,6))
			ppl.subplots_adjust(.1,.1,.9,.9)

		out.anchor_analyses, = ppl.plot(
			anchors_d,
			anchors_D,
			**kw_plot_anchors)
		out.unknown_analyses, = ppl.plot(
			unknowns_d,
			unknowns_D,
			**kw_plot_unknowns)
		out.anchor_avg = ppl.plot(
			*anchor_avg,
			**kw_plot_anchor_avg)
		out.unknown_avg = ppl.plot(
			*unknown_avg,
			**kw_plot_unknown_avg)
		if xylimits == 'constant':
			x = [r[f'd{self._4x}'] for r in self]
			y = [r[f'D{self._4x}'] for r in self]
			x1, x2, y1, y2 = np.min(x), np.max(x), np.min(y), np.max(y)
			w, h = x2-x1, y2-y1
			x1 -= w/20
			x2 += w/20
			y1 -= h/20
			y2 += h/20
			ppl.axis([x1, x2, y1, y2])
		elif xylimits == 'free':
			x1, x2, y1, y2 = ppl.axis()
		else:
			x1, x2, y1, y2 = ppl.axis(xylimits)

		if error_contour_interval != 'none':
			xi, yi = np.linspace(x1, x2), np.linspace(y1, y2)
			XI,YI = np.meshgrid(xi, yi)
			SI = np.array([[self.standardization_error(session, x, y, target = 'pooled') for x in xi] for y in yi])
			if error_contour_interval == 'auto':
				rng = np.max(SI) - np.min(SI)
				if rng <= 0.01:
					cinterval = 0.001
				elif rng <= 0.03:
					cinterval = 0.004
				elif rng <= 0.1:
					cinterval = 0.01
				elif rng <= 0.3:
					cinterval = 0.03
				elif rng <= 1.:
					cinterval = 0.1
				else:
					cinterval = 0.5
			else:
				cinterval = error_contour_interval

			cval = np.arange(np.ceil(SI.min() / .001) * .001, np.ceil(SI.max() / .001 + 1) * .001, cinterval)
			out.contour = ppl.contour(XI, YI, SI, cval, **kw_contour_error)
			out.clabel = ppl.clabel(out.contour)
			contour = (XI, YI, SI, cval, cinterval)

		if fig == None:
			return {
			'anchors':anchors,
			'unknowns':unknowns,
			'anchors_d':anchors_d,
			'anchors_D':anchors_D,
			'unknowns_d':unknowns_d,
			'unknowns_D':unknowns_D,
			'anchor_avg':anchor_avg,
			'unknown_avg':unknown_avg,
			'contour':contour,
			}

		ppl.xlabel(x_label)
		ppl.ylabel(y_label)
		ppl.title(session, weight = 'bold')
		ppl.grid(alpha = .2)
		out.ax = ppl.gca()

		return out

	def plot_residuals(
		self,
		kde = False,
		hist = False,
		binwidth = 2/3,
		dir = 'output',
		filename = None,
		highlight = [],
		colors = None,
		figsize = None,
		dpi = 100,
		yspan = None,
		savefig = True,
		):
		'''
		Plot residuals of each analysis as a function of time (actually, as a function of
		the order of analyses in the `D4xdata` object)

		+ `kde`: whether to add a kernel density estimate of residuals
		+ `hist`: whether to add a histogram of residuals (incompatible with `kde`)
		+ `histbins`: specify bin edges for the histogram
		+ `dir`: the directory in which to save the plot
		+ `highlight`: a list of samples to highlight
		+ `colors`: a dict of `{<sample>: (r, g, b)}` for all samples
		+ `figsize`: (width, height) of figure
		+ `dpi`: resolution for PNG output
		+ `yspan`: factor controlling the range of y values shown in plot
		  (by default: `yspan = 1.5 if kde else 1.0`)
		+ `savefig`: whether to export the figure to a file (the figure is then closed);
		  if `False`, return the `Figure()` instance instead of closing it.
		'''

		from matplotlib import ticker

		half_span_95CL = self.repeatability[f'r_D{self._4x}']*1000*self.standardization['latest']['t95']

		if yspan is None:
			if kde:
				yspan = 1.5
			else:
				yspan = 1.0

		# Layout
		fig = ppl.figure(figsize = (8,4) if figsize is None else figsize)
		if hist or kde:
			ppl.subplots_adjust(left = .08, bottom = .05, right = .98, top = .8, wspace = -0.72)
			ax1, ax2 = ppl.subplot(121), ppl.subplot(1,15,15)
		else:
			ppl.subplots_adjust(.08,.05,.78,.8)
			ax1 = ppl.subplot(111)

		# Colors
		N = len(self.anchors)
		if colors is None:
			if len(highlight) > 0:
				Nh = len(highlight)
				if Nh == 1:
					colors = {highlight[0]: (0,0,0)}
				elif Nh == 3:
					colors = {a: c for a,c in zip(highlight, [(0,0,1), (1,0,0), (0,2/3,0)])}
				elif Nh == 4:
					colors = {a: c for a,c in zip(highlight, [(0,0,1), (1,0,0), (0,2/3,0), (.75,0,.75)])}
				else:
					colors = {a: hls_to_rgb(k/Nh, .4, 1) for k,a in enumerate(highlight)}
			else:
				if N == 3:
					colors = {a: c for a,c in zip(self.anchors, [(0,0,1), (1,0,0), (0,2/3,0)])}
				elif N == 4:
					colors = {a: c for a,c in zip(self.anchors, [(0,0,1), (1,0,0), (0,2/3,0), (.75,0,.75)])}
				else:
					colors = {a: hls_to_rgb(k/N, .4, 1) for k,a in enumerate(self.anchors)}

		ppl.sca(ax1)

		ppl.axhline(0, color = 'k', alpha = .25, lw = 0.75)

		ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'${x:+.0f}$' if x else '$0$'))

		session = self[0]['Session']
		x1 = 0
# 		ymax = np.max([1e3 * (r['D47'] - self.samples[r['Sample']]['D47']) for r in self])
		x_sessions = {}
		one_or_more_singlets = False
		one_or_more_multiplets = False
		multiplets = set()
		for k,r in enumerate(self):
			if r['Session'] != session:
				x2 = k-1
				x_sessions[session] = (x1+x2)/2
				ppl.axvline(k - 0.5, color = 'k', lw = .5)
				session = r['Session']
				x1 = k
			singlet = len(self.samples[r['Sample']]['data']) == 1
			if not singlet:
				multiplets.add(r['Sample'])
			if r['Sample'] in self.unknowns:
				if singlet:
					one_or_more_singlets = True
				else:
					one_or_more_multiplets = True
			kw = dict(
				marker = 'x' if singlet else '+',
				ms = 4 if singlet else 5,
				ls = 'None',
				mec = colors[r['Sample']] if r['Sample'] in colors else (0,0,0),
				mew = 1,
				alpha = 0.2 if singlet else 1,
				)
			if highlight and r['Sample'] not in highlight:
				kw['alpha'] = 0.2
			ppl.plot(k, 1e3 * r[f'D{self._4x}_residual'], **kw)
		x2 = k
		x_sessions[session] = (x1+x2)/2

		ppl.axhspan(-self.repeatability[f'r_D{self._4x}']*1000, self.repeatability[f'r_D{self._4x}']*1000, color = 'k', alpha = .05, lw = 1)
		ppl.axhspan(-half_span_95CL, half_span_95CL, color = 'k', alpha = .05, lw = 1)
		if not (hist or kde):
			ppl.text(len(self), self.repeatability[f'r_D{self._4x}']*1000, f"   SD = {self.repeatability[f'r_D{self._4x}']*1000:.1f} ppm", size = 9, alpha = 1, va = 'center')
			ppl.text(len(self), half_span_95CL, f"   95% CL = ± {half_span_95CL:.1f} ppm", size = 9, alpha = 1, va = 'center')

		xmin, xmax, ymin, ymax = ppl.axis()
		if yspan != 1:
			ymin, ymax = (ymin + ymax)/2 - yspan * (ymax - ymin)/2, (ymin + ymax)/2 + yspan * (ymax - ymin)/2
		for s in x_sessions:
			ppl.text(
				x_sessions[s],
				ymax +1,
				s,
				va = 'bottom',
				**(
					dict(ha = 'center')
					if len(self.sessions[s]['data']) > (0.15 * len(self))
					else dict(ha = 'left', rotation = 45)
					)
				)

		if hist or kde:
			ppl.sca(ax2)

		for s in colors:
			kw['marker'] = '+'
			kw['ms'] = 5
			kw['mec'] = colors[s]
			kw['label'] = s
			kw['alpha'] = 1
			ppl.plot([], [], **kw)

		kw['mec'] = (0,0,0)

		if one_or_more_singlets:
			kw['marker'] = 'x'
			kw['ms'] = 4
			kw['alpha'] = .2
			kw['label'] = 'other (N$\\,$=$\\,$1)' if one_or_more_multiplets else 'other'
			ppl.plot([], [], **kw)

		if one_or_more_multiplets:
			kw['marker'] = '+'
			kw['ms'] = 4
			kw['alpha'] = 1
			kw['label'] = 'other (N$\\,$>$\\,$1)' if one_or_more_singlets else 'other'
			ppl.plot([], [], **kw)

		if hist or kde:
			leg = ppl.legend(loc = 'upper right', bbox_to_anchor = (1, 1), bbox_transform=fig.transFigure, borderaxespad = 1.5, fontsize = 9)
		else:
			leg = ppl.legend(loc = 'lower right', bbox_to_anchor = (1, 0), bbox_transform=fig.transFigure, borderaxespad = 1.5)
		leg.set_zorder(-1000)

		ppl.sca(ax1)

		ppl.ylabel(f'Δ$_{{{self._4x}}}$ residuals (ppm)')
		ppl.xticks([])
		ppl.axis([-1, len(self), None, None])

		if hist or kde:
			ppl.sca(ax2)
			X = 1e3 * np.array([r[f'D{self._4x}_residual'] for r in self if r['Sample'] in multiplets or r['Sample'] in self.anchors])

			if kde:
				from scipy.stats import gaussian_kde
				yi = np.linspace(ymin, ymax, 201)
				xi = gaussian_kde(X).evaluate(yi)
				ppl.fill_betweenx(yi, xi, xi*0, fc = (0,0,0,.15), lw = 1, ec = (.75,.75,.75,1))
# 				ppl.plot(xi, yi, 'k-', lw = 1)
			elif hist:
				ppl.hist(
					X,
					orientation = 'horizontal',
					histtype = 'stepfilled',
					ec = [.4]*3,
					fc = [.25]*3,
					alpha = .25,
					bins = np.linspace(-9e3*self.repeatability[f'r_D{self._4x}'], 9e3*self.repeatability[f'r_D{self._4x}'], int(18/binwidth+1)),
					)
			ppl.text(0, 0,
				f"   SD = {self.repeatability[f'r_D{self._4x}']*1000:.1f} ppm\n   95% CL = ± {half_span_95CL:.1f} ppm",
				size = 7.5,
				alpha = 1,
				va = 'center',
				ha = 'left',
				)

			ppl.axis([0, None, ymin, ymax])
			ppl.xticks([])
			ppl.yticks([])
# 			ax2.spines['left'].set_visible(False)
			ax2.spines['right'].set_visible(False)
			ax2.spines['top'].set_visible(False)
			ax2.spines['bottom'].set_visible(False)

		ax1.axis([None, None, ymin, ymax])


		if savefig:
			if not os.path.exists(dir):
				os.makedirs(dir)
			if filename is None:
				filename = f'D{self._4x}_residuals.pdf'
			ppl.savefig(f'{dir}/{filename}', dpi = dpi)
			print('done')
			ppl.close(fig)
		else:
			return fig


	def simulate(self, *args, **kwargs):
		'''
		Legacy function with warning message pointing to `virtual_data()`
		'''
		raise DeprecationWarning('D4xdata.simulate is deprecated and has been replaced by virtual_data()')

	def plot_anchor_residuals(
		self,
		dir = 'output',
		filename = '',
		figsize = None,
		subplots_adjust = (0.05, 0.1, 0.95, 0.98, .25, .25),
		dpi = 100,
		colors = None,
		):
		'''
		Plot a summary of the residuals for all anchors, intended to help detect systematic bias.

		**Parameters**

		+ `dir`: the directory in which to save the plot
		+ `filename`: the file name to save to.
		+ `dpi`: resolution for PNG output
		+ `figsize`: (width, height) of figure
		+ `subplots_adjust`: passed to the figure
		+ `dpi`: resolution for PNG output
		+ `colors`: a dict of `{<sample>: (r, g, b)}` for all samples
		'''

		# Colors
		N = len(self.anchors)
		if colors is None:
			if N == 3:
				colors = {a: c for a,c in zip(self.anchors, [(0,0,1), (1,0,0), (0,2/3,0)])}
			elif N == 4:
				colors = {a: c for a,c in zip(self.anchors, [(0,0,1), (1,0,0), (0,2/3,0), (.75,0,.75)])}
			else:
				colors = {a: hls_to_rgb(k/N, .4, 1) for k,a in enumerate(self.anchors)}

		if figsize is None:
			figsize = (4, 1.5*N+1)
		fig = ppl.figure(figsize = figsize)
		ppl.subplots_adjust(*subplots_adjust)
		axs = {}
		X = np.array([r[f'D{self._4x}_residual'] for a in self.anchors for r in self.anchors[a]['data']])*1000
		sigma = self.repeatability['r_D47a'] * 1000
		D = max(np.abs(X))

		for k,a in enumerate(self.anchors):
			color = colors[a]
			axs[a] = ppl.subplot(N, 1, 1+k)
			axs[a].text(
				0.02, 1-0.05, a,
				va = 'top',
				ha = 'left',
				weight = 'bold',
				size = 9,
				color = [_*0.75 for _ in color],
				transform = axs[a].transAxes,
			)
			X = np.array([r[f'D{self._4x}_residual'] for r in self.anchors[a]['data']])*1000
			axs[a].axvline(0, lw = 0.5, color = color)
			axs[a].plot(X, X*0, 'o', mew = 0.7, mec = (*color,.5), mfc = (*color, 0), ms = 7, clip_on = False)

			xi = np.linspace(-3*D, 3*D, 601)
			yi = np.array([np.exp(-0.5 * ((xi - x)/sigma)**2) for x in X]).sum(0)
			ppl.fill_between(xi, yi, yi*0, fc = (*color, .15), lw = 1, ec = color)

			axs[a].errorbar(
				X.mean(), yi.max()*.2, None, 1.96*sigma/len(X)**0.5,
				ecolor = color,
				marker = 's',
				ls = 'None',
				mec = color,
				mew = 1,
				mfc = 'w',
				ms = 8,
				elinewidth = 1,
				capsize = 4,
				capthick = 1,
			)

			axs[a].axis([xi[0], xi[-1], 0, yi.max()*1.05])
			ppl.yticks([])

		ppl.xlabel(f'$Δ_{{{self._4x}}}$ residuals (ppm)')

		if not os.path.exists(dir):
			os.makedirs(dir)
		if filename is None:
			return fig
		elif filename == '':
			filename = f'D{self._4x}_anchor_residuals.pdf'
		ppl.savefig(f'{dir}/{filename}', dpi = dpi)
		ppl.close(fig)


	def plot_distribution_of_analyses(
		self,
		dir = 'output',
		filename = None,
		vs_time = False,
		figsize = (6,4),
		subplots_adjust = (0.02, 0.13, 0.85, 0.8),
		output = None,
		dpi = 100,
		):
		'''
		Plot temporal distribution of all analyses in the data set.

		**Parameters**

		+ `dir`: the directory in which to save the plot
		+ `vs_time`: if `True`, plot as a function of `TimeTag` rather than sequentially.
		+ `dpi`: resolution for PNG output
		+ `figsize`: (width, height) of figure
		+ `dpi`: resolution for PNG output
		'''

		asamples = [s for s in self.anchors]
		usamples = [s for s in self.unknowns]
		if output is None or output == 'fig':
			fig = ppl.figure(figsize = figsize)
			ppl.subplots_adjust(*subplots_adjust)
		Xmin = min([r['TimeTag'] if vs_time else j for j,r in enumerate(self)])
		Xmax = max([r['TimeTag'] if vs_time else j for j,r in enumerate(self)])
		Xmax += (Xmax-Xmin)/40
		Xmin -= (Xmax-Xmin)/41
		for k, s in enumerate(asamples + usamples):
			if vs_time:
				X = [r['TimeTag'] for r in self if r['Sample'] == s]
			else:
				X = [x for x,r in enumerate(self) if r['Sample'] == s]
			Y = [-k for x in X]
			ppl.plot(X, Y, 'o', mec = None, mew = 0, mfc = 'b' if s in usamples else 'r', ms = 3, alpha = .75)
			ppl.axhline(-k, color = 'b' if s in usamples else 'r', lw = .5, alpha = .25)
			ppl.text(Xmax, -k, f'   {s}', va = 'center', ha = 'left', size = 7, color = 'b' if s in usamples else 'r')
		ppl.axis([Xmin, Xmax, -k-1, 1])
		ppl.xlabel('\ntime')
		ppl.gca().annotate('',
			xy = (0.6, -0.02),
			xycoords = 'axes fraction',
			xytext = (.4, -0.02),
			arrowprops = dict(arrowstyle = "->", color = 'k'),
		)


		x2 = -1
		for session in self.sessions:
			x1 = min([r['TimeTag'] if vs_time else j for j,r in enumerate(self) if r['Session'] == session])
			if vs_time:
				ppl.axvline(x1, color = 'k', lw = .75)
			if x2 > -1:
				if not vs_time:
					ppl.axvline((x1+x2)/2, color = 'k', lw = .75, alpha = .5)
			x2 = max([r['TimeTag'] if vs_time else j for j,r in enumerate(self) if r['Session'] == session])
# 			from xlrd import xldate_as_datetime
# 			print(session, xldate_as_datetime(x1, 0), xldate_as_datetime(x2, 0))
			if vs_time:
				ppl.axvline(x2, color = 'k', lw = .75)
				ppl.axvspan(x1,x2,color = 'k', zorder = -100, alpha = .15)
			ppl.text((x1+x2)/2, 1, f' {session}', ha = 'left', va = 'bottom', rotation = 45, size = 8)

		ppl.xticks([])
		ppl.yticks([])

		if output is None:
			if not os.path.exists(dir):
				os.makedirs(dir)
			if filename == None:
				filename = f'D{self._4x}_distribution_of_analyses.pdf'
			ppl.savefig(f'{dir}/{filename}', dpi = dpi)
			ppl.close(fig)
		elif output == 'ax':
			return ppl.gca()
		elif output == 'fig':
			return fig


	def plot_bulk_compositions(
		self,
		samples = None,
		dir = 'output/bulk_compositions',
		figsize = (6,6),
		subplots_adjust = (0.15, 0.12, 0.95, 0.92),
		show = False,
		sample_color = (0,.5,1),
		analysis_color = (.7,.7,.7),
		labeldist = 0.3,
		radius = 0.05,
		):
		'''
		Plot δ13C_VBDP vs δ18O_VSMOW (of CO2) for all analyses.

		By default, creates a directory `./output/bulk_compositions` where plots for
		each sample are saved. Another plot named `__all__.pdf` shows all analyses together.


		**Parameters**

		+ `samples`: Only these samples are processed (by default: all samples).
		+ `dir`: where to save the plots
		+ `figsize`: (width, height) of figure
		+ `subplots_adjust`: passed to `subplots_adjust()`
		+ `show`: whether to call `matplotlib.pyplot.show()` on the plot with all samples,
		allowing for interactive visualization/exploration in (δ13C, δ18O) space.
		+ `sample_color`: color used for replicate markers/labels
		+ `analysis_color`: color used for sample markers/labels
		+ `labeldist`: distance (in inches) from replicate markers to replicate labels
		+ `radius`: radius of the dashed circle providing scale. No circle if `radius = 0`.
		'''

		from matplotlib.patches import Ellipse

		if samples is None:
			samples = [_ for _ in self.samples]

		saved = {}

		for s in samples:

			fig = ppl.figure(figsize = figsize)
			fig.subplots_adjust(*subplots_adjust)
			ax = ppl.subplot(111)
			ppl.xlabel('$δ^{18}O_{VSMOW}$ of $CO_2$ (‰)')
			ppl.ylabel('$δ^{13}C_{VPDB}$ (‰)')
			ppl.title(s)


			XY = np.array([[_['d18O_VSMOW'], _['d13C_VPDB']] for _ in self.samples[s]['data']])
			UID = [_['UID'] for _ in self.samples[s]['data']]
			XY0 = XY.mean(0)

			for xy in XY:
				ppl.plot([xy[0], XY0[0]], [xy[1], XY0[1]], '-', lw = 1, color = analysis_color)

			ppl.plot(*XY.T, 'wo', mew = 1, mec = analysis_color)
			ppl.plot(*XY0, 'wo', mew = 2, mec = sample_color)
			ppl.text(*XY0, f'  {s}', va = 'center', ha = 'left', color = sample_color, weight = 'bold')
			saved[s] = [XY, XY0]

			x1, x2, y1, y2 = ppl.axis()
			x0, dx = (x1+x2)/2, (x2-x1)/2
			y0, dy = (y1+y2)/2, (y2-y1)/2
			dx, dy = [max(max(dx, dy), radius)]*2

			ppl.axis([
				x0 - 1.2*dx,
				x0 + 1.2*dx,
				y0 - 1.2*dy,
				y0 + 1.2*dy,
				])

			XY0_in_display_space = fig.dpi_scale_trans.inverted().transform(ax.transData.transform(XY0))

			for xy, uid in zip(XY, UID):

				xy_in_display_space = fig.dpi_scale_trans.inverted().transform(ax.transData.transform(xy))
				vector_in_display_space = xy_in_display_space - XY0_in_display_space

				if (vector_in_display_space**2).sum() > 0:

					unit_vector_in_display_space = vector_in_display_space / ((vector_in_display_space**2).sum())**0.5
					label_vector_in_display_space = vector_in_display_space + unit_vector_in_display_space * labeldist
					label_xy_in_display_space = XY0_in_display_space + label_vector_in_display_space
					label_xy_in_data_space = ax.transData.inverted().transform(fig.dpi_scale_trans.transform(label_xy_in_display_space))

					ppl.text(*label_xy_in_data_space, uid, va = 'center', ha = 'center', color = analysis_color)

				else:

					ppl.text(*xy, f'{uid}  ', va = 'center', ha = 'right', color = analysis_color)

			if radius:
				ax.add_artist(Ellipse(
					xy = XY0,
					width = radius*2,
					height = radius*2,
					ls = (0, (2,2)),
					lw = .7,
					ec = analysis_color,
					fc = 'None',
					))
				ppl.text(
					XY0[0],
					XY0[1]-radius,
					f'\n± {radius*1e3:.0f} ppm',
					color = analysis_color,
					va = 'top',
					ha = 'center',
					linespacing = 0.4,
					size = 8,
					)

			if not os.path.exists(dir):
				os.makedirs(dir)
			fig.savefig(f'{dir}/{s}.pdf')
			ppl.close(fig)

		fig = ppl.figure(figsize = figsize)
		fig.subplots_adjust(*subplots_adjust)
		ppl.xlabel('$δ^{18}O_{VSMOW}$ of $CO_2$ (‰)')
		ppl.ylabel('$δ^{13}C_{VPDB}$ (‰)')

		for s in saved:
			for xy in saved[s][0]:
				ppl.plot([xy[0], saved[s][1][0]], [xy[1], saved[s][1][1]], '-', lw = 1, color = analysis_color)
			ppl.plot(*saved[s][0].T, 'wo', mew = 1, mec = analysis_color)
			ppl.plot(*saved[s][1], 'wo', mew = 1.5, mec = sample_color)
			ppl.text(*saved[s][1], f'  {s}', va = 'center', ha = 'left', color = sample_color, weight = 'bold')

		x1, x2, y1, y2 = ppl.axis()
		ppl.axis([
			x1 - (x2-x1)/10,
			x2 + (x2-x1)/10,
			y1 - (y2-y1)/10,
			y2 + (y2-y1)/10,
			])


		if not os.path.exists(dir):
			os.makedirs(dir)
		fig.savefig(f'{dir}/__all__.pdf')
		if show:
			ppl.show()
		ppl.close(fig)


	def _save_D4x_correl(
		self,
		samples = None,
		dir = 'output',
		filename = None,
		D4x_precision = 4,
		correl_precision = 4,
		save_to_file = True,
		):
		'''
		Save D4x values along with their SE and correlation matrix.

		**Parameters**

		+ `samples`: Only these samples are output (by default: all samples).
		+ `dir`: the directory in which to save the faile (by defaut: `output`)
		+ `filename`: the name to the csv file to write to (by default: `D4x_correl.csv`)
		+ `D4x_precision`: the precision to use when writing `D4x` and `D4x_SE` values (by default: 4)
		+ `correl_precision`: the precision to use when writing correlation factor values (by default: 4)
		+ `save_to_file`: whether to write the output to a file factor values (by default: True). If `False`,
		returns the output as a string
		'''
		if samples is None:
			samples = sorted([s for s in self.unknowns])

		out = [['Sample']] + [[s] for s in samples]
		out[0] += [f'D{self._4x}', f'D{self._4x}_SE', f'D{self._4x}_correl']
		for k,s in enumerate(samples):
			out[k+1] += [f'{self.samples[s][f"D{self._4x}"]:.4f}', f'{self.samples[s][f"SE_D{self._4x}"]:.4f}']
			for s2 in samples:
				out[k+1] += [f'{self.sample_D4x_correl(s,s2):.4f}']

		if save_to_file:
			if not os.path.exists(dir):
				os.makedirs(dir)
			if filename is None:
				filename = f'D{self._4x}_correl.csv'
			with open(f'{dir}/{filename}', 'w') as fid:
				fid.write(make_csv(out))
		else:
			return make_csv(out)

	def pprint(self):
		from rich.pretty import pprint as pp
		pp(_pp({k:v for k,v in self.__dict__.items()}))

def _pp(x):
	if isinstance(x, np.float64):
		return float(x)
	if isinstance(x, np.ndarray):
		if len(x.shape) > 1:
			with np.printoptions(formatter={'float_kind': lambda x: f"{x: 7.2e}"}):
				return {k:str(x[k]) for k in range(x.shape[0])}
		elif x.size < 10:
			return x
		else:
			return '<large array>'
	if isinstance(x, dict):
		return {
			k: '<large list>' if k == 'data' else _pp(v)
			for k,v in x.items()
		}
	return x


class D47data(D4xdata):
	'''
	Store and process data for a large set of Δ47 analyses,
	usually comprising more than one analytical session.
	'''

	RMs = {
		'ETH-1':   0.2052,
		'ETH-2':   0.2085,
		'ETH-3':   0.6132,
		'ETH-4':   (0.4511, 0.0011),
		'IAEA-C1': (0.3018, 0.0011),
		'IAEA-C2': (0.6409, 0.0018),
		'MERCK':   (0.5135, 0.0024),
		} # I-CDES (Bernasconi et al., 2021)
	'''
	Nominal Δ47 values assigned to the Δ47 anchor samples, used by
	`D47data.standardize()` to normalize unknown samples to an absolute Δ47
	reference frame.

	By default equal to (after [Bernasconi et al. (2021)](https://doi.org/10.1029/2020GC009588)):
	```py
	{
		'ETH-1'   : 0.2052,
		'ETH-2'   : 0.2085,
		'ETH-3'   : 0.6132,
		'ETH-4'   : 0.4511,
		'IAEA-C1' : 0.3018,
		'IAEA-C2' : 0.6409,
		'MERCK'   : 0.5135,
	}
	```
	'''


	# @property
	# def Nominal_D47(self):
	# 	return self.RMs


	# @Nominal_D47.setter
	# def Nominal_D47(self, new):
	# 	self.RMs = dict(**new)
	# 	self.refresh()


	def __init__(self, l = [], **kwargs):
		'''
		**Parameters:** same as `D4xdata.__init__()`
		'''
		D4xdata.__init__(self, l = l, mass = '47', **kwargs)


	def D47fromTeq(self, fCo2eqD47 = 'petersen', priority = 'new'):
		'''
		Find all samples for which `Teq` is specified, compute equilibrium Δ47
		value for that temperature, and add treat these samples as additional anchors.

		**Parameters**

		+ `fCo2eqD47`: Which CO2 equilibrium law to use
		(`petersen`: [Petersen et al. (2019)](https://doi.org/10.1029/2018GC008127);
		`wang`: [Wang et al. (2019)](https://doi.org/10.1016/j.gca.2004.05.039)).
		+ `priority`: if `replace`: forget old anchors and only use the new ones;
		if `new`: keep pre-existing anchors but update them in case of conflict
		between old and new Δ47 values;
		if `old`: keep pre-existing anchors but preserve their original Δ47
		values in case of conflict.
		'''
		f = {
			'petersen': fCO2eqD47_Petersen,
			'wang': fCO2eqD47_Wang,
			}[fCo2eqD47]
		foo = {}
		for r in self:
			if 'Teq' in r:
				if r['Sample'] in foo:
					assert foo[r['Sample']] == f(r['Teq']), f'Different values of `Teq` provided for sample `{r["Sample"]}`.'
				else:
					foo[r['Sample']] = f(r['Teq'])
			else:
					assert r['Sample'] not in foo, f'`Teq` is inconsistently specified for sample `{r["Sample"]}`.'

		if priority == 'replace':
			self.RMs = {}
		for s in foo:
			if priority != 'old' or s not in self.RMs:
				self.RMs[s] = foo[s]

	def save_D47_correl(self, *args, **kwargs):
		return self._save_D4x_correl(*args, **kwargs)

	save_D47_correl.__doc__ = D4xdata._save_D4x_correl.__doc__.replace('D4x', 'D47')


class D48data(D4xdata):
	'''
	Store and process data for a large set of Δ48 analyses,
	usually comprising more than one analytical session.
	'''

	RMs = {
		'ETH-1':  0.138,
		'ETH-2':  0.138,
		'ETH-3':  0.270,
		'ETH-4':  0.223,
		'GU-1':  -0.419,
		} # (Fiebig et al., 2019, 2021)
	'''
	Nominal Δ48 values assigned to the Δ48 anchor samples, used by
	`D48data.standardize()` to normalize unknown samples to an absolute Δ48
	reference frame.

	By default equal to (after [Fiebig et al. (2019)](https://doi.org/10.1016/j.chemgeo.2019.05.019),
	[Fiebig et al. (2021)](https://doi.org/10.1016/j.gca.2021.07.012)):

	```py
	{
		'ETH-1' :  0.138,
		'ETH-2' :  0.138,
		'ETH-3' :  0.270,
		'ETH-4' :  0.223,
		'GU-1'  : -0.419,
	}
	```
	'''

	def __init__(self, l = [], **kwargs):
		'''
		**Parameters:** same as `D4xdata.__init__()`
		'''
		D4xdata.__init__(self, l = l, mass = '48', **kwargs)

	def save_D48_correl(self, *args, **kwargs):
		return self._save_D4x_correl(*args, **kwargs)

	save_D48_correl.__doc__ = D4xdata._save_D4x_correl.__doc__.replace('D4x', 'D48')


class D49data(D4xdata):
	'''
	Store and process data for a large set of Δ49 analyses,
	usually comprising more than one analytical session.
	'''

	RMs = {"1000C": 0.0, "25C": 2.228}  # Wang 2004
	'''
	Nominal Δ49 values assigned to the Δ49 anchor samples, used by
	`D49data.standardize()` to normalize unknown samples to an absolute Δ49
	reference frame.

	By default equal to (after [Wang et al. (2004)](https://doi.org/10.1016/j.gca.2004.05.039)):

	```py
	{
		"1000C": 0.0,
		"25C": 2.228
	}
	```
	'''

	# @property
	# def Nominal_D49(self):
	# 	return self.RMs

	# @Nominal_D49.setter
	# def Nominal_D49(self, new):
	# 	self.RMs = dict(**new)
	# 	self.refresh()

	def __init__(self, l=[], **kwargs):
		'''
		**Parameters:** same as `D4xdata.__init__()`
		'''
		D4xdata.__init__(self, l=l, mass='49', **kwargs)

	def save_D49_correl(self, *args, **kwargs):
		return self._save_D4x_correl(*args, **kwargs)

	save_D49_correl.__doc__ = D4xdata._save_D4x_correl.__doc__.replace('D4x', 'D49')

class _SessionPlot():
	'''
	Simple placeholder class
	'''
	def __init__(self):
		pass

_app = typer.Typer(
	add_completion = False,
	context_settings={'help_option_names': ['-h', '--help']},
	rich_markup_mode = 'rich',
	)

@_app.command()
def _cli(
	rawdata: Annotated[str, typer.Argument(help = "Specify the path of a rawdata input file")],
	exclude: Annotated[str, typer.Option('--exclude', '-e', help = 'The path of a file specifying UIDs and/or Samples to exclude')] = 'none',
	anchors: Annotated[str, typer.Option('--anchors', '-a', help = 'The path of a file specifying custom anchors')] = 'none',
	output_dir: Annotated[str, typer.Option('--output-dir', '-o', help = 'Specify the output directory')] = 'output',
	run_D48: Annotated[bool, typer.Option('--D48', help = 'Also standardize D48')] = False,
	):
	"""
	Process raw D47 data and return standardized results.

	See [b]https://mdaeron.github.io/D47crunch/#3-command-line-interface-cli[/b] for more details.

	Reads raw data from an input file, optionally excluding some samples and/or analyses, thean standardizes
	the data based either on the default [b]d13C_VPDB[/b], [b]d18O_VPDB[/b], [b]D47[/b], and [b]D48[/b] anchors or on different
	user-specified anchors. A new directory (named `output` by default) is created to store the results and
	the following sequence is applied:

	* [b]D47data.wg()[/b]
	* [b]D47data.crunch()[/b]
	* [b]D47data.standardize()[/b]
	* [b]D47data.summary()[/b]
	* [b]D47data.table_of_samples()[/b]
	* [b]D47data.table_of_sessions()[/b]
	* [b]D47data.plot_sessions()[/b]
	* [b]D47data.plot_residuals()[/b]
	* [b]D47data.table_of_analyses()[/b]
	* [b]D47data.plot_distribution_of_analyses()[/b]
	* [b]D47data.plot_bulk_compositions()[/b]
	* [b]D47data.save_D47_correl()[/b]

	Optionally, also apply similar methods for [b]]D48[/b].

	[b]Example CSV file for --anchors option:[/b]
	[i]
	Sample,  d13C_VPDB,  d18O_VPDB,     D47,    D48
	ETH-1,        2.02,      -2.19,  0.2052,  0.138
	ETH-2,      -10.17,     -18.69,  0.2085,  0.138
	ETH-3,        1.71,      -1.78,  0.6132,  0.270
	ETH-4,            ,           ,  0.4511,  0.223
	[/i]
	Except for [i]Sample[/i], none of the columns above are mandatory.

	[b]Example CSV file for --exclude option:[/b]
	[i]
	Sample,  UID
	 FOO-1,
	 BAR-2,
	      ,  A04
	      ,  A17
	      ,  A88
	[/i]
	This will exclude all analyses of samples [i]FOO-1[/i] and [i]BAR-2[/i],
	and the analyses with UIDs [i]A04[/i], [i]A17[/i], and [i]A88[/i].
	Neither column is mandatory.
	"""

	data = D47data()
	data.read(rawdata)

	if exclude != 'none':
		exclude = read_csv(exclude)
		exclude_uid = {r['UID'] for r in exclude if 'UID' in r}
		exclude_sample = {r['Sample'] for r in exclude if 'Sample' in r}
	else:
		exclude_uid = []
		exclude_sample = []

	data = D47data([r for r in data if r['UID'] not in exclude_uid and r['Sample'] not in exclude_sample])

	if anchors != 'none':
		anchors = read_csv(anchors)
		if len([_ for _ in anchors if 'd13C_VPDB' in _]):
			data.Nominal_d13C_VPDB = {
				_['Sample']: _['d13C_VPDB']
				for _ in anchors
				if 'd13C_VPDB' in _
				}
		if len([_ for _ in anchors if 'd18O_VPDB' in _]):
			data.Nominal_d18O_VPDB = {
				_['Sample']: _['d18O_VPDB']
				for _ in anchors
				if 'd18O_VPDB' in _
				}
		if len([_ for _ in anchors if 'D47' in _]):
			data.RMs = {
				_['Sample']: _['D47']
				for _ in anchors
				if 'D47' in _
				}

	data.refresh()
	data.wg()
	data.crunch()
	data.standardize()
	data.summary(dir = output_dir)
	data.plot_residuals(dir = output_dir, filename = 'D47_residuals.pdf', kde = True)
	data.plot_bulk_compositions(dir = output_dir + '/bulk_compositions')
	data.plot_sessions(dir = output_dir)
	data.save_D47_correl(dir = output_dir)

	if not run_D48:
		data.table_of_samples(dir = output_dir)
		data.table_of_analyses(dir = output_dir)
		data.table_of_sessions(dir = output_dir)


	if run_D48:
		data2 = D48data()
		print(rawdata)
		data2.read(rawdata)

		data2 = D48data([r for r in data2 if r['UID'] not in exclude_uid and r['Sample'] not in exclude_sample])

		if anchors != 'none':
			if len([_ for _ in anchors if 'd13C_VPDB' in _]):
				data2.Nominal_d13C_VPDB = {
					_['Sample']: _['d13C_VPDB']
					for _ in anchors
					if 'd13C_VPDB' in _
					}
			if len([_ for _ in anchors if 'd18O_VPDB' in _]):
				data2.Nominal_d18O_VPDB = {
					_['Sample']: _['d18O_VPDB']
					for _ in anchors
					if 'd18O_VPDB' in _
					}
			if len([_ for _ in anchors if 'D48' in _]):
				data2.RMs = {
					_['Sample']: _['D48']
					for _ in anchors
					if 'D48' in _
					}

		data2.refresh()
		data2.wg()
		data2.crunch()
		data2.standardize()
		data2.summary(dir = output_dir)
		data2.plot_sessions(dir = output_dir)
		data2.plot_residuals(dir = output_dir, filename = 'D48_residuals.pdf', kde = True)
		data2.plot_distribution_of_analyses(dir = output_dir)
		data2.save_D48_correl(dir = output_dir)

		table_of_analyses(data, data2, dir = output_dir)
		table_of_samples(data, data2, dir = output_dir)
		table_of_sessions(data, data2, dir = output_dir)

def __cli():
	_app()
