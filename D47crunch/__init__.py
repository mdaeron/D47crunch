'''
Standardization and analytical error propagation of Δ47 and Δ48 clumped-isotope measurements

Process and standardize carbonate and/or CO2 clumped-isotope analyses,
from low-level data out of a dual-inlet mass spectrometer to final, “absolute”
Δ47 and Δ48 values with fully propagated analytical error estimates
([Daëron, 2021](https://doi.org/10.1029/2020GC009592)).

The **tutorial** section takes you through a series of simple steps to import/process data and print out the results.
The **how-to** section provides instructions applicable to various specific tasks.

.. include:: ../docs/tutorial.md
.. include:: ../docs/howto.md
.. include:: ../docs/cli.md

# 4. API Documentation
'''

__docformat__ = "restructuredtext"
__author__    = 'Mathieu Daëron'
__contact__   = 'daeron@lsce.ipsl.fr'
__copyright__ = 'Copyright (c) 2023 Mathieu Daëron'
__license__   = 'Modified BSD License - https://opensource.org/licenses/BSD-3-Clause'
__date__      = '2023-07-22'
__version__   = '2.3.1'

import os
import numpy as np
import typer
from typing_extensions import Annotated
from statistics import stdev
from scipy.stats import t as tstudent
from scipy.stats import levene
from scipy.interpolate import interp1d
from numpy import linalg
from lmfit import Minimizer, Parameters, report_fit
from matplotlib import pyplot as ppl
from datetime import datetime as dt
from functools import wraps
from colorsys import hls_to_rgb
from matplotlib import rcParams

typer.rich_utils.STYLE_HELPTEXT = ''

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


def pretty_table(x, header = 1, hsep = '  ', vsep = '–', align = '<'):
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
	x = [['A', 'B', 'C'], ['1', '1.9999', 'foo'], ['10', 'x', 'bar']]
	print(pretty_table(x))
	```
	yields:	
	```
	--  ------  ---
	A        B    C
	--  ------  ---
	1   1.9999  foo
	10       x  bar
	--  ------  ---
	```
	
	'''
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
	Nominal_D47 = None,
	Nominal_D48 = None,
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
	+ `Nominal_D47`, `Nominal_D48`: where to lookup Δ47 and
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
	
	if Nominal_D47 is None:
		Nominal_D47 = D47data().Nominal_D47

	if Nominal_D48 is None:
		Nominal_D48 = D48data().Nominal_D48
	
	if d13C_VPDB is None:
		if sample in Nominal_d13C_VPDB:
			d13C_VPDB = Nominal_d13C_VPDB[sample]
		else:
			raise KeyError(f"Sample {sample} is missing d13C_VDP value, and it is not defined in Nominal_d13C_VDP.")

	if d18O_VPDB is None:
		if sample in Nominal_d18O_VPDB:
			d18O_VPDB = Nominal_d18O_VPDB[sample]
		else:
			raise KeyError(f"Sample {sample} is missing d18O_VPDB value, and it is not defined in Nominal_d18O_VPDB.")

	if D47 is None:
		if sample in Nominal_D47:
			D47 = Nominal_D47[sample]
		else:
			raise KeyError(f"Sample {sample} is missing D47 value, and it is not defined in Nominal_D47.")

	if D48 is None:
		if sample in Nominal_D48:
			D48 = Nominal_D48[sample]
		else:
			raise KeyError(f"Sample {sample} is missing D48 value, and it is not defined in Nominal_D48.")

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
	Nominal_D47 = None, Nominal_D48 = None,
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
	+ `Nominal_D47`, `Nominal_D48`: where to lookup Δ47 and Δ48 values
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
	.. include:: ../code_examples/virtual_data/example.py
	```
	
	This should output something like:
	
	```
	.. include:: ../code_examples/virtual_data/output.txt
	```
	'''
	
	kwargs = locals().copy()

	from numpy import random as nprandom
	if seed:
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
					'Nominal_D47', 'Nominal_D48', 'Nominal_d13C_VPDB', 'Nominal_d18O_VPDB',
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
			out47 = data47.table_of_samples(save_to_file = False, print_out = False, output = 'raw')
			out48 = data48.table_of_samples(save_to_file = False, print_out = False, output = 'raw')
			out = transpose_table(transpose_table(out47) + transpose_table(out48)[4:])

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
				output = output
				)
	else:
		if data48 is None:
			return data47.table_of_sessions(
				dir = dir,
				filename = filename,
				save_to_file = save_to_file,
				print_out = print_out,
				output = output
				)
		else:
			out47 = data47.table_of_sessions(save_to_file = False, print_out = False, output = 'raw')
			out48 = data48.table_of_sessions(save_to_file = False, print_out = False, output = 'raw')
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
					filename = f'D47D48_sessions.csv'
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

	LEVENE_REF_SAMPLE = 'ETH-3'
	'''
	After the Δ4x standardization step, each sample is tested to
	assess whether the Δ4x variance within all analyses for that
	sample differs significantly from that observed for a given reference
	sample (using [Levene's test](https://en.wikipedia.org/wiki/Levene%27s_test),
	which yields a p-value corresponding to the null hypothesis that the
	underlying variances are equal).

	`LEVENE_REF_SAMPLE` (by default equal to `'ETH-3'`) specifies which
	sample should be used as a reference for this test.
	'''

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
		self.Nf = None
		self.repeatability = {}
		self.refresh(session = session)


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
		self.anchors = {s: self.samples[s] for s in self.samples if s in self.Nominal_D4x}
		self.unknowns = {s: self.samples[s] for s in self.samples if s not in self.Nominal_D4x}


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
	def wg(self, samples = None, a18_acid = None):
		'''
		Compute bulk composition of the working gas for each session based on
		the carbonate standards defined in both `self.Nominal_d13C_VPDB` and
		`self.Nominal_d18O_VPDB`.
		'''

		self.msg('Computing WG composition:')

		if a18_acid is None:
			a18_acid = self.ALPHA_18O_ACID_REACTION
		if samples is None:
			samples = [s for s in self.Nominal_d13C_VPDB if s in self.Nominal_d18O_VPDB]

		assert a18_acid, f'Acid fractionation factor should not be zero.'

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
		
		for s in self.sessions:
			db = [r for r in self.sessions[s]['data'] if r['Sample'] in samples]
			assert db, f'No sample from {samples} found in session "{s}".'
# 			dbsamples = sorted({r['Sample'] for r in db})

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

			self.msg(f'Session {s} WG:   δ13C_VPDB = {d13Cwg_VPDB:.3f}   δ18O_VSMOW = {d18Owg_VSMOW:.3f}')

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
		constraints = {},
		):
		'''
		Compute absolute Δ4x values for all replicate analyses and for sample averages.
		If `method` argument is set to `'pooled'`, the standardization processes all sessions
		in a single step, assuming that all samples (anchors and unknowns alike) are homogeneous,
		i.e. that their true Δ4x value does not change between sessions,
		([Daëron, 2021](https://doi.org/10.1029/2020GC009592)). If `method` argument is set to
		`'indep_sessions'`, the standardization processes each session independently, based only
		on anchors analyses.
		'''

		self.standardization_method = method
		self.assign_timestamps()

		if method == 'pooled':
			if weighted_sessions:
				for session_group in weighted_sessions:
					if self._4x == '47':
						X = D47data([r for r in self if r['Session'] in session_group])
					elif self._4x == '48':
						X = D48data([r for r in self if r['Session'] in session_group])
					X.Nominal_D4x = self.Nominal_D4x.copy()
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
					if r['Sample'] in self.Nominal_D4x:
						R += [ (
							r[f'D{self._4x}raw'] - (
								p[f'a_{session}'] * self.Nominal_D4x[r['Sample']]
								+ p[f'b_{session}'] * r[f'd{self._4x}']
								+	p[f'c_{session}']
								+ r['t'] * (
									p[f'a2_{session}'] * self.Nominal_D4x[r['Sample']]
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
			self.Nf = result.nfree
			self.t95 = tstudent.ppf(1 - 0.05/2, self.Nf)
			new_names, new_covar, new_se = _fullcovar(result)[:3]
			result.var_names = new_names
			result.covar = new_covar

			for r in self:
				s = pf(r["Session"])
				a = result.params.valuesdict()[f'a_{s}']
				b = result.params.valuesdict()[f'b_{s}']
				c = result.params.valuesdict()[f'c_{s}']
				a2 = result.params.valuesdict()[f'a2_{s}']
				b2 = result.params.valuesdict()[f'b2_{s}']
				c2 = result.params.valuesdict()[f'c2_{s}']
				r[f'D{self._4x}'] = (r[f'D{self._4x}raw'] - c - b * r[f'd{self._4x}'] - c2 * r['t'] - b2 * r['t'] * r[f'd{self._4x}']) / (a + a2 * r['t'])
				

			self.standardization = result

			for session in self.sessions:
				self.sessions[session]['Np'] = 3
				for k in ['scrambling', 'slope', 'wg']:
					if self.sessions[session][f'{k}_drift']:
						self.sessions[session]['Np'] += 1

			if consolidate:
				self.consolidate(tables = consolidate_tables, plots = consolidate_plots)
			return result


		elif method == 'indep_sessions':

			if weighted_sessions:
				for session_group in weighted_sessions:
					X = D4xdata([r for r in self if r['Session'] in session_group], mass = self._4x)
					X.Nominal_D4x = self.Nominal_D4x.copy()
					X.refresh()
					# This is only done to assign r['wD47raw'] for r in X:
					X.standardize(method = method, weighted_sessions = [], consolidate = False)
					self.msg(f'D{self._4x}raw weights set to {1000*X[0][f"wD{self._4x}raw"]:.1f} ppm for sessions in {session_group}')
			else:
				self.msg('All weights set to 1 ‰')
				for r in self:
					r[f'wD{self._4x}raw'] = 1

			for session in self.sessions:
				s = self.sessions[session]
				p_names = ['a', 'b', 'c', 'a2', 'b2', 'c2']
				p_active = [True, True, True, s['scrambling_drift'], s['slope_drift'], s['wg_drift']]
				s['Np'] = sum(p_active)
				sdata = s['data']

				A = np.array([
					[
						self.Nominal_D4x[r['Sample']] / r[f'wD{self._4x}raw'],
						r[f'd{self._4x}'] / r[f'wD{self._4x}raw'],
						1 / r[f'wD{self._4x}raw'],
						self.Nominal_D4x[r['Sample']] * r['t'] / r[f'wD{self._4x}raw'],
						r[f'd{self._4x}'] * r['t'] / r[f'wD{self._4x}raw'],
						r['t'] / r[f'wD{self._4x}raw']
						]
					for r in sdata if r['Sample'] in self.anchors
					])[:,p_active] # only keep columns for the active parameters
				Y = np.array([[r[f'D{self._4x}raw'] / r[f'wD{self._4x}raw']] for r in sdata if r['Sample'] in self.anchors])
				s['Na'] = Y.size
				CM = linalg.inv(A.T @ A)
				bf = (CM @ A.T @ Y).T[0,:]
				k = 0
				for n,a in zip(p_names, p_active):
					if a:
						s[n] = bf[k]
# 						self.msg(f'{n} = {bf[k]}')
						k += 1
					else:
						s[n] = 0.
# 						self.msg(f'{n} = 0.0')

				for r in sdata :
					a, b, c, a2, b2, c2 = s['a'], s['b'], s['c'], s['a2'], s['b2'], s['c2']
					r[f'D{self._4x}'] = (r[f'D{self._4x}raw'] - c - b * r[f'd{self._4x}'] - c2 * r['t'] - b2 * r['t'] * r[f'd{self._4x}']) / (a + a2 * r['t'])
					r[f'wD{self._4x}'] = r[f'wD{self._4x}raw'] / (a + a2 * r['t'])

				s['CM'] = np.zeros((6,6))
				i = 0
				k_active = [j for j,a in enumerate(p_active) if a]
				for j,a in enumerate(p_active):
					if a:
						s['CM'][j,k_active] = CM[i,:]
						i += 1

			if not weighted_sessions:
				w = self.rmswd()['rmswd']
				for r in self:
						r[f'wD{self._4x}'] *= w
						r[f'wD{self._4x}raw'] *= w
				for session in self.sessions:
					self.sessions[session]['CM'] *= w**2

			for session in self.sessions:
				s = self.sessions[session]
				s['SE_a'] = s['CM'][0,0]**.5
				s['SE_b'] = s['CM'][1,1]**.5
				s['SE_c'] = s['CM'][2,2]**.5
				s['SE_a2'] = s['CM'][3,3]**.5
				s['SE_b2'] = s['CM'][4,4]**.5
				s['SE_c2'] = s['CM'][5,5]**.5

			if not weighted_sessions:
				self.Nf = len(self) - len(self.unknowns) - np.sum([self.sessions[s]['Np'] for s in self.sessions])
			else:
				self.Nf = 0
				for sg in weighted_sessions:
					self.Nf += self.rmswd(sessions = sg)['Nf']

			self.t95 = tstudent.ppf(1 - 0.05/2, self.Nf)

			avgD4x = {
				sample: np.mean([r[f'D{self._4x}'] for r in self if r['Sample'] == sample])
				for sample in self.samples
				}
			chi2 = np.sum([(r[f'D{self._4x}'] - avgD4x[r['Sample']])**2 for r in self])
			rD4x = (chi2/self.Nf)**.5
			self.repeatability[f'sigma_{self._4x}'] = rD4x

			if consolidate:
				self.consolidate(tables = consolidate_tables, plots = consolidate_plots)


	def standardization_error(self, session, d4x, D4x, t = 0):
		'''
		Compute standardization error for a given session and
		(δ47, Δ47) composition.
		'''
		a = self.sessions[session]['a']
		b = self.sessions[session]['b']
		c = self.sessions[session]['c']
		a2 = self.sessions[session]['a2']
		b2 = self.sessions[session]['b2']
		c2 = self.sessions[session]['c2']
		CM = self.sessions[session]['CM']

		x, y = D4x, d4x
		z = a * x + b * y + c + a2 * x * t + b2 * y * t + c2 * t
# 		x = (z - b*y - b2*y*t - c - c2*t) / (a+a2*t)
		dxdy = -(b+b2*t) / (a+a2*t)
		dxdz = 1. / (a+a2*t)
		dxda = -x / (a+a2*t)
		dxdb = -y / (a+a2*t)
		dxdc = -1. / (a+a2*t)
		dxda2 = -x * a2 / (a+a2*t)
		dxdb2 = -y * t / (a+a2*t)
		dxdc2 = -t / (a+a2*t)
		V = np.array([dxda, dxdb, dxdc, dxda2, dxdb2, dxdc2])
		sx = (V @ CM @ V.T) ** .5
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
		out += [['Student\'s 95% t-factor', f"{self.t95:.2f}"]]
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
			out += [[
				session,
				f"{self.sessions[session]['Na']}",
				f"{self.sessions[session]['Nu']}",
				f"{self.sessions[session]['d13Cwg_VPDB']:.3f}",
				f"{self.sessions[session]['d18Owg_VSMOW']:.3f}",
				f"{self.sessions[session]['r_d13C_VPDB']:.4f}",
				f"{self.sessions[session]['r_d18O_VSMOW']:.4f}",
				f"{self.sessions[session][f'r_D{self._4x}']:.4f}",
				f"{self.sessions[session]['a']:.3f} ± {self.sessions[session]['SE_a']:.3f}",
				f"{1e3*self.sessions[session]['b']:.3f} ± {1e3*self.sessions[session]['SE_b']:.3f}",
				f"{self.sessions[session]['c']:.3f} ± {self.sessions[session]['SE_c']:.3f}",
				]]
			if include_a2:
				if self.sessions[session]['scrambling_drift']:
					out[-1] += [f"{self.sessions[session]['a2']:.1e} ± {self.sessions[session]['SE_a2']:.1e}"]
				else:
					out[-1] += ['']
			if include_b2:
				if self.sessions[session]['slope_drift']:
					out[-1] += [f"{self.sessions[session]['b2']:.1e} ± {self.sessions[session]['SE_b2']:.1e}"]
				else:
					out[-1] += ['']
			if include_c2:
				if self.sessions[session]['wg_drift']:
					out[-1] += [f"{self.sessions[session]['c2']:.1e} ± {self.sessions[session]['SE_c2']:.1e}"]
				else:
					out[-1] += ['']

		if save_to_file:
			if not os.path.exists(dir):
				os.makedirs(dir)
			if filename is None:
				filename = f'D{self._4x}_sessions.csv'
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

		out = [['Sample','N','d13C_VPDB','d18O_VSMOW',f'D{self._4x}','SE','95% CL','SD','p_Levene']]
		for sample in self.anchors:
			out += [[
				f"{sample}",
				f"{self.samples[sample]['N']}",
				f"{self.samples[sample]['d13C_VPDB']:.2f}",
				f"{self.samples[sample]['d18O_VSMOW']:.2f}",
				f"{self.samples[sample][f'D{self._4x}']:.4f}",'','',
				f"{self.samples[sample][f'SD_D{self._4x}']:.4f}" if self.samples[sample]['N'] > 1 else '', ''
				]]
		for sample in self.unknowns:
			out += [[
				f"{sample}",
				f"{self.samples[sample]['N']}",
				f"{self.samples[sample]['d13C_VPDB']:.2f}",
				f"{self.samples[sample]['d18O_VSMOW']:.2f}",
				f"{self.samples[sample][f'D{self._4x}']:.4f}",
				f"{self.samples[sample][f'SE_D{self._4x}']:.4f}",
				f"± {self.samples[sample][f'SE_D{self._4x}'] * self.t95:.4f}",
				f"{self.samples[sample][f'SD_D{self._4x}']:.4f}" if self.samples[sample]['N'] > 1 else '',
				f"{self.samples[sample]['p_Levene']:.3f}" if self.samples[sample]['N'] > 2 else ''
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


	def plot_sessions(self, dir = 'output', figsize = (8,8), filetype = 'pdf', dpi = 100):
		'''
		Generate session plots and save them to disk.

		**Parameters**

		+ `dir`: the directory in which to save the plots
		+ `figsize`: the width and height (in inches) of each plot
		+ `filetype`: 'pdf' or 'png'
		+ `dpi`: resolution for PNG output
		'''
		if not os.path.exists(dir):
			os.makedirs(dir)

		for session in self.sessions:
			sp = self.plot_single_session(session, xylimits = 'constant')
			ppl.savefig(f'{dir}/D{self._4x}_plot_{session}.{filetype}', **({'dpi': dpi} if filetype.lower() == 'png' else {}))
			ppl.close(sp.fig)


	@make_verbal
	def consolidate_samples(self):
		'''
		Compile various statistics for each sample.

		For each anchor sample:

		+ `D47` or `D48`: the nominal Δ4x value for this anchor, specified by `self.Nominal_D4x`
		+ `SE_D47` or `SE_D48`: set to zero by definition

		For each unknown sample:

		+ `D47` or `D48`: the standardized Δ4x value for this unknown
		+ `SE_D47` or `SE_D48`: the standard error of Δ4x for this unknown

		For each anchor and unknown:

		+ `N`: the total number of analyses of this sample
		+ `SD_D47` or `SD_D48`: the “sample” (in the statistical sense) standard deviation for this sample
		+ `d13C_VPDB`: the average δ13C_VPDB value for this sample
		+ `d18O_VSMOW`: the average δ18O_VSMOW value for this sample (as CO2)
		+ `p_Levene`: the p-value from a [Levene test](https://en.wikipedia.org/wiki/Levene%27s_test) of equal
		variance, indicating whether the Δ4x repeatability this sample differs significantly from
		that observed for the reference sample specified by `self.LEVENE_REF_SAMPLE`.
		'''
		D4x_ref_pop = [r[f'D{self._4x}'] for r in self.samples[self.LEVENE_REF_SAMPLE]['data']]
		for sample in self.samples:
			self.samples[sample]['N'] = len(self.samples[sample]['data'])
			if self.samples[sample]['N'] > 1:
				self.samples[sample][f'SD_D{self._4x}'] = stdev([r[f'D{self._4x}'] for r in self.samples[sample]['data']])

			self.samples[sample]['d13C_VPDB'] = np.mean([r['d13C_VPDB'] for r in self.samples[sample]['data']])
			self.samples[sample]['d18O_VSMOW'] = np.mean([r['d18O_VSMOW'] for r in self.samples[sample]['data']])

			D4x_pop = [r[f'D{self._4x}'] for r in self.samples[sample]['data']]
			if len(D4x_pop) > 2:
				self.samples[sample]['p_Levene'] = levene(D4x_ref_pop, D4x_pop, center = 'median')[1]
			
		if self.standardization_method == 'pooled':
			for sample in self.anchors:
				self.samples[sample][f'D{self._4x}'] = self.Nominal_D4x[sample]
				self.samples[sample][f'SE_D{self._4x}'] = 0.
			for sample in self.unknowns:
				self.samples[sample][f'D{self._4x}'] = self.standardization.params.valuesdict()[f'D{self._4x}_{pf(sample)}']
				try:
					self.samples[sample][f'SE_D{self._4x}'] = self.sample_D4x_covar(sample)**.5
				except ValueError:
					# when `sample` is constrained by self.standardize(constraints = {...}),
					# it is no longer listed in self.standardization.var_names.
					# Temporary fix: define SE as zero for now
					self.samples[sample][f'SE_D4{self._4x}'] = 0.

		elif self.standardization_method == 'indep_sessions':
			for sample in self.anchors:
				self.samples[sample][f'D{self._4x}'] = self.Nominal_D4x[sample]
				self.samples[sample][f'SE_D{self._4x}'] = 0.
			for sample in self.unknowns:
				self.msg(f'Consolidating sample {sample}')
				self.unknowns[sample][f'session_D{self._4x}'] = {}
				session_avg = []
				for session in self.sessions:
					sdata = [r for r in self.sessions[session]['data'] if r['Sample'] == sample]
					if sdata:
						self.msg(f'{sample} found in session {session}')
						avg_D4x = np.mean([r[f'D{self._4x}'] for r in sdata])
						avg_d4x = np.mean([r[f'd{self._4x}'] for r in sdata])
						# !! TODO: sigma_s below does not account for temporal changes in standardization error
						sigma_s = self.standardization_error(session, avg_d4x, avg_D4x)
						sigma_u = sdata[0][f'wD{self._4x}raw'] / self.sessions[session]['a'] / len(sdata)**.5
						session_avg.append([avg_D4x, (sigma_u**2 + sigma_s**2)**.5])
						self.unknowns[sample][f'session_D{self._4x}'][session] = session_avg[-1]
				self.samples[sample][f'D{self._4x}'], self.samples[sample][f'SE_D{self._4x}'] = w_avg(*zip(*session_avg))
				weights = {s: self.unknowns[sample][f'session_D{self._4x}'][s][1]**-2 for s in self.unknowns[sample][f'session_D{self._4x}']}
				wsum = sum([weights[s] for s in weights])
				for s in weights:
					self.unknowns[sample][f'session_D{self._4x}'][s] += [self.unknowns[sample][f'session_D{self._4x}'][s][1]**-2 / wsum]

		for r in self:
			r[f'D{self._4x}_residual'] = r[f'D{self._4x}'] - self.samples[r['Sample']][f'D{self._4x}']



	def consolidate_sessions(self):
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

		if self.standardization_method == 'pooled':
			for session in self.sessions:

				self.sessions[session]['a'] = self.standardization.params.valuesdict()[f'a_{pf(session)}']
				i = self.standardization.var_names.index(f'a_{pf(session)}')
				self.sessions[session]['SE_a'] = self.standardization.covar[i,i]**.5

				self.sessions[session]['b'] = self.standardization.params.valuesdict()[f'b_{pf(session)}']
				i = self.standardization.var_names.index(f'b_{pf(session)}')
				self.sessions[session]['SE_b'] = self.standardization.covar[i,i]**.5

				self.sessions[session]['c'] = self.standardization.params.valuesdict()[f'c_{pf(session)}']
				i = self.standardization.var_names.index(f'c_{pf(session)}')
				self.sessions[session]['SE_c'] = self.standardization.covar[i,i]**.5

				self.sessions[session]['a2'] = self.standardization.params.valuesdict()[f'a2_{pf(session)}']
				if self.sessions[session]['scrambling_drift']:
					i = self.standardization.var_names.index(f'a2_{pf(session)}')
					self.sessions[session]['SE_a2'] = self.standardization.covar[i,i]**.5
				else:
					self.sessions[session]['SE_a2'] = 0.

				self.sessions[session]['b2'] = self.standardization.params.valuesdict()[f'b2_{pf(session)}']
				if self.sessions[session]['slope_drift']:
					i = self.standardization.var_names.index(f'b2_{pf(session)}')
					self.sessions[session]['SE_b2'] = self.standardization.covar[i,i]**.5
				else:
					self.sessions[session]['SE_b2'] = 0.

				self.sessions[session]['c2'] = self.standardization.params.valuesdict()[f'c2_{pf(session)}']
				if self.sessions[session]['wg_drift']:
					i = self.standardization.var_names.index(f'c2_{pf(session)}')
					self.sessions[session]['SE_c2'] = self.standardization.covar[i,i]**.5
				else:
					self.sessions[session]['SE_c2'] = 0.

				i = self.standardization.var_names.index(f'a_{pf(session)}')
				j = self.standardization.var_names.index(f'b_{pf(session)}')
				k = self.standardization.var_names.index(f'c_{pf(session)}')
				CM = np.zeros((6,6))
				CM[:3,:3] = self.standardization.covar[[i,j,k],:][:,[i,j,k]]
				try:
					i2 = self.standardization.var_names.index(f'a2_{pf(session)}')
					CM[3,[0,1,2,3]] = self.standardization.covar[i2,[i,j,k,i2]]
					CM[[0,1,2,3],3] = self.standardization.covar[[i,j,k,i2],i2]
					try:
						j2 = self.standardization.var_names.index(f'b2_{pf(session)}')
						CM[3,4] = self.standardization.covar[i2,j2]
						CM[4,3] = self.standardization.covar[j2,i2]
					except ValueError:
						pass
					try:
						k2 = self.standardization.var_names.index(f'c2_{pf(session)}')
						CM[3,5] = self.standardization.covar[i2,k2]
						CM[5,3] = self.standardization.covar[k2,i2]
					except ValueError:
						pass
				except ValueError:
					pass
				try:
					j2 = self.standardization.var_names.index(f'b2_{pf(session)}')
					CM[4,[0,1,2,4]] = self.standardization.covar[j2,[i,j,k,j2]]
					CM[[0,1,2,4],4] = self.standardization.covar[[i,j,k,j2],j2]
					try:
						k2 = self.standardization.var_names.index(f'c2_{pf(session)}')
						CM[4,5] = self.standardization.covar[j2,k2]
						CM[5,4] = self.standardization.covar[k2,j2]
					except ValueError:
						pass
				except ValueError:
					pass
				try:
					k2 = self.standardization.var_names.index(f'c2_{pf(session)}')
					CM[5,[0,1,2,5]] = self.standardization.covar[k2,[i,j,k,k2]]
					CM[[0,1,2,5],5] = self.standardization.covar[[i,j,k,k2],k2]
				except ValueError:
					pass

				self.sessions[session]['CM'] = CM

		elif self.standardization_method == 'indep_sessions':
			pass # Not implemented yet


	@make_verbal
	def repeatabilities(self):
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
	def consolidate(self, tables = True, plots = True):
		'''
		Collect information about samples, sessions and repeatabilities.
		'''
		self.consolidate_samples()
		self.consolidate_sessions()
		self.repeatabilities()

		if tables:
			self.summary()
			self.table_of_sessions()
			self.table_of_analyses()
			self.table_of_samples()

		if plots:
			self.plot_sessions()


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

		if key in ['D47', 'D48']:
			# Full disclosure: the definition of Nf is tricky/debatable
			G = [r for r in self if r['Sample'] in mysamples and r['Session'] in sessions]
			chisq = (np.array([r[f'{key}_residual'] for r in G])**2).sum()
			Nf = len(G)
# 			print(f'len(G) = {Nf}')
			Nf -= len([s for s in mysamples if s in self.unknowns])
# 			print(f'{len([s for s in mysamples if s in self.unknowns])} unknown samples to consider')
			for session in sessions:
				Np = len([
					_ for _ in self.standardization.params
					if (
						self.standardization.params[_].expr is not None
						and (
							(_[0] in 'abc' and _[1] == '_' and _[2:] == pf(session))
							or (_[0] in 'abc' and _[1:3] == '2_' and _[3:] == pf(session))
							)
						)
					])
# 				print(f'session {session}: {Np} parameters to consider')
				Na = len({
					r['Sample'] for r in self.sessions[session]['data']
					if r['Sample'] in self.anchors and r['Sample'] in mysamples
					})
# 				print(f'session {session}: {Na} different anchors in that session')
				Nf -= min(Np, Na)
# 			print(f'Nf = {Nf}')

# 			for sample in mysamples :
# 				X = [ r[key] for r in self if r['Sample'] == sample and r['Session'] in sessions ]
# 				if len(X) > 1 :
# 					chisq += np.sum([ (x-self.samples[sample][key])**2 for x in X ])
# 					if sample in self.unknowns:
# 						Nf += len(X) - 1
# 					else:
# 						Nf += len(X)
# 			if samples in ['anchors', 'all samples']:
# 				Nf -= sum([self.sessions[s]['Np'] for s in sessions])
			r = (chisq / Nf)**.5 if Nf > 0 else 0

		else: # if key not in ['D47', 'D48']
			chisq, Nf = 0, 0
			for sample in mysamples :
				X = [ r[key] for r in self if r['Sample'] == sample and r['Session'] in sessions ]
				if len(X) > 1 :
					Nf += len(X) - 1
					chisq += np.sum([ (x-np.mean(X))**2 for x in X ])
			r = (chisq / Nf)**.5 if Nf > 0 else 0

		self.msg(f'Repeatability of r["{key}"] is {1000*r:.1f} ppm for {samples}.')
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


	def sample_D4x_covar(self, sample1, sample2 = None):
		'''
		Covariance between Δ4x values of samples

		Returns the error covariance between the average Δ4x values of two
		samples. If if only `sample_1` is specified, or if `sample_1 == sample_2`),
		returns the Δ4x variance for that sample.
		'''
		if sample2 is None:
			sample2 = sample1
		if self.standardization_method == 'pooled':
			i = self.standardization.var_names.index(f'D{self._4x}_{pf(sample1)}')
			j = self.standardization.var_names.index(f'D{self._4x}_{pf(sample2)}')
			return self.standardization.covar[i, j]
		elif self.standardization_method == 'indep_sessions':
			if sample1 == sample2:
				return self.samples[sample1][f'SE_D{self._4x}']**2
			else:
				c = 0
				for session in self.sessions:
					sdata1 = [r for r in self.sessions[session]['data'] if r['Sample'] == sample1]
					sdata2 = [r for r in self.sessions[session]['data'] if r['Sample'] == sample2]
					if sdata1 and sdata2:
						a = self.sessions[session]['a']
						# !! TODO: CM below does not account for temporal changes in standardization parameters
						CM = self.sessions[session]['CM'][:3,:3]
						avg_D4x_1 = np.mean([r[f'D{self._4x}'] for r in sdata1])
						avg_d4x_1 = np.mean([r[f'd{self._4x}'] for r in sdata1])
						avg_D4x_2 = np.mean([r[f'D{self._4x}'] for r in sdata2])
						avg_d4x_2 = np.mean([r[f'd{self._4x}'] for r in sdata2])
						c += (
							self.unknowns[sample1][f'session_D{self._4x}'][session][2]
							* self.unknowns[sample2][f'session_D{self._4x}'][session][2]
							* np.array([[avg_D4x_1, avg_d4x_1, 1]])
							@ CM
							@ np.array([[avg_D4x_2, avg_d4x_2, 1]]).T
							) / a**2
				return float(c)

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

	def plot_single_session(self,
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
		
		if fig == 'new':
			out.fig = ppl.figure(figsize = (6,6))
			ppl.subplots_adjust(.1,.1,.9,.9)

		out.anchor_analyses, = ppl.plot(
			[r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in self.anchors],
			[r[f'D{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in self.anchors],
			**kw_plot_anchors)
		out.unknown_analyses, = ppl.plot(
			[r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in self.unknowns],
			[r[f'D{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] in self.unknowns],
			**kw_plot_unknowns)
		out.anchor_avg = ppl.plot(
			np.array([ np.array([
				np.min([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) - 1,
				np.max([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) + 1
				]) for sample in anchors]).T,
			np.array([ np.array([0, 0]) + self.Nominal_D4x[sample] for sample in anchors]).T,
			**kw_plot_anchor_avg)
		out.unknown_avg = ppl.plot(
			np.array([ np.array([
				np.min([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) - 1,
				np.max([r[f'd{self._4x}'] for r in self.sessions[session]['data'] if r['Sample'] == sample]) + 1
				]) for sample in unknowns]).T,
			np.array([ np.array([0, 0]) + self.unknowns[sample][f'D{self._4x}'] for sample in unknowns]).T,
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
			SI = np.array([[self.standardization_error(session, x, y) for x in xi] for y in yi])
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
		):
		'''
		Plot residuals of each analysis as a function of time (actually, as a function of
		the order of analyses in the `D4xdata` object)

		+ `kde`: whether to add a kernel density estimate of residuals
		+ `hist`: whether to add a histogram of residuals (incompatible with `kde`)
		+ `histbins`: specify bin edges for the histogram
		+ `dir`: the directory in which to save the plot
		+ `highlight`: a list of samples to highlight
		+ `colors`: a dict of `{<sample>: <color>}` for all samples
		+ `figsize`: (width, height) of figure
		+ `dpi`: resolution for PNG output
		+ `yspan`: factor controlling the range of y values shown in plot
		  (by default: `yspan = 1.5 if kde else 1.0`)
		'''
		
		from matplotlib import ticker

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
		ppl.axhspan(-self.repeatability[f'r_D{self._4x}']*1000*self.t95, self.repeatability[f'r_D{self._4x}']*1000*self.t95, color = 'k', alpha = .05, lw = 1)
		if not (hist or kde):
			ppl.text(len(self), self.repeatability[f'r_D{self._4x}']*1000, f"   SD = {self.repeatability[f'r_D{self._4x}']*1000:.1f} ppm", size = 9, alpha = 1, va = 'center')
			ppl.text(len(self), self.repeatability[f'r_D{self._4x}']*1000*self.t95, f"   95% CL = ± {self.repeatability[f'r_D{self._4x}']*1000*self.t95:.1f} ppm", size = 9, alpha = 1, va = 'center')

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
				f"   SD = {self.repeatability[f'r_D{self._4x}']*1000:.1f} ppm\n   95% CL = ± {self.repeatability[f'r_D{self._4x}']*1000*self.t95:.1f} ppm",
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

		if not os.path.exists(dir):
			os.makedirs(dir)
		if filename is None:
			return fig
		elif filename == '':
			filename = f'D{self._4x}_residuals.pdf'
		ppl.savefig(f'{dir}/{filename}', dpi = dpi)
		ppl.close(fig)
				

	def simulate(self, *args, **kwargs):
		'''
		Legacy function with warning message pointing to `virtual_data()`
		'''
		raise DeprecationWarning('D4xdata.simulate is deprecated and has been replaced by virtual_data()')

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
		):
		'''
		Save D4x values along with their SE and correlation matrix.

		**Parameters**

		+ `samples`: Only these samples are output (by default: all samples).
		+ `dir`: the directory in which to save the faile (by defaut: `output`)
		+ `filename`: the name to the csv file to write to (by default: `D4x_correl.csv`)
		+ `D4x_precision`: the precision to use when writing `D4x` and `D4x_SE` values (by default: 4)
		+ `correl_precision`: the precision to use when writing correlation factor values (by default: 4)
		'''
		if samples is None:
			samples = sorted([s for s in self.unknowns])
		
		out = [['Sample']] + [[s] for s in samples]
		out[0] += [f'D{self._4x}', f'D{self._4x}_SE', f'D{self._4x}_correl']
		for k,s in enumerate(samples):
			out[k+1] += [f'{self.samples[s][f"D{self._4x}"]:.4f}', f'{self.samples[s][f"SE_D{self._4x}"]:.4f}']
			for s2 in samples:
				out[k+1] += [f'{self.sample_D4x_correl(s,s2):.4f}']
		
		if not os.path.exists(dir):
			os.makedirs(dir)
		if filename is None:
			filename = f'D{self._4x}_correl.csv'
		with open(f'{dir}/{filename}', 'w') as fid:
			fid.write(make_csv(out))
		
		
		

class D47data(D4xdata):
	'''
	Store and process data for a large set of Δ47 analyses,
	usually comprising more than one analytical session.
	'''

	Nominal_D4x = {
		'ETH-1':   0.2052,
		'ETH-2':   0.2085,
		'ETH-3':   0.6132,
		'ETH-4':   0.4511,
		'IAEA-C1': 0.3018,
		'IAEA-C2': 0.6409,
		'MERCK':   0.5135,
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


	@property
	def Nominal_D47(self):
		return self.Nominal_D4x
	

	@Nominal_D47.setter
	def Nominal_D47(self, new):
		self.Nominal_D4x = dict(**new)
		self.refresh()


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
			self.Nominal_D47 = {}
		for s in foo:
			if priority != 'old' or s not in self.Nominal_D47:
				self.Nominal_D47[s] = foo[s]
	
	def save_D47_correl(self, *args, **kwargs):
		return self._save_D4x_correl(*args, **kwargs)

	save_D47_correl.__doc__ = D4xdata._save_D4x_correl.__doc__.replace('D4x', 'D47')


class D48data(D4xdata):
	'''
	Store and process data for a large set of Δ48 analyses,
	usually comprising more than one analytical session.
	'''

	Nominal_D4x = {
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
	Fiebig et al. (in press)):

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


	@property
	def Nominal_D48(self):
		return self.Nominal_D4x

	
	@Nominal_D48.setter
	def Nominal_D48(self, new):
		self.Nominal_D4x = dict(**new)
		self.refresh()


	def __init__(self, l = [], **kwargs):
		'''
		**Parameters:** same as `D4xdata.__init__()`
		'''
		D4xdata.__init__(self, l = l, mass = '48', **kwargs)

	def save_D48_correl(self, *args, **kwargs):
		return self._save_D4x_correl(*args, **kwargs)

	save_D48_correl.__doc__ = D4xdata._save_D4x_correl.__doc__.replace('D4x', 'D48')

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
	the data based either on the default [b]d13C_VDPB[/b], [b]d18O_VPDB[/b], [b]D47[/b], and [b]D48[/b] anchors or on different
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
			data.Nominal_D4x = {
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
				data2.Nominal_D4x = {
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
