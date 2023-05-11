#! /usr/bin/env python3

from lmfit import report_fit
from D47crunch import *

mydata1 = virtual_data(
	session = 'mysession1',
	samples = [
		dict(Sample = 'ETH-1', N = 4),
		dict(Sample = 'ETH-2', N = 4),
		dict(Sample = 'ETH-3', N = 4),
		dict(Sample = 'FOO', N = 4, D47 = 0.6, D48 = 0.1, d13C_VPDB = -4.0, d18O_VPDB = -12.0),
		dict(Sample = 'BAR', N = 4, D47 = 0.5, D48 = 0.1, d13C_VPDB = -14.0, d18O_VPDB = -22.0),
	])

mydata2 = virtual_data(
	session = 'mysession2',
	samples = [
		dict(Sample = 'ETH-1', N = 4),
		dict(Sample = 'ETH-2', N = 4),
		dict(Sample = 'ETH-3', N = 4),
		dict(Sample = 'FOO', N = 4, D47 = 0.6, D48 = 0.1, d13C_VPDB = -4.0, d18O_VPDB = -12.0),
		dict(Sample = 'BAR', N = 4, D47 = 0.5, D48 = 0.1, d13C_VPDB = -14.0, d18O_VPDB = -22.0),
	])

mydata = D47data(mydata1+mydata2, verbose = True)

mydata.refresh()
mydata.wg()
mydata.crunch()
report_fit(mydata.standardize(
	constraints = {'D47_FOO': 'D47_BAR + 0.1'}
	))
mydata.table_of_sessions()
mydata.table_of_samples()
mydata.plot_sessions()
