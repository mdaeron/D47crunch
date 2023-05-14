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
	], seed = 123)

mydata2 = virtual_data(
	session = 'mysession2',
	samples = [
		dict(Sample = 'ETH-1', N = 4),
		dict(Sample = 'ETH-2', N = 4),
		dict(Sample = 'ETH-3', N = 4),
		dict(Sample = 'FOO', N = 4, D47 = 0.6, D48 = 0.1, d13C_VPDB = -4.0, d18O_VPDB = -12.0),
		dict(Sample = 'BAR', N = 4, D47 = 0.5, D48 = 0.1, d13C_VPDB = -14.0, d18O_VPDB = -22.0),
	], seed = 456)

mydata = D47data(mydata1+mydata2, verbose = True)

mydata.refresh()
mydata.wg()
mydata.crunch()
mydata.plot_bulk_compositions()
