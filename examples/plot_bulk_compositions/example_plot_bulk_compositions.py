#! /usr/bin/env python3

from lmfit import report_fit
from D47crunch import *

mydata = D47data(virtual_data(
	session = 'mysession',
	samples = [
		dict(Sample = 'ETH-1', N = 4),
		dict(Sample = 'ETH-2', N = 4),
		dict(Sample = 'ETH-3', N = 4),
		dict(Sample = 'MYSAMPLE', N = 8, D47 = 0.6, D48 = 0.1, d13C_VPDB = -4.0, d18O_VPDB = -12.0),
	], seed = 123))

mydata.refresh()
mydata.wg()
mydata.crunch()
mydata.plot_bulk_compositions()
