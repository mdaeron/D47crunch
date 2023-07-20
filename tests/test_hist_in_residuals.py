#! /usr/bin/env python3

from D47crunch import *

mydata = virtual_data(
	session = 'mysession',
	samples = [
		dict(Sample = 'ETH-1', N = 40),
		dict(Sample = 'ETH-2', N = 40),
		dict(Sample = 'ETH-3', N = 40),
		dict(Sample = 'ETH-4', N = 40, d13C_VPDB = -10.2, d18O_VPDB = -18.8),
		dict(Sample = 'FOO', N = 80, D47 = 0.6, D48 = 0.1, d13C_VPDB = -4.0, d18O_VPDB = -12.0),
		dict(Sample = 'BAR', N = 80, D47 = 0.6, D48 = 0.1, d13C_VPDB = -4.0, d18O_VPDB = -12.0),
	], seed = 0)

mydata = D47data(mydata, verbose = True)

mydata.refresh()
mydata.wg()
mydata.crunch()
mydata.standardize()
mydata.plot_residuals(filename = 'residuals.pdf', hist = True)
