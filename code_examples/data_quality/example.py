from D47crunch import *
from random import shuffle

# generate virtual data:
args = dict(
	samples = [
		dict(Sample = 'ETH-1', N = 8),
		dict(Sample = 'ETH-2', N = 8),
		dict(Sample = 'ETH-3', N = 8),
		dict(Sample = 'FOO', N = 4,
			d13C_VPDB = -5., d18O_VPDB = -10.,
			D47 = 0.3, D48 = 0.15),
		dict(Sample = 'BAR', N = 4,
			d13C_VPDB = -15., d18O_VPDB = -15.,
			D47 = 0.5, D48 = 0.2),
		])

sessions = [
	virtual_data(session = f'Session_{k+1:02.0f}', seed = 123456+k, **args)
	for k in range(10)]

# shuffle the data:
data = [r for s in sessions for r in s]
shuffle(data)
data = sorted(data, key = lambda r: r['Session'])

# create D47data instance:
data47 = D47data(data)

# process D47data instance:
data47.crunch()
data47.standardize()