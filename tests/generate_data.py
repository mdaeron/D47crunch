from D47crunch import virtual_data, D47data

args = dict(
	samples = [
		dict(Sample = 'ETH-1', N = 2),
		dict(Sample = 'ETH-2', N = 2),
		dict(Sample = 'ETH-3', N = 2),
		dict(Sample = 'TAC-1', N = 4,
			d13C_VPDB = 6.5,
			d18O_VPDB = 2.0,
			D47 = 0.700,
			D48 = 0.3,
		),
		dict(Sample = 'FOO-1', N = 3,
			d13C_VPDB = -5.,
			d18O_VPDB = -10.,
			D47 = 0.3,
			D48 = 0.15,
		),
		dict(Sample = 'BAR-1', N = 3,
			d13C_VPDB = 0.,
			d18O_VPDB = -2.,
			D47 = 0.6,
			D48 = 0.2,
		),
		dict(Sample = 'BAZ-1', N = 3,
			d13C_VPDB = -30.,
			d18O_VPDB = -17.,
			D47 = 0.65,
			D48 = 0.2,
		),
		],
		rD47 = 0.010,
		rD48 = 0.030,
)

D = (
	virtual_data(session = 'Session_01', **args, seed = 12)
	+ virtual_data(session = 'Session_02', **args, seed = 123)
	+ virtual_data(session = 'Session_03', **args, seed = 1234)
)

out = [[
	'Session',
	'Sample',
	'd45',
	'd46',
	'd47',
	'd48',
	'd49',
]]

for r in D:
	out.append([str(r[k]) for k in out[0]])

with open('vdata.csv', 'w') as f:
	f.write('\n'.join([','.join(l) for l in out]))
