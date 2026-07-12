from rich.pretty import pprint

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
		dict(Sample = 'BAR-2', N = 3,
			d13C_VPDB = 0.,
			d18O_VPDB = -2.,
			D47 = 0.6,
			D48 = 0.2,
		),
		dict(Sample = 'BAR-3', N = 3,
			d13C_VPDB = 0.,
			d18O_VPDB = -2.,
			D47 = 0.6,
			D48 = 0.2,
		),
		],
		rD47 = 0.010,
		rD48 = 0.030,
)

D47data.RMs['TAC-1'] = (0.7, 0.01)

D = D47data(
	virtual_data(session = 'Session_01', **args, seed = 12)
	+ virtual_data(session = 'Session_02', **args, seed = 123)
	+ virtual_data(session = 'Session_03', **args, seed = 1234),
	verbose = False,
)

D.crunch()
D.standardize(
	method = 'pooled',
)

D.plot_sessions(dir = 'output/pooled')
D.table_of_sessions(verbose = True, save_to_file = True, dir = 'output/pooled')
D.table_of_samples(verbose = True, save_to_file = True, dir = 'output/pooled')

D.standardize(
	method = 'bayes',
	constraints = {
		# two different ways to specify that WG bulk composition remains constant:
		"c['Session_02']": "c['Session_01'] / a['Session_01'] * a['Session_02']",
		"c[2]": "c[0] / a[0] * a[2]",
		# also specify a known D47 offset between FOO and BAR:
		"D47['FOO-1']": "D47['BAR-1'] - 0.3",
	},
)

D.plot_sessions(dir = 'output/bayes')
D.table_of_sessions(verbose = True, save_to_file = True, dir = 'output/bayes')
D.table_of_samples(verbose = True, save_to_file = True, dir = 'output/bayes')

D.table_of_least_squares_vs_bayesian_results(dir = 'output/bayes', verbose = True)
D.plot_least_squares_vs_bayesian_results(dir = 'output/bayes')
