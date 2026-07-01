from D47crunch import virtual_data, D47data, ppl

args = dict(
	samples = [
		dict(Sample = 'ETH-1', N = 4),
		dict(Sample = 'ETH-2', N = 4),
		dict(Sample = 'ETH-3', N = 1),
		dict(Sample = 'TAC-1', N = 4,
			d13C_VPDB = 6.5, d18O_VPDB = 2.0,
			D47 = 0.700, D48 = 0.3),
		dict(Sample = 'FOO', N = 3,
			d13C_VPDB = -5., d18O_VPDB = -10.,
			D47 = 0.3, D48 = 0.15),
		dict(Sample = 'BAR', N = 3,
			d13C_VPDB = 0., d18O_VPDB = -2.,
			D47 = 0.6, D48 = 0.2),
		], rD47 = 0.010, rD48 = 0.030)

D = D47data(
	virtual_data(session = 'Session_01', **args, seed = 12)
	+ virtual_data(session = 'Session_02', **args, seed = 123)
)

D.crunch()
D.standardize()

D.plot_sessions(dir = 'output/ls')
D.table_of_sessions(verbose = True, save_to_file = False)
D.table_of_samples(verbose = True, save_to_file = False)
# D.table_of_analyses(verbose = True, save_to_file = False)

D.bayesian_standardization(weak_anchors = {'TAC-1': (0.700, 0.02)})
# print(D.bayes['summary'])

D.plot_bayesian_sessions(dir = 'output/bayes')

# D.table_of_least_squares_vs_bayesian_results()
# D.plot_least_squares_vs_bayesian_results()
