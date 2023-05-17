from D47crunch import virtual_data, D47data

args = dict(
	samples = [
		dict(Sample = 'ETH-1', N = 3),
		dict(Sample = 'ETH-2', N = 3),
		dict(Sample = 'ETH-3', N = 3),
		dict(Sample = 'FOO', N = 3,
			d13C_VPDB = -5., d18O_VPDB = -10.,
			D47 = 0.3, D48 = 0.15),
		dict(Sample = 'BAR', N = 3,
			d13C_VPDB = -15., d18O_VPDB = -2.,
			D47 = 0.6, D48 = 0.2),
		], rD47 = 0.010, rD48 = 0.030)

session1 = virtual_data(session = 'Session_01', **args, seed = 123)
session2 = virtual_data(session = 'Session_02', **args, seed = 1234)
session3 = virtual_data(session = 'Session_03', **args, seed = 12345)
session4 = virtual_data(session = 'Session_04', **args, seed = 123456)

D = D47data(session1 + session2 + session3 + session4)

D.crunch()
D.standardize()

D.table_of_sessions(verbose = True, save_to_file = False)
D.table_of_samples(verbose = True, save_to_file = False)
D.table_of_analyses(verbose = True, save_to_file = False)