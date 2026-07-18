def test_bayes():

	from D47crunch import D47data, pf

	print()
	print()

	D47data.RMs['TAC-1'] = (0.700, 0.01)

	D = D47data(verbose = True)
	D.read('tests/vdata.csv')
	for s in D.sessions:
		D.sessions[s]['wg_drift'] = True

	pooled_constraints = {}
	bayes_constraints = {}
	s0 = list(D.sessions)[0]
	for s in list(D.sessions)[1:]:
		pooled_constraints[f"c2_{pf(s)}"] = f"c2_{pf(s0)}"
		bayes_constraints[f"c2['{s}']"] = f"c2['{s0}']"

	assert len(D) == 57
	assert len(D.sessions) == 3
	assert len(D.samples) == 7
	assert len(D.unknowns) == 3
	assert len(D.anchors) == 4
	assert len(D.fixed_anchors) == 3
	assert len(D.loose_anchors) == 1

	D.wg()

	for s in D.sessions:
		assert round(D.sessions[s]['d13Cwg_VPDB'], 0) == -4
		assert round(D.sessions[s]['d18Owg_VSMOW'], 0) == 26

	D.crunch()

	for r in D:
		assert 'd13C_VPDB' in r
		assert 'd18O_VSMOW' in r

	D.standardize(method = 'pooled', constraints = pooled_constraints)

	D.table_of_sessions(dir = 'tests/output')
	D.table_of_samples(dir = 'tests/output')
	D.table_of_analyses(dir = 'tests/output')

	D.plot_sessions(dir = 'tests/output')
	D.plot_residuals(dir = 'tests/output')

	D.standardize(method = 'bayes', constraints = bayes_constraints)

	D.pprint()
