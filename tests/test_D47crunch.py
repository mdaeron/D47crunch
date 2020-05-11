import D47crunch
from math import isclose

def test_correlated_sum():
	X = [1., -1.]
	C = [[0.010, -0.005], [-0.005, 0.010]]
	Y, sY = D47crunch.correlated_sum(X, C)
	assert(
		Y == 0.
		)
	assert(
		sY == 0.1
		)
		
	f = [.75, .25]
	Y, sY = D47crunch.correlated_sum(X, C, f)
	assert(
		Y == 0.5
		)
	assert(
		isclose(sY, 0.06614378277661476)
		)


def test_make_csv():
	x = [['a', 'b'], ['c', 'd']]
	y = D47crunch.make_csv(x, hsep = '-', vsep = '+')
	assert(
		y == 'a-b+c-d'
		)
	

def test_pf():
	assert(
		D47crunch.pf('a.b-c') == 'a_b_c'
		)
	

def test_smart_type():
	assert(
		isinstance(D47crunch.smart_type('1'), int)
		)
	assert(
		D47crunch.smart_type('1') == 1
		)
	assert(
		isinstance(D47crunch.smart_type('1.'), float)
		)
	assert(
		D47crunch.smart_type('1.') == 1
		)
	assert(
		isinstance(D47crunch.smart_type('foo'), str)
		)
	assert(
		D47crunch.smart_type('foo') == 'foo'
		)

