import D47crunch
from math import isclose

def test_fCO2eq():
	assert(D47crunch.fCO2eqD47_Petersen(25.5) == 0.9169628945)
	assert(D47crunch.fCO2eqD47_Petersen(999.5) == 0.0266157117)
	assert(D47crunch.fCO2eqD47_Wang(25.5) == 0.91979)
	assert(D47crunch.fCO2eqD47_Wang(999.5) == 0.027241)


def test_correlated_sum():
	X = [1., -1.]
	C = [[0.010, -0.005], [-0.005, 0.010]]
	Y, sY = D47crunch.correlated_sum(X, C)
	assert(Y == 0.)
	assert(sY == 0.1)
	w = [.75, .25]
	Y, sY = D47crunch.correlated_sum(X, C, w)
	assert(Y == 0.5)
	assert(isclose(sY, 0.06614378277661476))


def test_make_csv():
	x = [['a', 'b'], ['c', 'd']]
	y = D47crunch.make_csv(x, hsep = '-', vsep = '+')
	assert(y == 'a-b+c-d')


def test_pf():
	assert(D47crunch.pf('a.b-c d') == 'a_b_c_d')


def test_smart_type():
	assert(isinstance(D47crunch.smart_type('1'), int))
	assert(D47crunch.smart_type('1') == 1)
	assert(isinstance(D47crunch.smart_type('1.'), float))
	assert(D47crunch.smart_type('1.') == 1)
	assert(isinstance(D47crunch.smart_type('foo'), str))
	assert(D47crunch.smart_type('foo') == 'foo')


def test_pretty_table():
	x = [['a','bb','ccc'],['ddd', 'ee', 'f'],['g', 'h', 'i']]
	assert(
		D47crunch.pretty_table(x, header = 1, hsep = '_', vsep = '+', align = '<')
		== '+++_++_+++\na  _bb_ccc\n+++_++_+++\nddd_ee_  f\ng  _ h_  i\n+++_++_+++\n'
		)
	assert(
		D47crunch.pretty_table(x, header = 1, hsep = '_', vsep = '+', align = '>')
		== '+++_++_+++\n  a_bb_ccc\n+++_++_+++\nddd_ee_  f\n  g_ h_  i\n+++_++_+++\n'
		)
	assert(
		D47crunch.pretty_table(x, header = 1, hsep = '_', vsep = '+', align = '')
		== '+++_++_+++\n  a_bb_ccc\n+++_++_+++\nddd_ee_  f\n  g_ h_  i\n+++_++_+++\n'
		)


def test_transpose_table():
	x = [['a','bb','ccc'],['ddd', 'ee', 'f'],['g', 'h', 'i']]
	y = [['a','ddd','g'],['bb', 'ee', 'h'],['ccc', 'f', 'i']]
	assert(D47crunch.transpose_table(x) == y)


def test_w_avg():
	x = [0., 10., 25.]
	sx = [1., 2., 3.]
	assert(D47crunch.w_avg(x, sx) == (3.877551020408163, 0.8571428571428571))


def test_D47data_init():
	N = 8
	l = [dict(Sample = '', d45 = 0., d46 = 0., d47 = 0.) for k in range(N)]
	for k in range(len(l)):
		l[k]['Sample'] = f'ETH-{1+k%4}'
	x = D47crunch.D47data(l)
	assert(len(x) == len(l))
	assert([s for s in x.samples] == [f'ETH-{k+1}' for k in range(4)])
	assert([s for s in x.anchors] == sorted({r['Sample'] for r in x if r['Sample'] in x.Nominal_D47}))
	assert([s for s in x.unknowns] == sorted({r['Sample'] for r in x if r['Sample'] not in x.Nominal_D47}))


def test_D47data_standardize():
	rawdata_input_str = '''UID\tSession\tSample\td45\td46\td47\td13Cwg_VPDB\td18Owg_VSMOW
A01\tSession01\tETH-1\t5.795017\t11.627668\t16.893512\t-3.75\t25.13
A02\tSession01\tIAEA-C1\t6.219070\t11.491072\t17.277490\t-3.75\t25.13
A03\tSession01\tETH-2\t-6.058681\t-4.817179\t-11.635064\t-3.75\t25.13
A04\tSession01\tIAEA-C2\t-3.861839\t4.941839\t0.606117\t-3.75\t25.13
A05\tSession01\tETH-3\t5.543654\t12.052277\t17.405548\t-3.75\t25.13
A06\tSession01\tMERCK\t-35.929352\t-2.087501\t-39.548484\t-3.75\t25.13
A07\tSession01\tETH-4\t-6.222218\t-5.194170\t-11.944111\t-3.75\t25.13
A08\tSession01\tETH-2\t-6.067055\t-4.877104\t-11.699265\t-3.75\t25.13
A09\tSession01\tMERCK\t-35.930739\t-2.080798\t-39.545632\t-3.75\t25.13
A10\tSession01\tETH-1\t5.788207\t11.559104\t16.801908\t-3.75\t25.13
A11\tSession01\tETH-4\t-6.217508\t-5.221407\t-11.987503\t-3.75\t25.13
A12\tSession01\tIAEA-C2\t-3.876921\t4.868892\t0.521845\t-3.75\t25.13
A13\tSession01\tETH-3\t5.539840\t12.013444\t17.368631\t-3.75\t25.13
A14\tSession01\tIAEA-C1\t6.219046\t11.447846\t17.234280\t-3.75\t25.13
A15\tSession01\tMERCK\t-35.932060\t-2.088659\t-39.531627\t-3.75\t25.13
A16\tSession01\tETH-3\t5.516658\t11.978320\t17.295740\t-3.75\t25.13
A17\tSession01\tETH-4\t-6.223370\t-5.253980\t-12.025298\t-3.75\t25.13
A18\tSession01\tETH-2\t-6.069734\t-4.868368\t-11.688559\t-3.75\t25.13
A19\tSession01\tIAEA-C1\t6.213642\t11.465109\t17.244547\t-3.75\t25.13
A20\tSession01\tETH-1\t5.789982\t11.535603\t16.789811\t-3.75\t25.13
A21\tSession01\tETH-4\t-6.205703\t-5.144529\t-11.909160\t-3.75\t25.13
A22\tSession01\tIAEA-C1\t6.212646\t11.406548\t17.187214\t-3.75\t25.13
A23\tSession01\tETH-3\t5.531413\t11.976697\t17.332700\t-3.75\t25.13
A24\tSession01\tMERCK\t-35.926347\t-2.124579\t-39.582201\t-3.75\t25.13
A25\tSession01\tETH-1\t5.786979\t11.527864\t16.775547\t-3.75\t25.13
A26\tSession01\tIAEA-C2\t-3.866505\t4.874630\t0.525332\t-3.75\t25.13
A27\tSession01\tETH-2\t-6.076302\t-4.922424\t-11.753283\t-3.75\t25.13
A28\tSession01\tIAEA-C2\t-3.878438\t4.818588\t0.467595\t-3.75\t25.13
A29\tSession01\tETH-3\t5.546458\t12.133931\t17.501646\t-3.75\t25.13
A30\tSession01\tETH-1\t5.802916\t11.642685\t16.904286\t-3.75\t25.13
A31\tSession01\tETH-2\t-6.069274\t-4.847919\t-11.677722\t-3.75\t25.13
A32\tSession01\tETH-3\t5.523018\t12.007363\t17.362080\t-3.75\t25.13
A33\tSession01\tETH-1\t5.802333\t11.616032\t16.884255\t-3.75\t25.13
A34\tSession01\tETH-3\t5.537375\t12.000263\t17.350856\t-3.75\t25.13
A35\tSession01\tETH-2\t-6.060713\t-4.893088\t-11.728465\t-3.75\t25.13
A36\tSession01\tETH-3\t5.532342\t11.990022\t17.342273\t-3.75\t25.13
A37\tSession01\tETH-3\t5.533622\t11.980853\t17.342245\t-3.75\t25.13
A38\tSession01\tIAEA-C2\t-3.867587\t4.893554\t0.540404\t-3.75\t25.13
A39\tSession01\tIAEA-C1\t6.201760\t11.406628\t17.189625\t-3.75\t25.13
A40\tSession01\tETH-1\t5.802150\t11.563414\t16.836189\t-3.75\t25.13
A41\tSession01\tETH-2\t-6.068598\t-4.897545\t-11.722343\t-3.75\t25.13
A42\tSession01\tMERCK\t-35.928359\t-2.098440\t-39.577150\t-3.75\t25.13
A43\tSession01\tETH-4\t-6.219175\t-5.168031\t-11.936923\t-3.75\t25.13
A44\tSession01\tIAEA-C2\t-3.871671\t4.871517\t0.518290\t-3.75\t25.13
B01\tSession02\tETH-1\t5.800180\t11.640916\t16.939044\t-3.74\t25.16
B02\tSession02\tETH-1\t5.799584\t11.631297\t16.917656\t-3.74\t25.16
B03\tSession02\tIAEA-C1\t6.225135\t11.512637\t17.335876\t-3.74\t25.16
B04\tSession02\tETH-2\t-6.030415\t-4.746444\t-11.525506\t-3.74\t25.16
B05\tSession02\tIAEA-C2\t-3.837017\t4.992780\t0.675292\t-3.74\t25.16
B06\tSession02\tETH-3\t5.536997\t12.048918\t17.420228\t-3.74\t25.16
B07\tSession02\tMERCK\t-35.928379\t-2.105615\t-39.594573\t-3.74\t25.16
B08\tSession02\tETH-4\t-6.218801\t-5.185168\t-11.964407\t-3.74\t25.16
B09\tSession02\tETH-2\t-6.068197\t-4.840037\t-11.686296\t-3.74\t25.16
B10\tSession02\tMERCK\t-35.926951\t-2.071047\t-39.546767\t-3.74\t25.16
B11\tSession02\tETH-1\t5.782634\t11.571818\t16.835185\t-3.74\t25.16
B12\tSession02\tETH-2\t-6.070168\t-4.877700\t-11.703876\t-3.74\t25.16
B13\tSession02\tETH-4\t-6.214873\t-5.190550\t-11.967040\t-3.74\t25.16
B14\tSession02\tIAEA-C2\t-3.853550\t4.919425\t0.584634\t-3.74\t25.16
B15\tSession02\tETH-3\t5.522265\t12.011737\t17.368407\t-3.74\t25.16
B16\tSession02\tIAEA-C1\t6.219374\t11.447014\t17.264258\t-3.74\t25.16
B17\tSession02\tMERCK\t-35.927733\t-2.103033\t-39.603494\t-3.74\t25.16
B18\tSession02\tETH-3\t5.527002\t11.984062\t17.332660\t-3.74\t25.16
B19\tSession02\tIAEA-C2\t-3.850358\t4.889230\t0.562794\t-3.74\t25.16
B20\tSession02\tETH-4\t-6.222398\t-5.263817\t-12.033650\t-3.74\t25.16
B21\tSession02\tETH-3\t5.525478\t11.970096\t17.340498\t-3.74\t25.16
B22\tSession02\tETH-2\t-6.070129\t-4.941487\t-11.773824\t-3.74\t25.16
B23\tSession02\tIAEA-C1\t6.217001\t11.434152\t17.232308\t-3.74\t25.16
B24\tSession02\tETH-1\t5.793421\t11.533191\t16.810838\t-3.74\t25.16
B25\tSession02\tETH-4\t-6.217740\t-5.198048\t-11.977179\t-3.74\t25.16
B26\tSession02\tIAEA-C1\t6.216912\t11.425200\t17.234224\t-3.74\t25.16
B27\tSession02\tETH-3\t5.522238\t11.932174\t17.286903\t-3.74\t25.16
B28\tSession02\tMERCK\t-35.914404\t-2.133955\t-39.614612\t-3.74\t25.16
B29\tSession02\tETH-1\t5.784156\t11.517244\t16.786548\t-3.74\t25.16
B30\tSession02\tIAEA-C2\t-3.852750\t4.884339\t0.551587\t-3.74\t25.16
B31\tSession02\tETH-2\t-6.068631\t-4.924103\t-11.764507\t-3.74\t25.16
B32\tSession02\tETH-4\t-6.220238\t-5.231375\t-12.009300\t-3.74\t25.16
B33\tSession02\tIAEA-C2\t-3.855245\t4.866571\t0.534914\t-3.74\t25.16
B34\tSession02\tETH-1\t5.788790\t11.544306\t16.809117\t-3.74\t25.16
B35\tSession02\tMERCK\t-35.935017\t-2.173682\t-39.664046\t-3.74\t25.16
B36\tSession02\tETH-3\t5.518320\t11.955048\t17.300668\t-3.74\t25.16
B37\tSession02\tETH-1\t5.790564\t11.521174\t16.781304\t-3.74\t25.16
B38\tSession02\tETH-4\t-6.218809\t-5.205256\t-11.979998\t-3.74\t25.16
B39\tSession02\tIAEA-C1\t6.204774\t11.391335\t17.181310\t-3.74\t25.16
B40\tSession02\tETH-2\t-6.076424\t-4.967973\t-11.815466\t-3.74\t25.16
C01\tSession03\tETH-3\t5.541868\t12.129615\t17.503738\t-3.74\t25.16
C02\tSession03\tETH-3\t5.534395\t12.034601\t17.391274\t-3.74\t25.16
C03\tSession03\tETH-1\t5.797568\t11.563575\t16.857871\t-3.74\t25.16
C04\tSession03\tETH-3\t5.529415\t11.969512\t17.342673\t-3.74\t25.16
C05\tSession03\tETH-1\t5.794026\t11.526540\t16.806934\t-3.74\t25.16
C06\tSession03\tETH-3\t5.527210\t11.937462\t17.294015\t-3.74\t25.16
C07\tSession03\tIAEA-C1\t6.220521\t11.430197\t17.242458\t-3.74\t25.16
C08\tSession03\tETH-2\t-6.064061\t-4.900852\t-11.732976\t-3.74\t25.16
C09\tSession03\tIAEA-C2\t-3.846482\t4.889242\t0.558395\t-3.74\t25.16
C10\tSession03\tETH-1\t5.789644\t11.520663\t16.795837\t-3.74\t25.16
C11\tSession03\tETH-4\t-6.219385\t-5.258604\t-12.036476\t-3.74\t25.16
C12\tSession03\tMERCK\t-35.936631\t-2.161769\t-39.693775\t-3.74\t25.16
C13\tSession03\tETH-2\t-6.076357\t-4.939912\t-11.803553\t-3.74\t25.16
C14\tSession03\tIAEA-C2\t-3.862518\t4.850015\t0.499777\t-3.74\t25.16
C15\tSession03\tETH-3\t5.515822\t11.928316\t17.287739\t-3.74\t25.16
C16\tSession03\tETH-4\t-6.216625\t-5.252914\t-12.033781\t-3.74\t25.16
C17\tSession03\tETH-1\t5.792540\t11.537788\t16.801906\t-3.74\t25.16
C18\tSession03\tIAEA-C1\t6.218853\t11.447394\t17.270859\t-3.74\t25.16
C19\tSession03\tETH-2\t-6.070107\t-4.944520\t-11.806885\t-3.74\t25.16
C20\tSession03\tMERCK\t-35.935001\t-2.155577\t-39.675070\t-3.74\t25.16
C21\tSession03\tETH-3\t5.542309\t12.082338\t17.471951\t-3.74\t25.16
C22\tSession03\tETH-4\t-6.209017\t-5.137393\t-11.920935\t-3.74\t25.16
C23\tSession03\tETH-1\t5.796781\t11.621197\t16.905496\t-3.74\t25.16
C24\tSession03\tMERCK\t-35.926449\t-2.053921\t-39.576918\t-3.74\t25.16
C25\tSession03\tETH-2\t-6.057158\t-4.797641\t-11.644824\t-3.74\t25.16
C26\tSession03\tIAEA-C1\t6.221982\t11.501725\t17.321709\t-3.74\t25.16
C27\tSession03\tETH-3\t5.535162\t12.023486\t17.396560\t-3.74\t25.16
C28\tSession03\tIAEA-C2\t-3.836934\t4.984196\t0.665651\t-3.74\t25.16
C29\tSession03\tETH-3\t5.531331\t11.991300\t17.353622\t-3.74\t25.16
C30\tSession03\tIAEA-C2\t-3.844008\t4.926554\t0.601156\t-3.74\t25.16
C31\tSession03\tETH-2\t-6.063163\t-4.907454\t-11.765065\t-3.74\t25.16
C32\tSession03\tMERCK\t-35.941566\t-2.163022\t-39.704731\t-3.74\t25.16
C33\tSession03\tETH-3\t5.523894\t11.992718\t17.363902\t-3.74\t25.16
C34\tSession03\tIAEA-C1\t6.220801\t11.462090\t17.282153\t-3.74\t25.16
C35\tSession03\tETH-1\t5.794369\t11.563017\t16.845673\t-3.74\t25.16
C36\tSession03\tETH-4\t-6.221257\t-5.272969\t-12.055444\t-3.74\t25.16
C37\tSession03\tETH-3\t5.517832\t11.957180\t17.312487\t-3.74\t25.16
C38\tSession03\tETH-2\t-6.053330\t-4.909476\t-11.740852\t-3.74\t25.16
C39\tSession03\tIAEA-C1\t6.217139\t11.440085\t17.244787\t-3.74\t25.16
C40\tSession03\tETH-1\t5.794091\t11.541948\t16.826158\t-3.74\t25.16
C41\tSession03\tIAEA-C2\t-3.803466\t4.894953\t0.624184\t-3.74\t25.16
C42\tSession03\tETH-3\t5.513788\t11.933062\t17.286883\t-3.74\t25.16
C43\tSession03\tETH-1\t5.793334\t11.569668\t16.844535\t-3.74\t25.16
C44\tSession03\tETH-2\t-6.064928\t-4.935031\t-11.786336\t-3.74\t25.16
C45\tSession03\tETH-4\t-6.216796\t-5.300373\t-12.075033\t-3.74\t25.16
C46\tSession03\tETH-3\t5.521772\t11.933713\t17.283775\t-3.74\t25.16
C47\tSession03\tMERCK\t-35.937762\t-2.181553\t-39.739636\t-3.74\t25.16
D01\tSession04\tETH-4\t-6.218867\t-5.242334\t-12.032129\t-3.74\t25.15
D02\tSession04\tIAEA-C1\t6.218458\t11.435622\t17.238776\t-3.74\t25.15
D03\tSession04\tETH-3\t5.522006\t11.946540\t17.300601\t-3.74\t25.15
D04\tSession04\tMERCK\t-35.931765\t-2.175265\t-39.716152\t-3.74\t25.15
D05\tSession04\tETH-1\t5.786884\t11.560397\t16.823187\t-3.74\t25.15
D06\tSession04\tIAEA-C2\t-3.846071\t4.861980\t0.534465\t-3.74\t25.15
D07\tSession04\tETH-2\t-6.072653\t-4.917987\t-11.786215\t-3.74\t25.15
D08\tSession04\tETH-3\t5.516592\t11.923729\t17.275641\t-3.74\t25.15
D09\tSession04\tETH-1\t5.789889\t11.531354\t16.804221\t-3.74\t25.15
D10\tSession04\tIAEA-C2\t-3.845074\t4.865635\t0.546284\t-3.74\t25.15
D11\tSession04\tETH-1\t5.795006\t11.507829\t16.772751\t-3.74\t25.15
D12\tSession04\tETH-1\t5.791371\t11.540606\t16.822704\t-3.74\t25.15
D13\tSession04\tETH-2\t-6.074029\t-4.937379\t-11.786614\t-3.74\t25.15
D14\tSession04\tETH-4\t-6.216977\t-5.273352\t-12.057294\t-3.74\t25.15
D15\tSession04\tIAEA-C1\t6.214304\t11.412869\t17.227005\t-3.74\t25.15
D16\tSession04\tETH-2\t-6.071021\t-4.966406\t-11.812116\t-3.74\t25.15
D17\tSession04\tETH-3\t5.543181\t12.065648\t17.455042\t-3.74\t25.15
D18\tSession04\tETH-1\t5.805793\t11.632212\t16.937561\t-3.74\t25.15
D19\tSession04\tIAEA-C1\t6.230425\t11.518038\t17.342943\t-3.74\t25.15
D20\tSession04\tETH-2\t-6.049292\t-4.811109\t-11.639895\t-3.74\t25.15
D21\tSession04\tIAEA-C2\t-3.829436\t4.967992\t0.665451\t-3.74\t25.15
D22\tSession04\tETH-3\t5.538827\t12.064780\t17.438156\t-3.74\t25.15
D23\tSession04\tMERCK\t-35.935604\t-2.092229\t-39.632228\t-3.74\t25.15
D24\tSession04\tETH-4\t-6.215430\t-5.166894\t-11.939419\t-3.74\t25.15
D25\tSession04\tETH-2\t-6.068214\t-4.868420\t-11.716099\t-3.74\t25.15
D26\tSession04\tMERCK\t-35.918898\t-2.041585\t-39.566777\t-3.74\t25.15
D27\tSession04\tETH-1\t5.786924\t11.584138\t16.861248\t-3.74\t25.15
D28\tSession04\tETH-2\t-6.062115\t-4.820423\t-11.664703\t-3.74\t25.15
D29\tSession04\tETH-4\t-6.210819\t-5.160997\t-11.943417\t-3.74\t25.15
D30\tSession04\tIAEA-C2\t-3.842542\t4.937635\t0.603831\t-3.74\t25.15
D31\tSession04\tETH-3\t5.527648\t11.985083\t17.353603\t-3.74\t25.15
D32\tSession04\tIAEA-C1\t6.221429\t11.481788\t17.284825\t-3.74\t25.15
D33\tSession04\tMERCK\t-35.922066\t-2.113682\t-39.642962\t-3.74\t25.15
D34\tSession04\tETH-3\t5.521955\t11.989323\t17.345179\t-3.74\t25.15
D35\tSession04\tIAEA-C2\t-3.838229\t4.937180\t0.617586\t-3.74\t25.15
D36\tSession04\tETH-4\t-6.215638\t-5.221584\t-11.999819\t-3.74\t25.15
D37\tSession04\tETH-2\t-6.067508\t-4.893477\t-11.754488\t-3.74\t25.15
D38\tSession04\tIAEA-C1\t6.214580\t11.440629\t17.254051\t-3.74\t25.15'''
	rawdata = D47crunch.D47data()
	rawdata.input(rawdata_input_str)
	rawdata.R13_VPDB = 0.01118
	rawdata.R18_VSMOW = 0.0020052
	rawdata.lambda_17 = 0.528
	rawdata.R17_VSMOW = 0.00038475
	rawdata.Nominal_D47 = {
		'ETH-1': 0.258,
		'ETH-2': 0.256,
		'ETH-3': 0.691,
		}
	rawdata.crunch()

	rawdata.standardize(method = 'pooled')
	for s,x,sx in [
		('ETH-4',   0.51410, 0.00337),
		('IAEA-C1', 0.36024, 0.00274),
		('IAEA-C2', 0.72459, 0.00341),
		('MERCK',   0.56752, 0.00573),
		]:
		assert(abs(rawdata.samples[s]['D47'] - x) < 1e-4)
		assert(abs(rawdata.samples[s]['SE_D47'] - sx) < 1e-4)

	rawdata.standardize(method = 'indep_sessions')
	for s,x,sx in [
		('ETH-4',   0.51420, 0.00337),
		('IAEA-C1', 0.36027, 0.00273),
		('IAEA-C2', 0.72444, 0.00340),
		('MERCK',   0.56719, 0.00572),
		]:
		assert(abs(rawdata.samples[s]['D47'] - x) < 1e-4)
		assert(abs(rawdata.samples[s]['SE_D47'] - sx) < 1e-4)


def test_virtual_data():
	args = dict(
		samples = [
			dict(Sample = 'ETH-1', N = 2),
			dict(Sample = 'ETH-2', N = 2),
			dict(Sample = 'ETH-3', N = 2),
			dict(Sample = 'FOO', N = 2,
				d13C_VPDB = -5., d18O_VPDB = -10.,
				D47 = 0.3, D48 = 0.15),
			], rD47 = 0.010, rD48 = 0.030)
	session1 = D47crunch.virtual_data(session = 'Session_01', **args, seed = 123)
	session2 = D47crunch.virtual_data(session = 'Session_02', **args, seed = 1234)
	
	out = str(session1 + session2).replace('}, {', '},\n{')
	expected = """[{'Sample': 'ETH-1', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 6.018962419832796, 'd46': 10.747025976136415, 'd47': 16.121901094852547, 'd48': 21.279952502419693, 'd49': 27.780042340821876, 'Session': 'Session_01'},
{'Sample': 'ETH-1', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 6.018962419832796, 'd46': 10.747025976136415, 'd47': 16.129740828798198, 'd48': 21.279754826864856, 'd49': 27.780042340821876, 'Session': 'Session_01'},
{'Sample': 'ETH-2', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -5.995859034201412, 'd46': -5.976075809403958, 'd47': -12.691522035559018, 'd48': -12.221351486037499, 'd49': -18.023380693805603, 'Session': 'Session_01'},
{'Sample': 'ETH-2', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -5.995859034201412, 'd46': -5.976075809403958, 'd47': -12.705325036691246, 'd48': -12.276730377250873, 'd49': -18.023380693805603, 'Session': 'Session_01'},
{'Sample': 'ETH-3', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 5.7423735352262195, 'd46': 11.16126986237509, 'd47': 16.682044554323745, 'd48': 22.303697873335423, 'd49': 28.306614367262117, 'Session': 'Session_01'},
{'Sample': 'ETH-3', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 5.7423735352262195, 'd46': 11.16126986237509, 'd47': 16.67771512407607, 'd48': 22.240124967341018, 'd49': 28.306614367262117, 'Session': 'Session_01'},
{'Sample': 'FOO', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -0.8404134452345557, 'd46': 2.828738359100358, 'd47': 1.3024427057226597, 'd48': 5.397915952831229, 'd49': 4.665655016163894, 'Session': 'Session_01'},
{'Sample': 'FOO', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -0.8404134452345557, 'd46': 2.828738359100358, 'd47': 1.3173114502005656, 'd48': 5.3684386774060195, 'd49': 4.665655016163894, 'Session': 'Session_01'},
{'Sample': 'ETH-1', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 6.018962419832796, 'd46': 10.747025976136415, 'd47': 16.12328304126664, 'd48': 21.248715348216987, 'd49': 27.780042340821876, 'Session': 'Session_02'},
{'Sample': 'ETH-1', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 6.018962419832796, 'd46': 10.747025976136415, 'd47': 16.134824966387452, 'd48': 21.299427708214044, 'd49': 27.780042340821876, 'Session': 'Session_02'},
{'Sample': 'ETH-2', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -5.995859034201412, 'd46': -5.976075809403958, 'd47': -12.70264564219742, 'd48': -12.237596504679146, 'd49': -18.023380693805603, 'Session': 'Session_02'},
{'Sample': 'ETH-2', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -5.995859034201412, 'd46': -5.976075809403958, 'd47': -12.706716415672554, 'd48': -12.191266291518088, 'd49': -18.023380693805603, 'Session': 'Session_02'},
{'Sample': 'ETH-3', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 5.7423735352262195, 'd46': 11.16126986237509, 'd47': 16.67641048241847, 'd48': 22.241315859156344, 'd49': 28.306614367262117, 'Session': 'Session_02'},
{'Sample': 'ETH-3', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 5.7423735352262195, 'd46': 11.16126986237509, 'd47': 16.690591779288617, 'd48': 22.27612966080137, 'd49': 28.306614367262117, 'Session': 'Session_02'},
{'Sample': 'FOO', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -0.8404134452345557, 'd46': 2.828738359100358, 'd47': 1.3002400472076234, 'd48': 5.331865982808294, 'd49': 4.665655016163894, 'Session': 'Session_02'},
{'Sample': 'FOO', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -0.8404134452345557, 'd46': 2.828738359100358, 'd47': 1.3170158915471402, 'd48': 5.309309108617985, 'd49': 4.665655016163894, 'Session': 'Session_02'}]"""
	
	assert(out == expected)

	
def test_D47_D48():

	rawdata = [
		{'Sample': 'ETH-1', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 6.018962419832796, 'd46': 10.747025976136415, 'd47': 16.121901094852547, 'd48': 21.279952502419693, 'd49': 27.780042340821876, 'Session': 'Session_01'},
		{'Sample': 'ETH-1', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 6.018962419832796, 'd46': 10.747025976136415, 'd47': 16.129740828798198, 'd48': 21.279754826864856, 'd49': 27.780042340821876, 'Session': 'Session_01'},
		{'Sample': 'ETH-2', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -5.995859034201412, 'd46': -5.976075809403958, 'd47': -12.691522035559018, 'd48': -12.221351486037499, 'd49': -18.023380693805603, 'Session': 'Session_01'},
		{'Sample': 'ETH-2', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -5.995859034201412, 'd46': -5.976075809403958, 'd47': -12.705325036691246, 'd48': -12.276730377250873, 'd49': -18.023380693805603, 'Session': 'Session_01'},
		{'Sample': 'ETH-3', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 5.7423735352262195, 'd46': 11.16126986237509, 'd47': 16.682044554323745, 'd48': 22.303697873335423, 'd49': 28.306614367262117, 'Session': 'Session_01'},
		{'Sample': 'ETH-3', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 5.7423735352262195, 'd46': 11.16126986237509, 'd47': 16.67771512407607, 'd48': 22.240124967341018, 'd49': 28.306614367262117, 'Session': 'Session_01'},
		{'Sample': 'FOO',   'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -0.8404134452345557, 'd46': 2.828738359100358, 'd47': 1.3024427057226597, 'd48': 5.397915952831229, 'd49': 4.665655016163894, 'Session': 'Session_01'},
		{'Sample': 'FOO',   'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -0.8404134452345557, 'd46': 2.828738359100358, 'd47': 1.3173114502005656, 'd48': 5.3684386774060195, 'd49': 4.665655016163894, 'Session': 'Session_01'},
		{'Sample': 'ETH-1', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 6.018962419832796, 'd46': 10.747025976136415, 'd47': 16.12328304126664, 'd48': 21.248715348216987, 'd49': 27.780042340821876, 'Session': 'Session_02'},
		{'Sample': 'ETH-1', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 6.018962419832796, 'd46': 10.747025976136415, 'd47': 16.134824966387452, 'd48': 21.299427708214044, 'd49': 27.780042340821876, 'Session': 'Session_02'},
		{'Sample': 'ETH-2', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -5.995859034201412, 'd46': -5.976075809403958, 'd47': -12.70264564219742, 'd48': -12.237596504679146, 'd49': -18.023380693805603, 'Session': 'Session_02'},
		{'Sample': 'ETH-2', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -5.995859034201412, 'd46': -5.976075809403958, 'd47': -12.706716415672554, 'd48': -12.191266291518088, 'd49': -18.023380693805603, 'Session': 'Session_02'},
		{'Sample': 'ETH-3', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 5.7423735352262195, 'd46': 11.16126986237509, 'd47': 16.67641048241847, 'd48': 22.241315859156344, 'd49': 28.306614367262117, 'Session': 'Session_02'},
		{'Sample': 'ETH-3', 'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': 5.7423735352262195, 'd46': 11.16126986237509, 'd47': 16.690591779288617, 'd48': 22.27612966080137, 'd49': 28.306614367262117, 'Session': 'Session_02'},
		{'Sample': 'FOO',   'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -0.8404134452345557, 'd46': 2.828738359100358, 'd47': 1.3002400472076234, 'd48': 5.331865982808294, 'd49': 4.665655016163894, 'Session': 'Session_02'},
		{'Sample': 'FOO',   'D17O': 0.0, 'd13Cwg_VPDB': -4.0, 'd18Owg_VSMOW': 26.0, 'd45': -0.8404134452345557, 'd46': 2.828738359100358, 'd47': 1.3170158915471402, 'd48': 5.309309108617985, 'd49': 4.665655016163894, 'Session': 'Session_02'}
		]

	data47 = D47crunch.D47data(rawdata)
	data47.crunch()
	data47.standardize()

	data48 =D47crunch.D48data(rawdata)
	data48.crunch()
	data48.standardize()

	table = D47crunch.table_of_samples(data47, data48, output = 'pretty')
	expected = """
––––––  –  –––––––––  ––––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––
Sample  N  d13C_VPDB  d18O_VSMOW     D47      SE    95% CL      SD  p_Levene     D48      SE    95% CL      SD  p_Levene
––––––  –  –––––––––  ––––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––
ETH-1   4       2.02       37.02  0.2052                    0.0054            0.1380                    0.0214          
ETH-2   4     -10.17       19.88  0.2085                    0.0057            0.1380                    0.0300          
ETH-3   4       1.71       37.45  0.6132                    0.0057            0.2700                    0.0268          
FOO     4      -5.00       28.91  0.2949  0.0044  ± 0.0100  0.0088     0.058  0.1457  0.0189  ± 0.0427  0.0304     0.853
––––––  –  –––––––––  ––––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––  ––––––  ––––––  ––––––––  ––––––  ––––––––
"""[1:]
	
	assert(table == expected)

if __name__ == '__main__':

# 	test_correlated_sum()
# 	test_fCO2eq()
# 	test_make_csv()
# 	test_pf()
# 	test_smart_type()
# 	test_pretty_table()
# 	test_w_avg()
# 	test_transpose_table()
# 	test_D47data_init()
# 	test_D47data_standardize()
# 	test_virtual_data()
# 	test_D47_D48()
	pass