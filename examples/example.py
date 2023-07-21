#! /usr/bin/env python3
'''
Example of using D47crunch in a script.
'''

import D47crunch

# create a D47data() object
clumpy = D47crunch.D47data(verbose = True)

clumpy.Nominal_D47 = {
	k: clumpy.Nominal_D47[k]
	for k in ['ETH-1', 'ETH-2', 'ETH-3']
	}

# read raw data from external csv file
clumpy.read('rawdata.csv')

# compute WG compositions for each session
clumpy.wg()

# compute δ13C, δ18O, and raw Δ47, Δ48, Δ49 values
clumpy.crunch()

# compute absolute Δ47 values
clumpy.standardize()

# print out a summary of the standardization results
clumpy.summary(print_out = True, save_to_file = False)

# print out some information about the sessions
clumpy.table_of_sessions(print_out = True, save_to_file = False)

# print out final δ13C, δ18O, and Δ47 values averaged by sample
clumpy.table_of_samples(print_out = True, save_to_file = False)

# print out all analyses after processing
clumpy.table_of_analyses(print_out = True, save_to_file = False)

# save the tables above as csv files
clumpy.summary(print_out = False, save_to_file = True, dir = 'csv', filename = 'summary.csv')
clumpy.table_of_sessions(print_out = False, save_to_file = True, dir = 'csv', filename = 'sessions_table.csv')
clumpy.table_of_samples(print_out = False, save_to_file = True, dir = 'csv', filename = 'samples_table.csv')
clumpy.table_of_analyses(print_out = False, save_to_file = True, dir = 'csv', filename = 'analyses_table.csv')

# generate plots for each session
clumpy.plot_sessions(dir = 'session_plots')

# save D47 data with correlations
clumpy.save_D47_correl()