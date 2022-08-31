#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Oct 2021
# Modified: Aug 2022 (changed font and output file format)
# Description: compute and visualize likelihood values from multiple structure-like runs to assess best K
##########################


import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# set the font output
matplotlib.rcParams['font.family'] = 'Arial'


# to avoid warnings when creating ragged arrays in np
# (e.g. when the number of reps per K is different)
np.warnings.filterwarnings('ignore', category = np.VisibleDeprecationWarning)


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to compute and visualize likelihoods from Structure-like runs for assessing best K')


# add arguments to parse
parser.add_argument('-l', type = str, dest = 'like_file', help = 'The likelihoods file for the K reps (multiple lines per K) in tab delimited form:' + 
																	'K_value likelihood')
parser.add_argument('-o', type = str, dest = 'output', help = 'Name of the output SVG without extension [default \"output\"]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
like_file = args.like_file
out_pre = args.output

if not like_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not out_pre:
	out_pre = 'output'


# load the likelihoods file as a dataframe
like_df = pd.read_csv(like_file, sep = '\t', header = None)
like_df.rename(columns = {	0: 'K',
							1: 'Likelihood'},
							inplace = True)


# deltaK
# Likelihood mean = L(K)
# L'(K) = L(K) - L(K - 1)
# |L''(K)| = |L'(K + 1) - L'(K)|
# dK	= mean [|L''(K)|] / std [L(K)]
# dK	= mean [|L(K + 1) - 2 L(K) + L(K - 1)|] / std [L(K)]

# compute mean and standard deviation, and values needed for deltaK
like_comps = like_df.groupby('K').agg([np.mean, np.std])
like_comps = like_comps['Likelihood']
like_comps.loc[:, 'Lprime'] = like_comps['mean'] - like_comps['mean'].shift(1)
like_comps.loc[:, 'absLdprime'] = np.abs(like_comps['Lprime'].shift(-1) - like_comps['Lprime'])
like_comps.loc[:, 'deltaK'] = like_comps['absLdprime'] / like_comps['std']


# plot visualizations of log probabilities and deltaK
fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, figsize = (7, 9))
plt.subplots_adjust(hspace = 0.2)
plt.rcParams['font.size'] = 20
like_comps.plot(y = 'mean', marker = '.', markersize = 12,
					title = 'Log Prob of Data L(K)', legend = False,
					linewidth = 0.75, color = 'steelblue', ax = ax1)
ax1.grid()
ax1.errorbar(x = like_comps.index, y = like_comps['mean'], yerr = like_comps['std'],
			fmt = 'none', color = 'black', capsize = 4, elinewidth = 1, zorder = 3)
ax1.get_figure().suptitle('')
ax1.set_title('Log Prob of Data L(K)')
ax1.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))
ax1.yaxis.get_offset_text().set_fontsize(18)
ax1.yaxis.set_tick_params(labelsize = 18)
like_comps.plot(y = 'deltaK', marker = '.', markersize = 12,
					title = '\u0394K', legend = False,
					color = 'steelblue', ax = ax2)
ax2.set_xlabel('K', fontsize = 22)
ax2.set_xticks(like_comps.index)
ax2.xaxis.set_tick_params(labelsize = 18)
ax2.yaxis.set_tick_params(labelsize = 18)
plt.grid()
plt.savefig(out_pre + '.svg')
