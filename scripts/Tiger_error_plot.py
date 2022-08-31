#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Nov 2021
# Description: provide summary plots from a Tiger genotyping error rates text file
##########################


import sys
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to summarise error rates calculated by Tiger')


# add arguments to parse
parser.add_argument('tiger_file', type = str, help = 'The Tiger error rates text file')
parser.add_argument('-o', type = str, dest = 'output', help = 'Name of the output pdf without extension [default \"output\"]')
parser.add_argument('-m', type = str, dest = 'maxd', help = 'Maximum potential read depth to graph [default 500]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
tiger_file = args.tiger_file
out_pre = args.output
maxd = int(args.maxd)

if not out_pre:
	out_pre = 'output'
if not maxd:
	maxd = 500


# load the Tiger file as a dataframe
error_df = pd.read_csv(tiger_file, sep = '\t')


# plot visualizations of depth and error rates
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex = True,
										figsize = (12, 12))
plt.subplots_adjust(hspace = 0.2)
## depth
obs_maxd = max(error_df['maxDepth'])
error_df.plot('maxDepth', 'Size', xlim = (0, min(maxd, obs_maxd)),
				title = 'Distribution of read depth', legend = None,
				ylabel = 'Sites', ax = ax1)
## overall error rate by depth
error_df.plot('maxDepth', 'ErrorRate', xlim = (0, min(maxd, obs_maxd)),
				title = 'Estimated overall error rate', legend = None,
				ylabel = 'Error rate', ylim = (-0.02, 0.502), ax = ax2)
## homo error rate by depth
error_df.plot('maxDepth', 'ErrorRateHom', xlim = (0, min(maxd, obs_maxd)),
				title = 'Estimated homozygous error rate', legend = None,
				ylabel = 'Error rate', ylim = (-0.02, 0.502), ax = ax3)
## het error rate by depth
error_df.plot('maxDepth', 'ErrorRateHet', xlim = (0, min(maxd, obs_maxd)),
				title = 'Estimated heterozygous error rate', legend = None,
				ylabel = 'Error rate', ylim = (-0.02, 0.502), ax = ax4)
ax4.set_xlabel('Depth')
plt.savefig(out_pre + '.pdf')
