#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Jan 2023
# Description: collect groups of clones from two columns of samples with distances less than a threshold
##########################


import sys
import argparse


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to collate clones from a comparison table with a given threshold of difference')


# add arguments to parse
parser.add_argument(type = str, dest = 'comp_file', help = 'The comparisons file with three columns (no header): ' + \
	'individual 1 (string), individual 2 (string), distance (float)')
parser.add_argument('-t', type = float, dest = 'thresh', help = 'The threshold of difference below which to consider the comparison clonal (default: 0.02)')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
comp_file = args.comp_file
thresh = args.thresh

if not comp_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not thresh:
	thresh = 0.02


# read in the file, line by line, processing and keeping comparisons less than the threshold
comparisons = []
with open(comp_file, 'r') as infile:
	for line in infile:
		parts = line.strip().split()
		comp = parts[0:2]
		dist = float(parts[2])
		if (dist < thresh):
			comparisons.append(comp)


# create sets iteratively by looping through the sets for each comparison
sets = []
sets.append(comparisons[0])
for comp in comparisons:
	found = False
	for thisset in sets:
		if comp[0] in thisset:
			if comp[1] in thisset:
				found = True
				break
			else:
				thisset.append(comp[1])
				found = True
				break
		elif comp[1] in thisset:
			thisset.append(comp[0])
			found = True
			break
		else:
			continue
	if not found:
		sets.append(comp)


# iterate through the sets to merge any overlaps
outsets = []
outsets.append(sets[0])
options = set([i for i in range(len(sets))])
merged_indices = []
merged_indices.append(0)
while len(options) > 0:
	for index,thisset in enumerate(sets):
		if index in merged_indices:
			continue
		else:
			end_loop = False
			while not end_loop:
				for thatset in outsets:
					for item in thisset:
						if item in thatset:		# overlap
							thatset.extend(thisset)
							merged_indices.append(index)
							end_loop = True
					else:
						continue
				end_loop = True
	options = options.difference(set(merged_indices))
	option_list = list(options)
	if len(option_list) > 0:
		outsets.append(sets[option_list[0]])
		merged_indices.append(option_list[0])


# Now report the sets
outlist = []
for thisset in outsets:
	outlist.append('\n'.join(set(thisset)))

print('\n\n'.join(outlist))
