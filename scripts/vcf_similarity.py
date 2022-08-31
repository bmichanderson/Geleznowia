#!/usr/bin/env python3

##########################
# Author: B. Anderson
# Date: May 2022
# Description: assess the similarity between samples in a VCF file
##########################


import sys
import argparse
import copy
import numpy
import pandas
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# set up the parser and arguments
parser = argparse.ArgumentParser(description = 'A script to assess similarity between samples in a VCF file')
parser.add_argument('-v', type = str, dest = 'vcf_file', help = 'The VCF file to use')
parser.add_argument('-o', type = str, dest = 'out_pre', help = 'Prefix of the output files [default \"output\"]')


# parse the command line
if len(sys.argv[1:]) == 0:
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
vcf_file = args.vcf_file
out_pre = args.out_pre

if not vcf_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not out_pre:
	out_pre = 'output'


# read in the VCF file and parse it
geno_list = []
with open(vcf_file, 'r') as infile:
	for line in infile:
		if line.startswith('#'):	# the first block of the VCF
			if line.startswith('#CHROM'):	# the line with column names
				parts = line.rstrip().split()
				samples = parts[9: ]
				for sample in samples:
					geno_list.append([])
		else:	# the data lines
			parts = line.rstrip().split()
			calls = parts[9: ]
			index = 0
			for call in calls:
				genotype = call.split(':')[0].split('/')
				if genotype == ['.', '.']:
					geno_list[index].append('NA')
				else:
					geno_list[index].append(''.join(genotype))
				index = index + 1


# turn the list into an array (sample genotypes as rows)
myarray = numpy.array(geno_list)
print('Created an array of ' + str(myarray.shape[0]) + ' samples and ' +
	str(myarray.shape[1]) + ' SNPs\n')


# check all the possible genotypes
sets = []
for element in myarray:
	sets.append(set(element))

print('Genotypes: ' + str(sets[0].union(*sets[1:])) + '\n')


# initialise a list to store the distances
# each list element will be a list representing a row in the distance matrix
# e.g. 	[(1v1), (1v2), (1v3), (1v4),..]
#		[2v1, (2v2), (2v3), (2v4),..]
# where () will be missing/NA
mydists = []
for iteration in range(len(samples)):
	mydists.append([])


# initialise another list to store each pairwise comparison 
mycomps = []
comparisons = int(len(samples) * (len(samples) - 1) / 2)
for comp in range(comparisons):
	mycomps.append([])


# compare rows (samples) pairwise
print('Assessing ' + str(len(samples)) + ' samples\n')
print('Comparing sample ', end = '', flush = True)
comp_index = 0
for index in range(len(myarray) - 1):
	print(str(index + 1) + ' ', end = '', flush = True)
	samparray = myarray[index]
	for index2 in range(index + 1, len(myarray)):
		mydists[index].append(numpy.nan)
		otherarray = myarray[index2]
		con1 = samparray != 'NA'
		con2 = otherarray != 'NA'
		con3 = samparray != otherarray
		dist = sum((con1) & (con2) & (con3)) / sum((con1) & (con2))
		mydists[index2].append(dist)
		mycomps[comp_index].append(samples[index])
		mycomps[comp_index].append(samples[index2])
		mycomps[comp_index].append(dist)
		comp_index = comp_index + 1
print('\nDone comparing samples! Plotting...')

# add a last column with all NA
for index in range(len(myarray)):
	mydists[index].append(numpy.nan)


# create a dataframe
mydf = pandas.DataFrame(mydists)
mydf.columns = samples
mydf.index = samples


# write an output friendly copy of the distances to a file
printdf = copy.deepcopy(mydf)
printdf[numpy.isnan(printdf)] = 'NA'
printdf.to_csv(out_pre + '_dists.tab', sep = '\t')


# write the comparisons, sorted from largest to smallest
mycomps.sort(key = lambda x: x[2], reverse = True)
with open(out_pre + '_comps.txt', 'w') as outfile:
	for comp in mycomps:
		outfile.write('\t'.join(str(i) for i in comp) + '\n')


# plot a boxplot of a square dataframe and save to file
# see https://stackoverflow.com/questions/36250729/...
# how-to-convert-triangle-matrix-to-square-in-numpy
mydistarray = numpy.array(mydists)
mydistarray[numpy.isnan(mydistarray)] = 0.0
squarearray = numpy.triu(mydistarray.T, 1) + mydistarray
numpy.fill_diagonal(squarearray, numpy.nan)
plotdf = pandas.DataFrame(squarearray)
plotdf.columns = samples
plotdf.index = samples
myfig = plotdf.boxplot(rot = 90, figsize = (18, 12)).get_figure()
myfig.tight_layout()
myfig.savefig(out_pre + '_boxplot.png', dpi = 300)


# plot a histogram and save to file
values = [i[2] for i in mycomps]
myfig2, ax = plt.subplots(1, 1, figsize = (12, 8))
plt.hist(values, bins = 50)
start, end = ax.get_xlim()
ax.xaxis.set_ticks(numpy.arange(0, round(end, 2), 0.02))
plt.title('Histogram of Site Differences', fontsize = 18)
plt.ylabel('Count', fontsize = 14)
plt.xlabel('Proportion different sites', fontsize = 14)
ax.xaxis.set_tick_params(labelsize = 12)
ax.yaxis.set_tick_params(labelsize = 12)
myfig2.savefig(out_pre + '_hist.png', dpi = 300)
