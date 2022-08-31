#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Oct 2021
# Description: convert a VCF file (VCF 4.0) to the Structure input format (and optionally fastStructure too)
##########################


import sys
import os
import argparse


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to convert a VCF file to the input format required by Structure;' +
												' it is up to the user to ensure the VCF has the SNPs of interest (e.g. already one per locus)')


# add arguments to parse
parser.add_argument('vcf_file', type = str, help = 'The VCF file to convert')
parser.add_argument('-p', type = str, dest = 'pops_file', help = 'A pops file of sample identifiers (same as in the VCF) and numerical ' +
																'population designation, one per line')
parser.add_argument('-f', type = str, dest = 'fast', help = 'Whether to output a fastStructure format file, yes or no [default]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file
pops_file = args.pops_file
fast = args.fast

mode = 'Structure'

if fast:
	if fast.lower() == 'yes':
		output_fast = True
		mode = 'FastStructure'
else:
	output_fast = False



# process the pops file
if pops_file:
	samples = []
	pop_dict = {}
	with open(pops_file, 'r') as input_file:
		for line in input_file:
			elements = line.rstrip().split()
			samples.append(elements[0])
			pop_dict[elements[0]] = elements[1]
else:
	parser.print_help(sys.stderr)
	sys.exit(1)


# process the VCF to grab sample labels and genotypes for each SNP locus
with open(vcf_file, 'r') as vcf:
	snps = []
	count_snps = 0
	for line in vcf:
		if line.startswith('#'):        # a header INFO line
			if line.startswith('#CHROM'):		# the line with sample names
				sample_labels = line.rstrip().split()[9:]
		else:
			calls = line.rstrip().split()[9:]
			genotypes = []
			for call in calls:
				'''	The format is GT:DP:CATG, so 0/0:85:0,85,0,0 for a homozygous AA with 85 depth
					We want to grab the genotype (GT) 0/0 then split that into the two alleles 0 and 0 as a list
				'''
				genotypes.append(call.split(':')[0].split('/'))
			snps.append(genotypes)
			count_snps = count_snps + 1
	print('Read in a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) + ' SNP loci')


# only retain non-monomorphic SNP loci
keep_indices = []
index = 0
for genotypes in snps:
	gts = []
	mono = True
	for gt in genotypes:
		if gt == ['.', '.']:		# ignore missing
			continue
		else:
			if gt not in gts:
				gts.append(gt)
		if len(gts) > 1:
			mono = False
			break
	if not mono:
		keep_indices.append(index)
	index = index + 1
new_snps = [snps[i] for i in keep_indices]
snps = new_snps
print('Retained ' + str(len(snps)) + ' loci after filtering for monomorphic')


# create the output lines as ordered in the pops file
lines_out = []
for sample in samples:
	if sample in sample_labels:
		this_index = sample_labels.index(sample)
		gt_a = []
		gt_b = []
		for genotypes in snps:
			gt = genotypes[this_index]
			if '.' in gt:	# missing data
				gt_a.append('-9')
				gt_b.append('-9')
			else:
				gt_a.append(str(int(gt[0]) + 1))
				gt_b.append(str(int(gt[1]) + 1))
		if output_fast:
			line_out_a = 5 * (str(pop_dict[sample]) + '\t') + sample + '\t' + '\t'.join(gt_a)
			line_out_b = 5 * (str(pop_dict[sample]) + '\t') + sample + '\t' + '\t'.join(gt_b)
		else:
			line_out_a = sample + '\t' + str(pop_dict[sample]) + '\t' + '\t'.join(gt_a)
			line_out_b = sample + '\t' + str(pop_dict[sample]) + '\t' + '\t'.join(gt_b)
		lines_out.append(line_out_a)
		lines_out.append(line_out_b)
	else:
		print('Sample ' + sample + ' was not present in the VCF but was present in your samples file')


# create the Structure format file
with open(os.path.basename(vcf_file) + '.str', 'w') as outfile:
	for line_out in lines_out:
		outfile.write(line_out + '\n')


# report completion
print('Converted a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) +
	' SNP loci to ' + mode + ' format with ' + str(int(len(lines_out)/2)) +	' samples and ' +
	str(len(snps)) + ' SNP loci')
