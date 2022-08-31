#!/usr/bin/env python

##########################
# Author: B. Anderson, based on https://github.com/ksil91/Ostrea_PopStructure/blob/master/Scripts/subsetSNPs.py
# Date: Sep 2021
# Description: Select a single SNP per locus (highest sample coverage) from a VCF file produced by ipyrad (VCF 4.0)
##########################


import sys			# allows access to command line arguments
import argparse
import random


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to select a single SNP per locus (highest sample coverage) from a VCF file produced by ipyrad and ' +
                                                'output a VCF with the same name but prefixed with "mod_"')


# add arguments to parse
parser.add_argument('vcf_file', type = str, help = 'The VCF file to select SNPs from')
parser.add_argument('-r', type = str, dest = 'sample_random', help = 'Whether to sample a random SNP in a locus (rather than first) in the case of a tie ' +
                                                                'in coverage ("yes" or "no" [default])')
parser.add_argument('-d', type = int, dest = 'rand', help = 'Randomly downsample to this many SNPs')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file
sample_random = args.sample_random
random_down = args.rand

if not sample_random:
	sample_random = False
elif sample_random.lower() == 'yes':
    sample_random = True
else:
    sample_random = False


# filter the VCF
with open(vcf_file, 'r') as vcf, open('mod_' + vcf_file, 'w') as outvcf:
	this_locus = ''
	count_loci = 0
	count_snps = 0
	lines = []
	NS_tally = []
	outlines = []
	for line in vcf:
		if line.startswith('#'):        # a header INFO line
			outvcf.write(line)
		else:
			line_fields = line.rstrip().split()
			locus = line_fields[0]
			snp_NS = line_fields[7].split(';')[0].split('=')[1]           # column 8 has NS=xx;DP=xx
			if locus != this_locus:
				if count_loci == 0:     # the first locus
					this_locus = locus
					lines.append(line)
					NS_tally.append(int(snp_NS))
					count_loci = count_loci + 1
					count_snps = count_snps + 1
				else:                   # we have finished iterating through SNPs for the last locus
					# select from the lines that have been stored
					if sample_random:   # pick max NS SNP (randomly if there are ties)
						indices = [index for index, val in enumerate(NS_tally) if val == max(NS_tally)]
						outlines.append(lines[random.choice(indices)])
						#outvcf.write(lines[random.choice(indices)])
					else:               # pick max NS SNP (first if there are ties)
						outlines.append(lines[NS_tally.index(max(NS_tally))])
						#outvcf.write(lines[NS_tally.index(max(NS_tally))])
					# reset lists and start tallying next locus
					this_locus = locus
					lines = []
					NS_tally = []
					lines.append(line)
					NS_tally.append(int(snp_NS))
					count_loci = count_loci + 1
					count_snps = count_snps + 1
			else:           # we are still iterating through SNPs for the locus
				lines.append(line)
				NS_tally.append(int(snp_NS))
				count_snps = count_snps + 1
	# after going through all the loci, need to select from the last one
	if sample_random:   # pick max NS SNP (randomly if there are ties)
		indices = [index for index, val in enumerate(NS_tally) if val == max(NS_tally)]
		outlines.append(lines[random.choice(indices)])
		#outvcf.write(lines[random.choice(indices)])
	else:               # pick max NS SNP (first if there are ties)
		outlines.append(lines[NS_tally.index(max(NS_tally))])
		#outvcf.write(lines[NS_tally.index(max(NS_tally))])
	# randomly downsample the SNPs to the desired number if requested
	if random_down:
		if len(outlines) > random_down:
			new_lines = random.sample(outlines, random_down)
			outlines = new_lines
			print('Randomly downsampled ' + str(count_loci) + ' single SNPs to ' + str(random_down) + ' SNPs')
			count_loci = random_down
	# write the output lines
	for line in outlines:
		outvcf.write(line)


print('Wrote a new VCF with ' + str(count_loci) + ' SNPs sampled from ' + str(count_snps) + ' SNPs')
if sample_random:
    print('In the case of ties for sample coverage in a locus, SNPs were selected randomly')
else:
    print('In the case of ties for sample coverage in a locus, the first SNP was selected')
