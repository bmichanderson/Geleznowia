#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Dec 2021
# Updated: April 2022
# Description: filter a VCF file by samples and/or snps
##########################


import sys
import argparse
import statistics
import random


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to filter a VCF file')


# add arguments to parse
parser.add_argument('vcf_file', type = str, help = 'The VCF file to filter')
parser.add_argument('-o', type = str, dest = 'outpre', help = 'Output file prefix [default \"output\"]')
parser.add_argument('-s', type = str, dest = 'samps', help = 'A file with sample names to be dropped, one per line')
parser.add_argument('-p', type = str, dest = 'pops', help = 'A file with samples and population identifiers (for pop filtering), one per line and tab separated')
parser.add_argument('--minperpop', type = int, dest = 'minperpop', help = 'The minimum number of samples per population that must have a locus to keep it')
parser.add_argument('--minpops', type = int, dest = 'minpops', help = 'The minimum number of populations that must have a locus to keep it')
parser.add_argument('--inpop', type = str, dest = 'inpop', help = 'A file with a list of populations that must have a locus to keep it, one per line')
parser.add_argument('--bial', type = str, dest = 'bial', help = 'Whether a locus must be biallelic to be kept (\"yes\" or \"no\" [default])')
parser.add_argument('--mincov', type = float, dest = 'mincov', help = 'The minimum cover as a proportion of samples a locus must be found in to be kept')
parser.add_argument('--mac', type = int, dest = 'mac', help = 'The minimum minor allele count (how many times the minor allele must be called to keep the locus)')
parser.add_argument('--maxd', type = int, dest = 'maxd', help = 'The maximum call depth found anywhere in a locus to keep it')
parser.add_argument('--maxmd', type = int, dest = 'maxmd', help = 'The maximum mean call depth for a locus to keep it')
parser.add_argument('--mind', type = int, dest = 'mind', help = 'The minimum call depth found anywhere in a locus to keep it')
parser.add_argument('--minmd', type = int, dest = 'minmd', help = 'The minimum mean call depth for a locus to keep it')
parser.add_argument('-d', type = int, dest = 'rand', help = 'Randomly downsample to this many SNPs')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file
out_pre = args.outpre
sample_file = args.samps
pops_file = args.pops
minperpop = args.minperpop
minpops = args.minpops
inpop_file = args.inpop
bial = args.bial
mincov = args.mincov
minmac = args.mac
maxd = args.maxd
maxmd = args.maxmd
mind = args.mind
minmd = args.minmd
random_down = args.rand

if not vcf_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not out_pre:
	out_pre = 'output'


# parse the samples file, if present
if sample_file:
	samples_remove = []
	with open(sample_file, 'r') as samples:
		for sample in samples:
			samples_remove.append(sample.rstrip())


# load the pops file and create the pop dictionary
if pops_file:
	pops = []
	pop_dict = {}
	with open(pops_file, 'r') as sample_input:
		for line in sample_input:
			sampleID = line.rstrip().split()[0]
			pop = line.rstrip().split()[1]
			if pop not in pops:
				pops.append(pop)
				pop_dict[pop] = [sampleID]
			else:
				pop_dict[pop].append(sampleID)


# parse the inpop file, if present
if inpop_file:
	req_pops = []
	with open(inpop_file, 'r') as input:
		for pop in input:
			req_pops.append(pop.rstrip())


# parse and filter the VCF
with open(vcf_file, 'r') as vcf, open(out_pre + '.vcf', 'w') as outfile:
	in_snps = 0
	in_samples = 0
	out_snps = 0
	not_bial = 0
	outlines = []
	for line in vcf:
		if line.startswith('#'):		# a header INFO line
			if line.startswith('#CHROM'):		# the line with sample names
				pieces = line.rstrip().split()
				sample_labels = pieces[9:]
				in_samples = len(sample_labels)
				if sample_file:			# if we need to remove any
					indices = []
					for index, sample in enumerate(sample_labels):
						if sample not in samples_remove:
							indices.append(index)
					newline = '\t'.join(pieces[:9] + [pieces[9:][index] for index in indices])
					sample_labels = [sample_labels[index] for index in indices]
					outfile.write(newline + '\n')
				else:
					outfile.write(line)
				if pops_file:		# if we need to check pop membership
					pop_labels = []
					for sample_label in sample_labels:
						for pop in pops:
							if sample_label in pop_dict[pop]:
								pop_labels.append(pop)
			else:
				outfile.write(line)
		else:
			in_snps = in_snps + 1
			pieces = line.rstrip().split()
			if bial:
				if bial.lower() == 'yes':
					if len(pieces[4].split(',')) > 1:
						not_bial = not_bial + 1
						continue
			if sample_file:		# if samples were removed
				calls = [pieces[9:][index] for index in indices]
			else:
				calls = pieces[9:]
			depths = []
			genotypes = []
			for call in calls:
				'''	The format is GT:DP:CATG, so 0/0:85:0,85,0,0 for a homozygous AA with 85 depth
					We want to grab the genotype (GT) and depth (DP)
				'''
				gt = call.split(':')[0].split('/')		# a list, e.g. ['0', '0']
				genotypes.append(gt)		
				if gt != ['.', '.']:		# if not an N (shouldn't be needed, but currently is)
					depths.append(int(call.split(':')[1]))
			# calculate stats depending on what filters need to be applied
			keep_locus = True
			while keep_locus:
				if inpop_file and pops_file:
					for pop in req_pops:
						pop_genos = [genotypes[ind] for ind, item in enumerate(pop_labels) if item == pop]
						count = len([geno for geno in pop_genos if geno != ['.', '.']])
						if count == 0:
							keep_locus = False
							break
					if not keep_locus:
						break
				if minperpop and pops_file:
					for pop in set(pop_labels):
						pop_genos = [genotypes[ind] for ind, item in enumerate(pop_labels) if item == pop]
						count = len([geno for geno in pop_genos if geno != ['.', '.']])
						if count < minperpop:
							keep_locus = False
							break
					if not keep_locus:
						break
				if minpops and pops_file:
					pops_have = 0
					for pop in set(pop_labels):
						pop_genos = [genotypes[ind] for ind, item in enumerate(pop_labels) if item == pop]
						count = len([geno for geno in pop_genos if geno != ['.', '.']])
						if count > 0:
							pops_have = pops_have + 1
						if pops_have >= minpops:
							break
					if pops_have < minpops:
						keep_locus = False
						break
				if mincov:
					count_present = len([geno for geno in genotypes if geno != ['.', '.']])
					if count_present / len(sample_labels) < mincov:
						keep_locus = False
						break
				if minmac:
					count_0 = ''.join(''.join(geno) for geno in genotypes).count('0')
					count_1 = ''.join(''.join(geno) for geno in genotypes).count('1')
					if any([count_0 < minmac, count_1 < minmac]):
						keep_locus = False
						break
				if maxd:
					if max(depths) > maxd:
						keep_locus = False
						break
				if maxmd:
					if sum(depths) > 0:
						if statistics.mean([depth for depth in depths if depth > 0]) > maxmd:
							keep_locus = False
							break
					else:
						keep_locus = False
						break
				if mind:
					if min(depths) < mind:
						keep_locus = False
						break
				if minmd:
					if sum(depths) > 0:
						if statistics.mean([depth for depth in depths if depth > 0]) < minmd:
							keep_locus = False
							break
					else:
						keep_locus = False
						break
				# retain the locus if it should be kept
				if keep_locus:
					newline = '\t'.join(pieces[:9] + calls)
					outlines.append(newline + '\n')
					out_snps = out_snps + 1
					break
	# randomly downsample the SNPs to the desired number if requested
	if random_down:
		if len(outlines) > random_down:
			new_lines = random.sample(outlines, random_down)
			outlines = new_lines
			out_snps = random_down
			print('Randomly downsampled filtered SNPs to ' + str(random_down) + ' SNPs')
	# report filtering
	if bial:
		if bial.lower() == 'yes':
			print('Dropped ' + str(not_bial) + ' SNPs that were not biallelic')
	# write the output lines
	for line in outlines:
		outfile.write(line)
	print('Read a VCF file with ' + str(in_samples) + ' samples and ' + str(in_snps) + ' SNPs' +
		' and filtered it to ' + str(len(sample_labels)) + ' samples and ' + str(out_snps) + ' SNPs')
