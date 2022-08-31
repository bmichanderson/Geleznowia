#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Dec 2021
# Updated: May 2022
# Description: convert a VCF file (VCF 4.0) to a Nexus file for SVDquartets
##########################


import sys
import argparse
import os


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to convert a VCF file to a Nexus input format required by SVDquartets;' +
	' it is up to the user to ensure the VCF has the SNPs of interest (e.g. one per locus)')


# add arguments to parse
parser.add_argument('-v', type = str, dest = 'vcf_file', help = 'The VCF file to convert')
parser.add_argument('-s', type = str, dest = 'sample_file',
	help = 'A tab-delimited file of sample name and taxon, one per line [required if specifying higher groups]')
parser.add_argument('-f', type = str, dest = 'out_format', help = 'Output format as two haplotypes ("hap"; default)' +
	' or single sequence ("seq") with ambiguities')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file
sample_file = args.sample_file
out_format = args.out_format

if not vcf_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not out_format or out_format.lower() == 'hap':
	out_format = 'hap'
	print('Output format is two haplotypes per sample')
elif out_format.lower() == 'seq':
	out_format = 'seq'
	print('Output format is a single sequence per sample')
else:
	print('Incorrectly specified format type')
	parser.print_help(sys.stderr)
	sys.exit(1)


# create an ambiguity dictionary
amb_dict = {
	'CT': 'Y',
	'TC': 'Y',
	'AG': 'R',
	'GA': 'R',
	'AT': 'W',
	'TA': 'W',
	'CG': 'S',
	'GC': 'S',
	'GT': 'K',
	'TG': 'K',
	'AC': 'M',
	'CA': 'M'
}


# process the VCF to grab sample labels and genotypes for each SNP locus
with open(vcf_file, 'r') as vcf:
	snps = []
	count_snps = 0
	for line in vcf:
		if line.startswith('#'):        # a header INFO line
			if line.startswith('#CHROM'):		# the line with sample names
				sample_labels = line.rstrip().split()[9:]
		else:
			parts = line.rstrip().split()
			all_bases = []
			all_bases.append(parts[3])	# the ref allele
			if ',' in parts[4]:		# a multi-allele
				for allele in parts[4].split(','):
					all_bases.append(allele)
			else:
				all_bases.append(parts[4])
			calls = parts[9:]
			genotypes = []
			for call in calls:
				'''
					The format is GT:DP:CATG, so 0/0:85:0,85,0,0 for a homozygous AA with 85 depth
					We want to grab the genotype (GT) 0/0 then split that into the two alleles 0 and 0 as a list
					Then convert that into bases
				'''
				genotype = call.split(':')[0].split('/')
				allele1 = genotype[0]
				allele2 = genotype[1]
				if allele1 == '.':		# missing data
					genotypes.append(['?', '?'])
				else:
					genotypes.append([all_bases[int(allele1)], all_bases[int(allele2)]])
			snps.append(genotypes)
			count_snps = count_snps + 1


# process the sample file (if present) into a dictionary
# e.g., taxon1: [sample1, sample2], taxon2: [sample3], taxon3: [...
if sample_file:
	taxon_dict = {}
	taxa = []
	samples = []
	with open(sample_file, 'r') as infile:
		for line in infile:
			parts = line.strip().split()
			samples.append(parts[0])
			if parts[1] not in taxa:
				taxa.append(parts[1])
				taxon_dict[parts[1]] = [parts[0]]
			else:
				taxon_dict[parts[1]].append(parts[0])


# determine if any samples in the VCF were not in the sample list
samples_missed = []
if sample_file:
	for sample in sample_labels:
		if sample not in samples:
			samples_missed.append(sample)

	if len(samples_missed) > 0:
		print('\nSamples present in the VCF but not the sample list (so not included):')
		print(samples_missed)


# determine if any samples in the sample list weren't in the VCF, and 
# whether any taxa should be excluded if all their samples weren't
samples_excl = []
taxa_excl = []
if sample_file:
	for sample in samples:
		if sample not in sample_labels:
			samples_excl.append(sample)
		
	if len(samples_excl) > 0:
		print('\nSamples present in the sample list but not the VCF (so not included):')
		print(samples_excl)

	for taxon in taxa:
		keep_samples = []
		keep_taxon = False
		for sample in taxon_dict[taxon]:
			if sample in sample_labels:
				keep_taxon = True
				keep_samples.append(sample)
		
		if not keep_taxon:
			taxa_excl.append(taxon)
		else:		# update dictionary
			taxon_dict.update({taxon: keep_samples})

	if len(taxa_excl) > 0:
		print('\nTaxa lacking samples in the VCF (so not included):')
		print(taxa_excl)


# create the output samples list
samples_out = [item for item in sample_labels if item not in samples_missed]


# create the output taxa if present
if sample_file:
	taxa_out_dict = {taxon: taxon_dict[taxon] for taxon in taxon_dict if taxon not in taxa_excl}


# create the data lines to write to the Nexus file
lines_out = []
for sample in samples_out:
	for index, sample_label in enumerate(sample_labels):
		if sample_label == sample:
			this_index = index
			break
	bases1 = []
	bases2 = []
	sequence = []
	# sample names can't have dashes in PAUP, so substitute with periods
	sample = sample.replace('-', '.')
	for genotypes in snps:
		genotype = genotypes[this_index]
		bases1.append(genotype[0])
		bases2.append(genotype[1])
		if genotype[0] == genotype[1]:
			sequence.append(genotype[0])
		else:
			sequence.append(amb_dict[genotype[0] + genotype[1]])
	hap1_line = sample + '_A' + '\t' + ''.join(bases1)
	hap2_line = sample + '_B' + '\t' + ''.join(bases2)
	seq_line = sample + '\t' + ''.join(sequence)
	if out_format == 'seq':
		lines_out.append(seq_line)
	else:
		lines_out.append(hap1_line)
		lines_out.append(hap2_line)


# create the sets lines to write (if sample_file present)
if sample_file:
	sets_lines = []
	for taxon in taxa_out_dict:
		if out_format == 'seq':
			out_line = taxon + ': ' + ' '.join(sample for sample in taxa_out_dict[taxon])
		else:
			out_line = taxon + ': ' + ' '.join(sample + '_A ' + sample + '_B' for sample in taxa_out_dict[taxon])
		out_line = out_line.replace('-', '.')		# to avoid error in PAUP
		sets_lines.append(out_line)


# create the Nexus file
num_taxa = len(lines_out)
num_chars = len(snps)
with open(os.path.basename(vcf_file).replace('.vcf', '') + '.nex', 'w') as outfile:
	outfile.write('#Nexus\n\n')
	data_block = ('begin data;\n\tdimensions ntax=' + str(num_taxa) + ' nchar=' + str(num_chars) + ';\n\t' +
				'format datatype=nucleotide missing=?;\nmatrix\n')
	outfile.write(data_block)
	for line_out in lines_out:
		outfile.write(line_out + '\n')
	outfile.write('\t;\nend;\n')
	if sample_file:
		sets_block = ('begin sets;\ntaxpartition svdtaxa =\n\t' + ',\n\t'.join(sets_lines) + '\n\t;\nend;\n')
		outfile.write(sets_block)


# report completion
print('\nConverted a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) +
	' SNPs to Nexus format with ' +	str(len(lines_out)) + ' taxa (lineages) and ' + str(num_chars) + ' characters')
if sample_file:
	print('Created a sets block with ' + str(len(sets_lines)) + ' higher taxa designated as taxpartition svdtaxa')
