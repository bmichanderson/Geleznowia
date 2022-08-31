#!/usr/bin/env python

# Launch this script from within the working folder, with number of cores as first arg, samples file as second

import sys
import ipyrad as ip
import ipyrad.analysis as ipa
import ipyparallel as ipp
from subprocess import call
import time


## Set parameter values for the run
project = "Geleznowia"

# steps 1 and 2 (already done, but have to repeat because I removed json files)
datatype = "ddrad"
restrict_overhang = ("TGCAG", "CCG")        # overhang after enzyme digestion; in this case PstI and MspI
barcode_mismatch = "1"                      # max number of mismatches in the barcode
min_trim_len = "50"                         # minimum read length to keep after trimming
max_low_qual = "2"                          # max number of low quality bases (Ns) in a read after trimming
phred = "43"                                # converts to min qscore = 30, for trimming at 3' end of reads
filter = "2"                                # level of filtering 2 = adapters + quality; 1 = quality; 0 = Ns

# steps 3 to 7
clust = "0.93"                              # the chosen optimum clustering threshold
min_depth = "6"                             # minimum read depth to make a statistical base call (6 during optimization)
min_samps = "3"                             # minimum samples in the full data set having a locus to keep that locus
maxAllcon = "2"                             # maximum number of alleles in individual consensus sequences (2 for diploid)
maxNscon = "0.05"                           # maximum proportion of bases that are Ns in individual consensus sequences
maxHscon = "0.05"                           # maximum proportion of bases that are heterozygous in individual consensus sequences
maxSNPsloc = "0.2"                          # maximum proportion of bases that are SNPs in a locus
maxIndels = "8"                             # maximum number of indels in a locus
maxHshare = "0.5"                           # maximum proportion of samples sharing heterozygous positions for a locus


## Read in the number of cores available as the first command argument, samples list as the second (one per line)
cores = str(sys.argv[1])
print(str(cores) + " CPU cores specified")
sample_list = []
with open(sys.argv[2], "r") as sample_file:
	for sample in sample_file:
		sample_list.append(sample.rstrip())

print("Read in " + str(len(sample_list)) + " samples for analysis")


## start the ipcluster engines and associate them with the ipyclient
# --daemonize puts them into the background; set up a profile name for associating the client and later shutting down
call(["ipcluster", "start", "--n=" + cores, "--daemonize", "--profile=local"])
time.sleep(10)
ipc = ipp.Client(profile="local")
time.sleep(10)
print(ipc.ids)
print(len(ipc), "cores")
ip.cluster_info(ipc)
sys.stdout.flush()


#######
# ipyrad run
#######

# initialise the ipyrad datasets, set parameters needed for demultiplexing and trimming, then demultiplex
# our reads came in four index groups, with barcodes specific to each group
# 1
data1 = ip.Assembly("data1")
data1.set_params("raw_fastq_path", "*ACAGTG*.gz")
data1.set_params("barcodes_path", "barcodes_ACAGTG.tab")
data1.set_params("project_dir", project)
data1.set_params("datatype", datatype)
data1.set_params("restriction_overhang", restrict_overhang)
data1.set_params("max_barcode_mismatch", barcode_mismatch)
data1.run(steps="1", ipyclient=ipc, show_cluster=True)

# 2
data2 = ip.Assembly("data2")
data2.set_params("raw_fastq_path", "*CTTGTA*.gz")
data2.set_params("barcodes_path", "barcodes_CTTGTA.tab")
data2.set_params("project_dir", project)
data2.set_params("datatype", datatype)
data2.set_params("restriction_overhang", restrict_overhang)
data2.set_params("max_barcode_mismatch", barcode_mismatch)
data2.run(steps="1", ipyclient=ipc, show_cluster=True)

# 3
data3 = ip.Assembly("data3")
data3.set_params("raw_fastq_path", "*GCCAAT*.gz")
data3.set_params("barcodes_path", "barcodes_GCCAAT.tab")
data3.set_params("project_dir", project)
data3.set_params("datatype", datatype)
data3.set_params("restriction_overhang", restrict_overhang)
data3.set_params("max_barcode_mismatch", barcode_mismatch)
data3.run(steps="1", ipyclient=ipc, show_cluster=True)

# 4
data4 = ip.Assembly("data4")
data4.set_params("raw_fastq_path", "*GTGAAA*.gz")
data4.set_params("barcodes_path", "barcodes_GTGAAA.tab")
data4.set_params("project_dir", project)
data4.set_params("datatype", datatype)
data4.set_params("restriction_overhang", restrict_overhang)
data4.set_params("max_barcode_mismatch", barcode_mismatch)
data4.run(steps="1", ipyclient=ipc, show_cluster=True)


# merge datasets and run step 2
merge = ip.merge("merge", [data1, data2, data3, data4])
merge.set_params("max_low_qual_bases", max_low_qual)
merge.set_params("phred_Qscore_offset", phred)
merge.set_params("filter_adapters", filter)
merge.set_params("filter_min_trim_len", min_trim_len)
merge.run(steps="2", ipyclient=ipc, show_cluster=True)
print(merge.stats)
sys.stdout.flush()


# branch the dataset to include only samples in the sample list
full = merge.branch("full", subsamples = sample_list)


# set parameters for steps 3 to 7 and run the steps
full.set_params("clust_threshold", clust)
full.set_params("max_alleles_consens", maxAllcon)
full.set_params("max_Ns_consens", maxNscon)
full.set_params("mindepth_statistical", min_depth)
full.set_params("mindepth_majrule", min_depth)
full.set_params("min_samples_locus", min_samps)
full.set_params("max_SNPs_locus", maxSNPsloc)
full.set_params("max_Indels_locus", maxIndels)
full.set_params("max_shared_Hs_locus", maxHshare)
full.set_params("output_formats", ["v"])                        # v: vcf
full.run(steps="34567", ipyclient=ipc, show_cluster=True)
print(full.stats)
sys.stdout.flush()


# run extra step 7 for 50% phylip file
tree_br = full.branch("50_tree")
tree_br.set_params("min_samples_locus", round(len(sample_list)/2))
tree_br.set_params("output_formats", ["p"])                     # p: phylip
tree_br.run(steps="7", ipyclient=ipc, show_cluster=True)


# stop the cluster
call(["ipcluster", "stop", "--profile=local"])
