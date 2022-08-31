***Geleznowia* ddRAD analyses**  
Author: B.M. Anderson  

These notes outline the analyses run on ddRAD data to generate results for the paper "[Title]"  
All scripts are assumed to be in an accessible folder for relative calls, indicated here as a directory in home called `~/scripts/`  


# Data
ddRAD enzymes: PstI and MspI  
Illumina NovaSeq 6000, 150 bp single end reads, four indices: ACAGTG, CTTGTA, GCCAAT, GTGAAA (2 files per index)  

Barcodes are per index, so there will be four sets to demultiplex (and four barcode files needed)  
Create a tab-delimited text file (`barcodes.tab`) with columns corresponding to (no header): sampleID seqID index barcode
```s
for index in ACAGTG CTTGTA GCCAAT GTGAAA
do
# grab lines with the index
grep "$index" barcodes.tab > temp
# join the first two columns with a "-" and add a column with the barcode
paste <(cut -f 1,2 temp | tr "\t" "-") <(cut -f 4 temp) > barcodes_"$index".tab
done
rm temp
```
This will produce four barcodes files, one per index  


# Assembly in ipyrad
ipyrad v. 0.9.81 running in a Singularity container (ipyrad installed with miniconda)  
Depending on the high performance computing setup, the method of executing the run will change  
The run was executed on the Zeus cluster on Pawsey (no longer running), using an sbatch script to call the container to execute a Python script that specified the run parameters (`ipyrad_full.py`)  

An initial run showed that four samples had too few reads to be retained: MAR06-06-343713, ZPP01-05_R-343786, ZPS01-04_R-343787, ZPS01-06-343778  
Based on preliminary optimisation runs on subsets of the data, the optimal clustering threshold was determined to be 0.93  

The full run therefore used the optimal clustering threshold and was subset to exclude the above four samples  
The run used 28 CPUs and took roughly 12.5 hours  
The command to run was roughly as follows, with the read and barcode files in the working directory
```s
singularity exec -B "$(pwd)":/home/"$USER" ipyrad.sif python ipyrad_full.py 28 samples_full.txt
```
The primary output for downstream analyses is the VCF file: `full.vcf`  
The VCF file is largely unfiltered, keeping all loci found in at least three samples  
Another output, `full.loci`, can be used to access all sites (not just SNPs)  


# Filtering
Datasets were filtered for each analysis differently, and were informed by an initial estimate of error rate based on technical replicates  
After estimating error and clonality (see below), one of each of the technical replicates can be removed (chosen manually based on more recovered loci); these sampleIDs can be put in a `reps.to_remove` text file, one per line for later filtering steps  


## Tiger genotyping error estimate
Create a subset of the VCF for the replicate pairs (one of the reps in each pair is named with a "_R" subscript for sampleID) 
First, create a list of the samples
```s
grep "_R" samples_full.txt | cut -f 1 -d "_" > temp
grep -f temp samples_full.txt | sort > reps.txt
rm temp
```

Now filter the VCF for just those samples, as well as for a minimum sample coverage of 50%, a minimum mean depth of 10 and only biallelic SNPs using the script `filter_vcf.py`  
Then use the script `single_snp.py` to select a single SNP per locus (max depth or random in the case of ties in coverage)  
```s
grep -v -f reps.txt samples_full.txt > samples.to_remove 
python ~/scripts/filter_vcf.py -o error --mincov 0.5 --minmd 10 --bial yes -s samples.to_remove full.vcf
python ~/scripts/single_snp.py -r yes error.vcf
mv mod_error.vcf error.vcf
```

Using the technical replicates and the software Tiger (https://bitbucket.org/wegmannlab/tiger/wiki/Home), estimate error rates  
Create an input samples file for Tiger (needs reps assigned to groups)  
```s
mkdir tiger && cd tiger
num=$(($(wc -l < reps.txt) / 2))
for ((i=1; i<=num; i++)); do echo $i >> temp && echo $i >> temp; done
echo -e "Sample\tGroup" > rep_groups.txt
paste reps.txt temp >> rep_groups.txt
rm temp
```

Now run error estimation and visualise it with the script `Tiger_error_plot.py`  
```s
tiger task=estimateIndReps outname=error vcf=error.vcf groups=rep_groups.txt
python ~/scripts/Tiger_error_plot.py -o error -m 500 error_errorRates.txt
```
This will create an `error.pdf` file  
Based on the error rate estimates, mean read depth cutoffs of min 17 and max 100 were chosen  


## Clonality
To assess whether samples are indistinguishable from clones, technical replicates can again be used to establish what level of similarity is expected for two identical samples  
This can be based on the proportion of sites called in both that differ  

First, filter the VCF for only ingroup samples and for SNPs found in at least half the samples (to keep the dataset reasonable in size)  
In our naming convention, outgroup samples have a "Z" in their sampleIDs  
The script to assess similarity is `vcf_similarity.py`  
```s
grep "Z" samples_full.txt > samples.to_remove
python ~/scripts/filter_vcf.py -o clones --mincov 0.5 -s samples.to_remove full.vcf
python ~/scripts/vcf_similarity.py -v clones.vcf -o clones
```
The output `clones_comps.txt` can be imported into a spreadsheet program and sorted to manually check  
The output `clones_hist.png` can be useful to see the distribution of comparisons and whether there is a clear break between the replicate comparisions and other comparisons in the dataset  

Using the average between the pairs of technical replicates (0.004) and doubling it (0.008), all comparisons less than that were considered essentially clones  
Create a text file (`clones.txt`) with groups of sampleIDs that are clones (one sampleID per line), each group separated by a blank line  
Note: ensure that any technical reps expected to be removed (in the `reps.to_remove` file) are not in this `clones.txt` file  

Now randomly select a single representative per clone group to keep  
```s
awk -v RS= '{print > ("clone_" NR ".txt")}' clones.txt
for file in clone_*; do shuf -n 1 "$file" >> clones.to_keep; done
for file in clone_*; do grep -v -f clones.to_keep "$file" >> clones.to_remove; done
rm clone_*
```


## Set 1: ingroup (strict filtering)
For analyses sensitive to missing data (PCA, Structure), the VCF was filtered to only ingroup samples and biallelic SNPs present in at least 90% of samples (one per locus), with a minimum minor allele count of 3  
Create a list of samples that should be excluded (including clones, replicates and outgroups)  
```s
cat clones.to_remove > samples.to_remove
cat reps.to_remove >> samples.to_remove
grep Z samples_full.txt >> samples.to_remove
sort samples.to_remove | uniq > temp
mv temp samples.to_remove
```

Filter
```s
python ~/scripts/filter_vcf.py -o set1 --mincov 0.9 --minmd 17 --maxmd 100 --mac 3 --bial yes -s samples.to_remove full.vcf
python ~/scripts/single_snp.py -r yes set1.vcf
mv mod_set1.vcf set1.vcf
```


## Set 2: ingroup (lax filtering)
For the distance network analysis, there was less sensitivity to missing data, so the VCF was filtered in a similar way to Set 1 but only requiring SNPs be present in 25% of samples (still one per locus), with no restrictions on allele count  
```s
python ~/scripts/filter_vcf.py -o set2 --mincov 0.25 --minmd 17 --maxmd 100 --mac 1 -s samples.to_remove full.vcf
python ~/scripts/single_snp.py -r yes set2.vcf
mv mod_set2.vcf set2.vcf
```


## Set 3: phylo (concatenation)
To generate alignments, loci were selected after removing putative hybrids (from results of network and Structure analyses) and some of the possible outgroups (retaining two), keeping loci present in at least 50% of samples  
Create a list of samples to remove  
```s
cat clones.to_remove > samples.to_remove
cat reps.to_remove >> samples.to_remove
grep "BIN01-02" samples_full.txt >> samples.to_remove
grep -E "BIN02-0[1,3,5]" samples_full.txt >> samples.to_remove
grep "VER04" samples_full.txt >> samples.to_remove
grep ZP[S,P] samples_full.txt >> samples.to_remove
sort samples.to_remove | uniq > temp
mv temp samples.to_remove
```

Filter
```s
python ~/scripts/filter_vcf.py -o phylo --mincov 0.5 --minmd 17 --maxmd 100 --mac 1 -s samples.to_remove full.vcf
```

From the filtered VCF, loci names can be extracted and used to extract the full alignments from `full.loci` using the script `loci_extract.py`  
```s
grep -v "#" phylo.vcf | cut -f 1 | uniq | sed 's/RAD_//g' > loci.txt
mkdir loci_extract && cd loci_extract
python ~/scripts/loci_extract.py -l ../loci.txt ../full.loci
```

Now, the loci files need to have any taxa removed (using the script `remove_fastas.py`) that were removed in the VCF filtering  
```s
for locus in locus*.fasta; do
python ~/scripts/remove_fastas.py -f ../samples.to_remove "$locus"
mv mod_"$locus" "$locus"
done
```

All the loci can be combined into a single alignment using the script `combine_alignments.py`  
In addition, the alignment needs to then be filtered to remove any positions with more than 50% missing data (using the script `clean_alignment.py`)
```s
python ~/scripts/combine_alignments.py -f single *.fasta
python ~/scripts/clean_alignment.py -p 50 combine_out.fasta
mv combine_out_clean.fasta ../concat_loci.fasta
cd .. && rm -r loci_extract
```
The resulting file is `concat_loci.fasta`  


## Set 4: phylo (coalescent)
For running SVDquartets, extract a single SNP per locus from the same loci as were used for concatenation (present in at least 50% of samples)
```s
python ~/scripts/single_snp.py -r yes phylo.vcf
```
This will create the file `mod_phylo.vcf`  


# Analyses
Scripts to assess each of the input VCF files for missing data and read depth are `missing_plot.R` and `vcf_depth.py`  

Population numbering is based on the south-north distribution of the populations  
```
P#	Sample designator
P1	MAR10
P2	MAR09
P3	MAR08
P4	MAR07
P5	VER09
P6	MAR06
P7	MAR01
P8	VER08
P9	WHT01
P10	VER07
P11	MAR05
P12	VER06
P13	VER05
P14	BIN02
P15	VER04
P16	BIN01
P17	VER03
P18	AMA03
P19	AMA02
P20	AMA01
P21	VER02
P22	VER01
P23	MAR04
P24	MAR03
P25	MAR02
```


## PCA
With dataset 1 (strict filtering) as input, the PCA can be run interactively using `pca.rmd`  
Copy or link the `set1.vcf` into the working directory where the PCA is run  


## Distances
With dataset 2 (lax filtering) as input, run the script `distances.R` to generate a distances matrix  
Use the option `-d M` for MATCHSTATES distances from the package `pofadinr`  
```s
Rscript ~/scripts/distances.R -d M -o set2 -v set2.vcf
```
This will produce a number of files, but the key output is `set2_distMATCHSTATES.nex`, which can be input into SplitsTree4 v. 4.17.1 to create a NeighborNet network  


## Structure
Using dataset 1 (strict filtering), convert the VCF to an input format appropriate for Structure, including the generation of a population mapping file to order groups (based on results of PCA and Distances)  
```s
grep "#CHROM" set1.vcf | cut -f 1-9 --complement | tr -s "\t" "\n" > temp
paste temp <(cut -f 1 -d "-" temp) > str_pops.tab && rm temp
sed -i 's/VER0[3-8]$/1/g' str_pops.tab
sed -i 's/BIN0[1-2]$/2/g' str_pops.tab
sed -i 's/MAR0[1,6,9]$/3/g' str_pops.tab
sed -i 's/VER0[1-2]$/3/g' str_pops.tab
sed -i 's/VER09$/3/g' str_pops.tab
sed -i 's/MAR08$/4/g' str_pops.tab
sed -i 's/MAR10$/4/g' str_pops.tab
sed -i 's/MAR0[2-5]$/5/g' str_pops.tab
sed -i 's/AMA0[1-3]$/6/g' str_pops.tab
sed -i 's/MAR07$/7/g' str_pops.tab
sed -i 's/WHT01$/8/g' str_pops.tab
```

Use the mapping file and generate the Structure input file with the script `vcf_to_structure.py`  
```s
python ~/scripts/vcf_to_structure.py -p str_pops.tab set1.vcf
mv set1.vcf.str set1.str
```

Use the `mainparams` and `extraparams` files from Structure (see https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html), and alter specific lines to (as appropriate for the input file just generated):
```s
# mainparams
BURNIN	100000
NUMREPS	100000
NUMINDS	84
NUMLOCI	1683
MARKERNAMES	0

#extraparams
UPDATEFREQ	1000
RANDOMIZE	0
```

Running Structure can be done on a cluster, submitting multiple jobs, one per K, each run using the input file `set1.str` and the two params files  
In this case, Structure was run for 40 replicates (in parallel) per K from K1 to K15 (15 40-core jobs)  
Each rep in each job called Structure like so:
```s
structure -m mainparams -e extraparams -K "$k_val" -i "$str_file" -o output_"$k_val"_"$rep" -D "$rep""$k_val""$RANDOM"
```
The jobs took from 24 to 74 minutes  

It is best to download the results to a separate folder, as there will be many files  

For reducing file size, remove all the extraneous "Locus" information from the output files  
```s
# from within the output folder, having files named in the pattern "output_x_x_f"
for file in *_f; do
sed -n '/Locus/q;p' $file > temp && mv temp $file
done
```

Evaluate the likelihoods of the replicates and keep the top 10 per K
```s
# grab likelihoods from the output files
for K in {1..15}; do
for rep in {1..40}; do
grep "Estimated Ln Prob" output_"$K"_"$rep"_f | cut -f 2 -d "=" | sed "s/^[[:space:]]*/$K\t/" >> k${K}_likes.txt
done
done

for K in {1..15}; do
# number the likelihoods to keep track of which reps they refer to
num=$(wc -l < k"$K"_likes.txt)
for ((i=1; i<=num; i++)); do
echo $i >> temp.txt
done
paste k"$K"_likes.txt temp.txt > temp_lines.txt && rm temp.txt
# sort the likelihoods and keep the top 10
sort -k2 -nr temp_lines.txt | cut -f 3 | head -n 10 > temp_select.txt
# grab those reps and lines from the likelihood files for plotting
for row in $(cat temp_select.txt); do
cp output_"$K"_"$row"_f select_"$K"_"$row"_f
done
cut -f 1,2 <(sort -k2 -nr temp_lines.txt) | head -n 10 > select_k"$K"_likes.txt
done
rm temp_lines.txt temp_select.txt
```

Evaluate the best K based on the Evanno method using the script `bestK_Evanno.py`  
```s
cat k*likes.txt > all_likes.txt
cat select*likes.txt > all_likes_selected.txt
python ~/scripts/bestK_Evanno.py -l all_likes.txt -o set1_bestK_all
python ~/scripts/bestK_Evanno.py -l all_likes_selected.txt -o set1_bestK_selected
```

Upload the top 10 per K to CLUMPAK (http://clumpak.tau.ac.il/index.html) to run with default settings
```s
zip results.zip select_*_f
# upload results.zip
```

Download, then create qfiles from the major modes for plotting  
Note: in the download, there should be a summary logfile indicating how many minor modes there are per K  
Wherever the downloaded folder is (has subfolders for each K), store that path in a variable e.g. `output=path/to/output`  
```s
for K in {1..15}; do
cut -f 2 -d ':' ${output}/K\=${K}/MajorCluster/CLUMPP.files/ClumppIndFile.output | awk '{$1=$1; print}' > qfile"$K".txt
done
```

If there are minor modes of interest, they can also be kept as qfiles for plotting  
```s
for K in {1..15}; do
for minor_num in {1..4}; do
if [ -d ${output}/K\=${K}/MinorCluster${minor_num} ]; then
cut -f 2 -d ':' ${output}/K\=${K}/MinorCluster${minor_num}/CLUMPP.files/ClumppIndFile.output | awk '{$1=$1; print}' > qfile"$K"_minor"$minor_num".txt
fi
done
done
```

Plot modes of interest using the script `structure_barplots.py`  
Colours for the bars are specified in a `colours.txt` file, one colour (in hex code) per line  
Note: change the qfile and `-o` arguments for the minor mode files, if desired  
```s
for num in {2..8}; do
python ~/scripts/structure_barplots.py -o K"$num" -p str_pops.tab -q qfile"$num".txt -c colours.txt
done
```
The resulting barplots can be combined manually into a single figure in Inkscape  


## ML concatenation
Run the concatenated loci (`concat_loci.fasta`) in IQ-TREE v. 2.1.3, running a search for the best model, 1000 ultrafast bootstrap replicates and Shimodaira-Hasegawa-like approximate likelihood ratio tests  
Here it is run using a maximum of 8 cores on a laptop (took about 6 hours)  
```s
iqtree -s concat_loci.fasta --prefix phylo --threads-max 8 -T AUTO --seqtype DNA -m MFP --merit BIC --ufboot 1000 --alrt 1000
```

The resulting tree (`phylo.treefile`) can be interactively plotted using `plot_tree.rmd`  
Create an `outgroup.txt` file for the process  
```s
grep ">Z" concat_loci.fasta | cut -f 1 -d " " | sed 's/>//g' > outgroup.txt
```

If wanting to collapse clades, make a `clades.txt` file for clades of samples (search terms) that should be collapsed, one per line, tab separated if more than one search term per clade  

It can also be insightful to look at splits present in the bootstrap trees  
Use the `phylo.splits.nex` from the run in SplitsTree4 v. 4.17.1 to create a network, using the weights of the splits (bootstrap occurrences), which can be displayed to show other relationships recovered in the bootstrap trees  


## SVDquartets
Use the script `vcf_to_svd.py` to convert the filtered VCF (`mod_phylo.vcf`) into the correct format (and as a single sequence with ambiguities rather than two haplotypes)  
First, set up the hypothesised species to run a grouping of lineages if desired  
The sampleIDs are such that pop designators come first followed by a "-", so they can be modified to putative groups (based on other analyses)  
```s
grep "#CHROM" mod_phylo.vcf | cut -f 1-9 --complement | tr -s "\t" "\n" > temp
paste temp <(cut -f 1 -d "-" temp) > spp.tab && rm temp
sed -i 's/VER0[3-8]$/G1/g' spp.tab
sed -i 's/BIN0[1-2]$/G2/g' spp.tab
sed -i 's/VER0[1-2]$/G3/g' spp.tab
sed -i 's/VER09$/G3/g' spp.tab
sed -i 's/MAR0[1,6,9]$/G3/g' spp.tab
sed -i 's/MAR08$/G4/g' spp.tab
sed -i 's/MAR10$/G4/g' spp.tab
sed -i 's/MAR0[2-5]$/G5/g' spp.tab
sed -i 's/AMA0[1-3]$/G6/g' spp.tab
sed -i 's/MAR07$/P4/g' spp.tab
sed -i 's/WHT01$/P9/g' spp.tab
sed -i 's/ZDE0[1-2]$/Dericoides/g' spp.tab
sed -i 's/ZPC01$/Psericea/g' spp.tab
```

Use that for creating the input NEXUS file (`mod_phylo.nex`)  
```s
python ~/scripts/vcf_to_svd.py -v mod_phylo.vcf -f seq -s spp.tab
```

Now create batch files (text files called `batch.nex`) for running PAUP, the first for lineages  
```s
#NEXUS
BEGIN PAUP;
	LOG start file=log.txt;
	SET autoclose=yes;
	exe mod_phylo.nex;
	svdq nthreads=8 evalQuartets=all treeInf=QFM bootstrap=no ambigs=distribute;
	describeTrees;
	saveTrees file=qfm.tre format=Newick;
	svdq nthreads=8 evalQuartets=all treeInf=QFM bootstrap=standard nreps=100 treefile=boots.tre ambigs=distribute;
	describeTrees;
	saveTrees file=bootcon.tre supportValues=nodeLabels format=Newick;
	LOG stop;
	QUIT;
END;
```

and the second for hypothesised species  
```s
#NEXUS
BEGIN PAUP;
	LOG start file=log.txt;
	SET autoclose=yes;
	exe mod_phylo.nex;
	svdq nthreads=8 taxpartition=svdtaxa evalQuartets=all treeInf=QFM bootstrap=no ambigs=distribute;
	describeTrees;
	saveTrees file=qfm.tre format=Newick;
	svdq nthreads=8 taxpartition=svdtaxa evalQuartets=all treeInf=QFM bootstrap=standard nreps=100 treefile=boots.tre ambigs=distribute;
	describeTrees;
	saveTrees file=bootcon.tre supportValues=nodeLabels format=Newick;
	LOG stop;
	QUIT;
END;
```

Each batch file should be moved to its own directory along with a copy of (or link to) the nexus input file  
This is for running PAUP* 4.0a (build 168) for Unix/Linux in a Singularity container (details elsewhere), but it gives the rough idea for running PAUP  
```s
singularity exec paup.sif paup batch.nex
```
(this would be run twice for the lineages and species runs)  

The resulting trees can be plotted using the same `plot_tree.rmd`  

For the species run, create an `outgroup.txt` manually with the lines `Dericoides` and `Psericea`  
Plot the `qfm.tre`, then also the `bootcon.tre` to get bootstrap values; if any clades differ, the bootstrap support can be obtained from the log file  

For the lineages run, create an `outgroup.txt` and `clades.txt` for collapsing groups  
```s
cut -f 1 mod_phylo.nex | grep "Z" > outgroup.txt
```
Note that bootstrap support is not in the `qfm.tre` file, so it can be input using the values from the log file or the `bootcon.tre` file (though sometimes the branching will differ from the inferred tree)  

If wanting to show site concordance factors in addition to bootstrap values, they can be calculated with IQ-TREE (for the lineages run)  
Note that the output trees need to be adjusted to have branch lengths to run properly (just set to 1 in this case)  
```s
sed "s/,/:1,/g" qfm.tre | sed "s/)/:1)/g" > qfm_mod.tre
iqtree -s mod_phylo.nex -t qfm_mod.tre --scf 100000 -T 8 --prefix concord
```
The output can be extracted into a tree form usable by the `plot_tree.rmd` using the script `concord_to_newick.py`  
```s
python ~/scripts/concord_to_newick.py -t concord.cf.tree.nex -o concord
```
This produces the file `concord_scf.tre` to input as a treefile to `plot_tree.rmd`  


It can be insightful to look at the bootstrap tree splits, which can be calculated using the `boots.tre` files from both runs as input to SplitsTree4 v. 4.17.1 to construct Consensus Networks showing splits present in at least 10 trees and weighted by count  


## Population diversity?
