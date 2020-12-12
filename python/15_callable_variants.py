#! /usr/bin/python3

""" Count per-region callable variants in each set

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

Outputs pileup (different then the pileup files used to calculate
coverage; the unortunate nomenclature is not mine) to investigate
if any variants can be called. The laater is particularly problematic
in the case of the AH set.

Sources:
https://samtools.github.io/bcftools/howtos/variant-calling.html

This this script assumes the directory tree of the format
.
├── alignments
├── clean_fastq
├── fastqc
├── pileup
├── python
└── task_dna

where task_dna contains the fastq files (the name chosen to correspond
with the original task folder).  Other dependencies shoud be set in config.py.
"""

import matplotlib.pyplot as plt
from pileup import *

# fastq-specific functions - imports utils
from fastqc import *


def read_regions(merged_regions_fnm):
	regions = {}
	inf = open(merged_regions_fnm)
	for line in inf:
		fields = line.strip().split()[:3]
		if len(fields) < 3: continue
		[chrom, start, end]  = [int_cast(field.replace("chr",""), merged_regions_fnm) for field in fields]
		if chrom not in regions: regions[chrom] = []
		regions[chrom].append([start, end])
	inf.close()
	return regions


def talk(regions, phred_likelihoods):
	for chrom in range(1,23): # I don have sex chroms in this toy problem
		if chrom not in regions: continue
		print("chrom", chrom)
		for [start, end] in regions[chrom]:
			calls_count = [0, 0]
			for i in range(2):
				if not chrom in phred_likelihoods[i]: continue
				if not start in phred_likelihoods[i][chrom]: continue
				calls_count[i] = len(phred_likelihoods[i][chrom][start])

				print(f"\t {start}     {calls_count}  ")


def write(vcf_dir, regions, phred_likelihoods):
	with open(f"{vcf_dir}/calls_per_interval.tsv", "w") as outf:
		for chrom in range(1,23): # I don have sex chroms in this toy problem
			if chrom not in regions: continue
			for [start, end] in regions[chrom]:
				calls_count = [0, 0]
				for i in range(2):
					if not chrom in phred_likelihoods[i]: continue
					if not start in phred_likelihoods[i][chrom]: continue
					calls_count[i] = len(phred_likelihoods[i][chrom][start])
				print(f"{chrom}\t{start}\t{calls_count[0]}\t{calls_count[1]}", file=outf)


def main():

	home_path = get_home_path()

	bcftools    = Config.bcftools
	alnmts_dir  = f"{home_path}/alignments"
	cvg_dir     = f"{home_path}/pileup/coverage"
	vcf_dir     = f"{home_path}/pileup/variants"
	dna_dir     = f"{home_path}/task_dna"
	rootnames   = Config.paired_reads_rootnames
	merged_regions_fnm = f"{cvg_dir}/merged_target_regions.bed"

	vcf_files   = [f"{vcf_dir}/{rootnm}.dedup.vcf.gz" for rootnm in rootnames]
	check_exist([bcftools, alnmts_dir, dna_dir, vcf_dir, cvg_dir, merged_regions_fnm] + vcf_files)

	regions = read_regions(merged_regions_fnm)

	# summarize coverage/depth for each region
	phred_likelihoods = []
	for i in range(2):
		phred_likelihoods.append(summarize_variant_pileup(bcftools, vcf_files[i], regions))

	talk(regions, phred_likelihoods)

	# let's write this out because it is rather slow - we'll use the info in the xlsx table
	write(vcf_dir, regions, phred_likelihoods)

	return


if __name__ == "__main__":
	main()
