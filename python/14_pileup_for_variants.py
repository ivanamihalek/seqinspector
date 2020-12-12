#! /usr/bin/python3

""" Pileup for variants

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
with the original task folder). This script will add directory called
"pileup/variants."  Other dependencies shoud be set in config.py.

"""

import matplotlib.pyplot as plt
from pileup import *

# fastq-specific functions - imports utils
from fastqc import *


def bcf_index(samtools, ref_genome):
	if not os.path.exists(f"{ref_genome}.fai"):
		cmd = f"{samtools}  faidx  {ref_genome} 2> /dev/null"
		print(cmd)
		subprocess.call(["bash", "-c", cmd])
	else:
		print(f"found {ref_genome}.fai")


def main():

	home_path = get_home_path()

	samtools    = Config.samtools
	bcftools    = Config.bcftools
	ref_genome  = Config.reference_fasta
	alnmts_dir  = f"{home_path}/alignments"
	dna_dir     = f"{home_path}/task_dna"
	rootnames   = Config.paired_reads_rootnames
	bamfiles    = [f"{alnmts_dir}/{rootnm}.dedup.bam" for rootnm in rootnames]
	region_fnms = [f"{dna_dir}/{rootnm.replace('L001','target')}.txt" for rootnm in rootnames]
	check_exist([samtools, bcftools, alnmts_dir, dna_dir] + bamfiles + region_fnms)

	pileup_dir = f"{home_path}/pileup/variants"
	os.makedirs(pileup_dir, exist_ok=True)

	# bcftools wants reference file as an input, with its own indexing (faidx)
	# which is, of course, produced by samtools
	bcf_index(samtools, ref_genome)

	# run bcftools pileup on each of the regions
	output_pileup_vcf(bcftools, ref_genome, pileup_dir, bamfiles)

	return


if __name__ == "__main__":
	main()
