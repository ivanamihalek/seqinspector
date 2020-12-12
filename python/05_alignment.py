#! /usr/bin/python3

""" Alignment to reference fasta using bwa

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

This script alignes the fastq to the reference genome. Needed for
the following step which is duplicate removal. It wil use
trimmed fastq if available, otherwose fall back on the original
fastq.

Sources:
https://https://github.com/lh3/bwa

This this script assumes the directory tree of the format
.
├── clean_fastq
├── fastqc
├── python
└── task_dna

where task_dna contains the fastq files (the name chosen to correspond
with the original task folder). This script will ad directory called
"alignments."  Other dependencies shoud be set in config.py.

"""



import subprocess

from utils import *


def find_fastq(rootnames, read_labels, dna_dir, clean_dir):
	clean_fastq = {}
	# do we have trimmed version for any of the sequences?
	for root in rootnames:
		clean_fastq[root] = []
		for label in read_labels:
			# let's call this convention over configuration
			orig_version    = f"{dna_dir}/{root}_{label}.fastq"
			trimmed_version = f"{clean_dir}/{root}_trimmed_{label}.fastq"
			if os.path.exists(trimmed_version):
				clean_fastq[root].append(trimmed_version)
			elif os.path.exists(orig_version):
				clean_fastq[root].append(orig_version)
			else:
				print(f"fastq file for {root} {label} not found. I looked in {clean_dir} and in {dna_dir}.")
				exit()
	return clean_fastq


def check_reference_indexed(reference_fasta):
	# this is a rough check - it does not take care of  situations like .sa file exists but is broken
	# (.sa file is the last file created)
	if not os.path.exists(reference_fasta+".sa"):
		print(f"index file for {reference_fasta} not found - indexing")
		cmd = f"{bwa} index -a bwtsw  {reference_fasta}"
		subprocess.call(["bash", "-c", cmd])


def align(bwa, reference_fasta, clean_fastq, alignments, root):
	samfile = f"{alignments}/{root}.sam"
	# bwa will look for input_reference_fasta.sa etc for its indexed input
	cmd = f"{bwa} mem {reference_fasta} {clean_fastq[root][0]} {clean_fastq[root][1]} > {samfile} 2> /dev/null"
	print(cmd)
	subprocess.call(["bash", "-c", cmd])
	# TODO: check output
	return samfile


def sam2bam(samtools, samfile):
	# note that bwa does  have @SQ in the header, otherwise this command is not quite appropriate
	# see http://www.htslib.org/doc/samtools-view.html
	bamfile = samfile[:-3] + "bam"
	cmd = f"{samtools} view -bS {samfile} > {bamfile} 2> /dev/null"
	print(cmd)
	subprocess.call(["bash", "-c", cmd])
	# TODO: check output
	return bamfile


def main():
	home_path = get_home_path()

	bwa = Config.bwa
	samtools = Config.samtools
	reference_fasta = Config.reference_fasta
	dna_dir     = f"{home_path}/task_dna"
	clean_dir   = f"{home_path}/clean_fastq"
	rootnames   = Config.paired_reads_rootnames
	read_labels = Config.read_labels
	alnmts_dir  = f"{home_path}/alignments"

	dependencies = [bwa, samtools, reference_fasta, dna_dir]
	check_exist(dependencies)
	check_reference_indexed(reference_fasta)  # e.g. hg19, was it indexed?

	# check whether we have trimmed files
	clean_fastq = find_fastq(rootnames, read_labels, dna_dir, clean_dir)
	if not os.path.exists(alnmts_dir): os.mkdir(alnmts_dir)

	for root in rootnames:

		# make good ol sam
		samfile = align(bwa, reference_fasta, clean_fastq, alnmts_dir, root)

		# convert sam to bam = compress
		bamfile = sam2bam(samtools, samfile)
		os.remove(samfile)


	return

if __name__ == "__main__":
	main()