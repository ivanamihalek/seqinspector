#! /usr/bin/python3

""" Getting rid of PCR duplicates using samtools

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

The  order of processsing is the recipe suggested by samtools
http://www.htslib.org/doc/samtools-markdup.html
1) collate
2) fixmate
3) sort (in position order)
4) mark duplicates (confusingly, this also removes PCR duplicates)

Sources:
http://www.htslib.org/doc/samtools.html

This this script assumes the directory tree of the format
.
├── alignments
├── clean_fastq
├── fastqc
├── python
└── task_dna

"""

import subprocess
# setting home path, checing deps, importing config too
# from utils import *
# fastq-specific functions - imports utils
from fastqc import *


def get_bamfiles(alnmts_dir, rootnames):
	bamfiles = []
	for rootnm in rootnames:
		bamfile = f"{alnmts_dir}/{rootnm}.bam"
		if not os.path.exists(bamfile):
			print(f"{rootnm}.bam not found in {alnmts_dir}")
			exit()
		bamfiles.append(bamfile)
	return bamfiles


def collate(samtools, bamfile):
	colltfile = bamfile[:-3] + "collt.bam"
	# the flag -n here indicates sorting by name
	cmd = f"{samtools} collate  -o {colltfile} {bamfile}   2> /dev/null"
	print(cmd)
	subprocess.call(["bash", "-c", cmd])
	# TODO: check output
	return colltfile


def fixmate(samtools, colltfile):
	old_extension = "collt.bam"
	fixmate_bamfile = colltfile[:-len(old_extension)] + "fixmate.bam"
	# -m flag is needed if the output is to be used in markdup
	# note that this one for a change does not use -o option
	cmd = f"{samtools} fixmate -m  {colltfile} {fixmate_bamfile}  2> /dev/null"
	print(cmd)
	subprocess.call(["bash", "-c", cmd])
	# TODO: check output
	return fixmate_bamfile


def coordinate_sort(samtools, fixmate_bamfile):
	old_extension = "fixmate.bam"
	sortfile = fixmate_bamfile[:-len(old_extension)] + "sort.bam"
	cmd = f"{samtools} sort  -o {sortfile} {fixmate_bamfile}   2> /dev/null"
	print(cmd)
	subprocess.call(["bash", "-c", cmd])
	# TODO: check output
	return sortfile


def markdup(samtools, sortfile):
	old_extension = "sort.bam"
	dedup_bamfile = sortfile[:-len(old_extension)] + "dedup.bam"
	cmd = f"{samtools} markdup {sortfile}  {dedup_bamfile}    2> /dev/null"
	print(cmd)
	subprocess.call(["bash", "-c", cmd])
	# TODO: check output
	return dedup_bamfile


def index(samtools, dedupfile):
	cmd = f"{samtools} index {dedupfile}   2> /dev/null"
	print(cmd)
	# TODO: check output
	subprocess.call(["bash", "-c", cmd])


def main():

	home_path = get_home_path()

	samtools   = Config.samtools
	alnmts_dir = f"{home_path}/alignments"
	rootnames  = Config.paired_reads_rootnames

	check_exist([samtools, alnmts_dir])

	# note that the bamfile has to be sorted before we go into dedup
	bamfiles = get_bamfiles(alnmts_dir, rootnames)

	print("bamfiles:", bamfiles)
	for bamfile in bamfiles:
		print()
		# pipeline:
		prev_file = bamfile
		for step in [collate, fixmate, coordinate_sort, markdup, index]:
			print(step.__name__)
			new_file  = step(samtools, prev_file)
			# leave the orignal bam for now; the last step is indexing which does not change the penultimate ba,
			if step.__name__ not in ["collate", "index"]:  os.remove(prev_file)
			prev_file = new_file


	return


if __name__ == "__main__":
	main()
