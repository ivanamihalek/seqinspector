#! /usr/bin/python3

""" Read quality check using FastQC

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

This script runs fastqc and stores the results in the fastqc/first_pass
directory, then prints out the FAIL and WAR fields from the summary. The
html reports can also be found in the same directory after this is run.

Sources:
https://github.com/s-andrews/FastQC
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

This this script assumes the directory tree of the format
.
├── python
└── task_dna

where task_dna contains the fastq files (the name chosen to correspond
with the original task folder).  This script will add directory called
"fastqc." Other dependencies shoud be set in config.py.

"""

# setting home path, checing deps, importing config too:
# from utils import *
# fastq-specific functions - imports utils:
from fastqc import *


def main():
	home_path = get_home_path()

	fastqc = Config.fastqc
	unzip  = Config.unzip
	dna_dir = f"{home_path}/task_dna"
	out_dir = f"{home_path}/fastqc/first_pass"
	rootnames   = Config.paired_reads_rootnames
	read_labels = Config.read_labels

	if not os.path.exists(out_dir): os.makedirs(out_dir)

	dependencies = [fastqc, unzip, dna_dir, out_dir]
	check_exist(dependencies)
	check_fastq_exist(dna_dir)

	run_fastqc(dependencies, rootnames, read_labels)

	# check the summary file for each fastq
	# print issues to stdout
	talk_fastqc(rootnames, read_labels, out_dir)

	# for the report
	# fastqc_latex_table(rootnames, read_labels, out_dir)

	return


if __name__ == "__main__":
	main()
