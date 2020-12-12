#! /usr/bin/python3

""" Rerun FastQC on deduplicated sequences

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

Sources:
https://github.com/s-andrews/FastQC
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

This this script assumes the directory tree of the format
.
├── alignments
├── clean_fastq
├── fastqc
├── python
└── task_dna

It will make and remove a scratch directtry. It unpacks the
de-duplicated sam into two fastq file, to maintain the
format of input/outpud used on raw and trimmed seqeunces.

"""

# setting home path, checing deps, importing config too
# from utils import *
# fastq-specific functions - imports utils
from fastqc import *


def main():

	home_path = get_home_path()

	samtools        = Config.samtools
	fastqc          = Config.fastqc
	unzip           = Config.unzip
	alignments      = f"{home_path}/alignments"
	fastqc_this_dir = f"{home_path}/fastqc/dedup"  # the output in this round
	scratch         = f"{home_path}/scratch"

	rootnames   = Config.paired_reads_rootnames
	read_labels = Config.read_labels
	dedup_bams  = [f"{alignments}/{r}.dedup.bam" for r in rootnames]

	check_exist([samtools, alignments] + dedup_bams)

	os.makedirs(fastqc_this_dir, exist_ok=True)
	os.makedirs(scratch, exist_ok=True)

	for rootname in rootnames:
		cmd  = f"{samtools} fastq  "
		cmd += f"-1 {scratch}/{rootname}_{read_labels[0]}.fastq -2 {scratch}/{rootname}_{read_labels[1]}.fastq "
		cmd += f"-0 /dev/null -s /dev/null  -n  {alignments}/{rootname}.dedup.bam"
		print(cmd)
		subprocess.call(["bash", "-c", cmd])

	dependencies = [fastqc, unzip, scratch, fastqc_this_dir]
	run_fastqc(dependencies, rootnames, read_labels)
	talk_fastqc(rootnames, read_labels, fastqc_this_dir)
	# for the report
	# fastqc_latex_table(rootnames, read_labels, fastqc_this_dir)

	shutil.rmtree(scratch)

	return


if __name__ == "__main__":
	main()
