

""" Small utilities

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

"""

import os, sys
from config import Config

def get_home_path():
	script_path = sys.argv[0]
	return os.path.sep.join(os.path.abspath(script_path).split(os.path.sep)[:-2])


def check_exist(dependencies):
	for dep in dependencies:
		if not os.path.exists(dep):
			print(f"{dep} not found. Check the config.py file to set the paths.")
			print("Also note the expected directory tree for this project:", Config.dirtree)
			print("The pipeline should be run in the enumerated order ")
			print("- some files might be missing because they were not created yet.")
			exit()


def check_fastq_exist(dna_dir):
	for root in Config.paired_reads_rootnames:
		for label in Config.read_labels:
			fastq = f"{dna_dir}/{root}_{label}.fastq"
			if not os.path.exists(fastq):
				print(f"{fastq} not found.")
				exit()


def fastqc_dir(fastqc_out_dir, root, label):
	return f"{fastqc_out_dir}/{root}_{label}_fastqc"


def int_cast(x, fnm):
	try:
		return int(x)
	except:
		print(f"unexpected non-integer in {fnm}")
