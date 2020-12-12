

""" Utilities for fastqc

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

"""

from utils import *
import subprocess, os, shutil


def run_fastqc(dependencies, paired_reads_root, read_labels):
	[fastqc, unzip, dna_dir, out_dir] = dependencies
	print("running fastqc:")
	for root in paired_reads_root:
		for label in read_labels:
			fastq = f"{dna_dir}/{root}_{label}.fastq"

			# run fastqc analyzer; options: q is for quiet, o outdir
			cmd = f"{fastqc} {fastq}  -q -o {out_dir}"
			print(cmd)
			subprocess.call(["bash", "-c", cmd])
			# unzip the directory which we produced to access the summary file
			# TODO check that the zip file was indeed produced
			# remove old version - the name is the convention that fastqc uses
			fqcdir = fastqc_dir(out_dir, root, label)
			# there are some security issues with this function; outse of toy context organize differently
			# https://docs.python.org/3/library/shutil.html#shutil.rmtree
			if os.path.exists(fqcdir): shutil.rmtree(fqcdir)
			# qq must come before the input (zip) file
			cmd = f"{unzip} -qq {fqcdir}.zip  -d {out_dir}"
			print(cmd)
			subprocess.call(["bash", "-c", cmd])
	print("fastqc done\n")


def talk_fastqc(rootnames, read_labels, out_dir):
	for root in rootnames:
		for label in read_labels:
			print()
			print(root, label)
			fqcdir = fastqc_dir(out_dir, root, label)
			with open(f"{fqcdir}/summary.txt") as inf:
				flaglines  = sorted([line for line in inf if "WARN" in line or "FAIL" in line])
			if not flaglines:
				print("fastq reports no issues")
			else:
				[print(line, end="") for line in flaglines]


def fastqc_latex_table(rootnames, read_labels, out_dir):
	flags = {}
	issues = set()
	for root in rootnames:
		for label in read_labels:
			fqcdir = fastqc_dir(out_dir, root, label)
			with open(f"{fqcdir}/summary.txt") as inf:
				flaglines  = sorted([line for line in inf if "WARN" in line or "FAIL" in line])
			if not flaglines: continue
			for line in flaglines:
				[flag, issue] = line.split("\t")[:2]
				issues.add(issue)
				if root not in flags: flags[root] = {}
				if label not in flags[root]: flags[root][label] = {}
				flags[root][label][issue] = flag

	for issue in sorted(issues):
		print(issue, end=" & ")
		columns  = []
		for root in rootnames:
			for label in read_labels:
				if root in flags and label in flags[root] and issue in flags[root][label]:
					flag = flags[root][label][issue]
				else:
					flag = "OK"
				columns.append(flag)
		print("  & ".join(columns), end=" \\\\ \n")
