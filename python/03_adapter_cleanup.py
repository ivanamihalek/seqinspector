#! /usr/bin/python3

""" Adapter removal using cutadapt

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

This script removes adapters from the fastq files that were flagged
in the previous step (01_fastq_quality_check.py), and re-runs FastQC
on the trimmed seqeunces. The trimmed seqeunces are stored in
clean_fastq directory, and the new FastQC report in the fastqc/trimmed.

Sources:
https://github.com/marcelm/cutadapt/
https://cutadapt.readthedocs.io/en/stable/

This this script assumes the directory tree of the format
.
├── fastqc
├── python
└── task_dna

where task_dna contains the fastq files (the name chosen to correspond
with the original task folder). This script will add directory called
"clean_fastq."  Other dependencies shoud be set in config.py.

"""


# setting home path, checing deps, importing config too
# from utils import *
# fastq-specific functions - imports utils
from fastqc import *


def parse_fastqc_output(fastqc_data_file):
	retmsg = "ok"
	inf = open(fastqc_data_file)
	module_start = ">>Adapter Content"
	module_end   = ">>END_MODULE"
	reading = False
	adapter_content = {}
	column_name = []
	for line in inf:
		if line[:len(module_start)] == module_start:
			if not "fail" in line: break # we're ok
			reading = True
			continue
		if reading:
			if module_end in line: break  # we're done
			if line[0]=="#": # this is the header of the table
				column_name = 	line[1:].strip().split("\t")
				# the first column is position
				for cn in column_name[1:]: adapter_content[cn] = 0
			else:
				if len(column_name)==0:
					print(f"No header for adapter content table found in{fastqc_data_file}.")
					exit()
				field = line.strip().split("\t")
				for i in range(1,len(column_name)):
					adapter_content[column_name[1]] += float(field[i])

	if adapter_content:
		# TODO what's with this cutoff - what if more than one adapter shows up?
		notable_adapters = [k for k in adapter_content.keys() if adapter_content[k]>1]
		if len(notable_adapters)>1:
			print(f"More than one adapter content notable in {fastqc_data_file}: {notable_adapters}")
			print("I'm not equipped to handle this.")
			exit()
		retmsg  = notable_adapters[0]

	return retmsg


def run_cutadapt(cutadapt, adapter_seq, adapter_name, dna_dir, root, labels, out_dir):
	if adapter_name not in adapter_seq or not adapter_seq[adapter_name]:
		print(f"This piece of code is a stump: need to find systematic way to find adapter seqs.")
		exit()
	aseq = adapter_seq[adapter_name]
	fastq   = [f"{dna_dir}/{root}_{label}.fastq" for label in labels]
	trimmed = [f"{out_dir}/{root}_trimmed_{label}.fastq" for label in labels]

	if not os.path.exists(trimmed[0]) or not os.path.exists(trimmed[1]):
		print(f"running cutadapt for {fastq}")
		cmd  = f"{cutadapt} --quiet -a {aseq} -A {aseq} -o {trimmed[0]} -p {trimmed[1]}  {fastq[0]} {fastq[1]}"
		subprocess.call(["bash", "-c", cmd])
		# TODO check nonempty output present after the run
	else:
		print(f"found {trimmed}")

	return f"{root}_trimmed"


def adapter_cleanup(cutadapt, adapter_seq, rootnames, read_labels, dna_dir, fastqc_out_dir, out_dir):
	if not os.path.exists(out_dir): os.mkdir(out_dir)
	flagged = {}
	for root in rootnames:
		for label in read_labels:
			fastqc_data_file = f"{fastqc_dir(fastqc_out_dir, root, label)}/fastqc_data.txt"
			check_exist([fastqc_data_file])
			retmsg = parse_fastqc_output(fastqc_data_file)
			if retmsg == "ok":
				print(f"\t {root} {label} ok")
				continue
			else:
				adapter = retmsg
				print(f"\t {root} {label} problematic adapter: {adapter}")
				if not root in flagged: flagged[root] = set()
				flagged[root].add(adapter)

	trimmed_root_names = []
	for root in flagged:
		if len(flagged[root])>1:
			print(f"{root} flagged for the presence of multiple addapters, ({flagged[root]}).")
			print("I'm not equipped to deal with it - bailing out.")
			exit()
		new_root_name = run_cutadapt(cutadapt, adapter_seq, flagged[root].pop(), dna_dir, root, read_labels, out_dir)
		trimmed_root_names.append(new_root_name)

	return trimmed_root_names


def adapter_sanity_check(trimmed_dir,  trimmed_rootnames, read_labels, fastqc_out_dir):
	if not trimmed_rootnames: return
	if not os.path.exists(fastqc_out_dir): os.mkdir(fastqc_out_dir)
	fastqc = Config.fastqc
	unzip = Config.unzip
	fastqc_deps = [fastqc, unzip, trimmed_dir, fastqc_out_dir]
	run_fastqc(fastqc_deps, trimmed_rootnames, read_labels)
	# check the summary file for each fastq
	# print issues to stdout
	talk_fastqc(trimmed_rootnames, read_labels, fastqc_out_dir)


def main():

	home_path = get_home_path()

	cutadapt     = Config.cutadapt
	adapter_seq  = Config.adapter_seq
	dna_dir          = f"{home_path}/task_dna"
	fastqc_first_dir = f"{home_path}/fastqc/first_pass"  # should exitst previously
	clean_fastq_dir  = f"{home_path}/clean_fastq"
	fastqc_this_dir  = f"{home_path}/fastqc/trimmed"  # the output in this round

	rootnames   = Config.paired_reads_rootnames
	read_labels = Config.read_labels

	check_exist([cutadapt, dna_dir, fastqc_first_dir])

	check_fastq_exist(dna_dir)

	trimmed_root_names = adapter_cleanup(cutadapt, adapter_seq, rootnames, read_labels,
											dna_dir, fastqc_first_dir, clean_fastq_dir)
	# what does fastqc have to say now?
	adapter_sanity_check(clean_fastq_dir, trimmed_root_names, read_labels, fastqc_this_dir)

	return


if __name__ == "__main__":
	main()
