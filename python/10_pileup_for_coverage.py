#! /usr/bin/python3

""" Run samtools pileup as input for depth and coverage calc

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

This script will produce pielup output, while the following script,
11_depth_and_coverage.py, will actually crunch the numbers. The target
regions for the two sets are merged if overlaping, in the attempt
to provide compact comparison.

Sources:
http://www.htslib.org/doc/samtools-mpileup.html

This this script assumes the directory tree of the format
.
├── alignments
├── clean_fastq
├── fastqc
├── python
└── task_dna

where task_dna contains the fastq files (the name chosen to correspond
with the original task folder). This script will add directory called
"pileup/coverage."  Other dependencies shoud be set in config.py.

"""


import matplotlib.pyplot as plt
from pileup import *

# fastq-specific functions - imports utils
from fastqc import *


def add_interval(regions, new_interval):
	place_found = False
	for i in range(len(regions)):
		interval = regions[i]
		if new_interval[1]<interval[0]:
			# insert at i
			regions.insert(i, new_interval)
			place_found = True
			break
		elif new_interval[0]<=interval[1]:
			# merge
			regions[i] = [min(new_interval[0], interval[0]), max(new_interval[1], interval[1])]
			place_found = True
			break
	# append
	if not place_found: regions.append(new_interval)


def find_place(regions, qry, sanity_check=False):
	target_interval = None
	for interval in regions:
		if interval[0]<=qry[0] and qry[1]<=interval[1]:
			if sanity_check:
				if interval[0]< qry[0] or qry[1] < interval[1]:
					print(f"{qry}  belongs to {interval}")
				else:
					print(f"{qry} is equal to {interval}")
			target_interval = interval
			break
	if not target_interval:
		print(f"bug in interval manipulation: place not found for {qry}")
		exit()
	return target_interval


def int_cast(x, fnm):
	try:
		return int(x)
	except:
		print(f"unexpected non-integer in {fnm}")


def regions_merge(region_fnms, sanity_check=False):
	merged_regions = {}
	regions = {}
	for regionfnm in region_fnms:
		regions[regionfnm] = {}
		inf = open(regionfnm)
		for line in inf:
			if "seqname" in line.lower(): continue  # header
			fields = line.strip().split()[:3]
			if len(fields) < 3: continue
			[chrom, start, end]  = [int_cast(field, regionfnm) for field in fields]
			if chrom not in merged_regions: merged_regions[chrom] = []
			add_interval(merged_regions[chrom], [start, end])
			if chrom not in regions[regionfnm]: regions[regionfnm][chrom] = []
			regions[regionfnm][chrom].append([start, end])
		inf.close()

	if sanity_check:
		print(f"number of merged regions: {sum([len(v) for v in merged_regions.values()])}")
		for regionfnm in region_fnms:
			# print(regionfnm)
			for chrom, intervals in regions[regionfnm].items():
				# print("\t", chrom)
				for intv in intervals:
					# will exit on failure
					find_place(merged_regions[chrom], intv, sanity_check=False)

	return merged_regions


def output_merged_regions(cvg_dir, merged_regions):
	with open(f"{cvg_dir}/merged_target_regions.bed", "w") as outf:
		for chrom in sorted(merged_regions.keys()):
			for intv in merged_regions[chrom]:
				print(f"chr{chrom}\t{intv[0]}\t{intv[1]}", file=outf)


def main():

	home_path = get_home_path()

	samtools    = Config.samtools
	alnmts_dir  = f"{home_path}/alignments"
	dna_dir     = f"{home_path}/task_dna"
	rootnames   = Config.paired_reads_rootnames
	bamfiles    = [f"{alnmts_dir}/{rootnm}.dedup.bam" for rootnm in rootnames]
	region_fnms = [f"{dna_dir}/{rootnm.replace('L001','target')}.txt" for rootnm in rootnames]
	check_exist([samtools, alnmts_dir, dna_dir] + bamfiles + region_fnms)

	merged_regions = regions_merge(region_fnms, sanity_check=False)

	pileup_dir = f"{home_path}/pileup/coverage"
	os.makedirs(pileup_dir, exist_ok=True)

	# output the merged regions for input in IGV
	output_merged_regions(pileup_dir, merged_regions)

	# run samtools pileup on each of the regions
	output_pileup_cvg(pileup_dir, samtools, merged_regions, bamfiles)

	return


if __name__ == "__main__":
	main()
