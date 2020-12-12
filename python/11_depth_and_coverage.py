#! /usr/bin/python3

""" Calculate depth and coverage

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

Calculate depth and coverage for the pileup files generated in the
previous script, 10_pileup_for_coverage.py, and create figure
for the report.

Sources:
http://www.htslib.org/doc/samtools-mpileup.html
https://stackoverflow.com/questions/50571287/how-to-create-upside-down-bar-graphs-with-shared-x-axis-with-matplotlib-seabor

This this script assumes the directory tree of the format
.
├── alignments
├── clean_fastq
├── fastqc
├── pileup
├── python
└── task_dna

where task_dna contains the fastq files (the name chosen to correspond
with the original task folder). Other dependencies shoud be set in config.py.

"""

import matplotlib.pyplot as plt
from pileup import *

# fastq-specific functions - imports utils
from fastqc import *


def make_graph(coverage, avg_depth):
	number_of_regions = sum([len(v) for v in coverage.values()])
	height = [[], []]  # we will set the bar height to be proportional to the sequencing depth in each region
	alpha  = [[], []]  # and transparency to the coverage of the region
	for chrom in sorted(coverage.keys()):
		for start in sorted(coverage[chrom].keys()):
			for i in range(2):
				dpth = avg_depth[chrom][start][i]
				bar_height = dpth  # max(dpth, 0.1), no need, matplotlib know how to handle 0s on the log scale
				height[i].append(bar_height*(-1)**(1-i))  # invert the first set of depths (correspoding to 'AH' set)
				# alpha[i].append(1-coverage[chrom][start][i])
				alpha[i].append(0.5 if coverage[chrom][start][i]<0.3 else 1)

	# sort by height in the first set which I happen to know has more points
	sorted_height = [[], []]
	idx = range(number_of_regions)
	sorted_idx = sorted(idx, key=lambda i: height[1][i], reverse=True)
	sorted_height[0] = [height[0][i] for i in sorted_idx]
	sorted_height[1] = [height[1][i] for i in sorted_idx]

	fig, ax = plt.subplots()

	color = [(0.0, 0.0, 1.0, a) for a in alpha[0]]  # blue with alpha inverse prop to coverage
	ax.bar(idx, sorted_height[0], color=color, label="AH")

	color = [(1.0, 0.0, 0.0, a) for a in alpha[0]]  # red with alpha inverse prop to coverage
	ax.bar(idx, sorted_height[1], color=color, label="CH")

	# tickamrs: logarithmic and positive
	plt.yscale("symlog")  # symmetric log scale (log in negative direction) - do before tkaing abs value of ticks
	ticks = ax.get_yticks()
	ax.set_yticklabels([int(abs(tick)) for tick in ticks])  # the lower part of the graph is not really negative

	# axis labels
	plt.xlabel('regions, sorted by depth in the CH set', fontsize=16)
	plt.ylabel('depth', fontsize=16)

	# legend
	ax.legend()

	fig.tight_layout()  # adjusts the padding between and around subplots.


	plt.show()


def main():

	home_path = get_home_path()

	samtools    = Config.samtools
	alnmts_dir  = f"{home_path}/alignments"
	cvg_dir     = f"{home_path}/pileup/coverage"
	dna_dir     = f"{home_path}/task_dna"
	rootnames   = Config.paired_reads_rootnames
	bamfiles    = [f"{alnmts_dir}/{rootnm}.dedup.bam" for rootnm in rootnames]
	region_fnms = [f"{dna_dir}/{rootnm.replace('L001','target')}.txt" for rootnm in rootnames]
	check_exist([samtools, alnmts_dir, dna_dir, cvg_dir] + bamfiles + region_fnms)

	# summarize coverage/depth for each region
	[coverage, avg_depth] = summarize_coverage_pileup(cvg_dir)[:2]

	# visual summary of depth and coverage for the two samples
	make_graph(coverage, avg_depth)

	return


if __name__ == "__main__":
	main()
