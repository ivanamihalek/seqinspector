

""" Utilities for samtools/bcftools pileup manipulation

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

"""

import subprocess, re, os


def output_pileup_vcf(bcftools, ref_genome, pileup_dir, bamfiles):
	for bamfile in bamfiles:
		root = bamfile.split("/")[-1][:-4]
		# https://samtools.github.io/bcftools/howtos/variant-calling.html
		# z: compressed vcf
		# bcftools insists these should be indexed - not sure why it means they should be compressed
		# https://github.com/samtools/bcftools/issues/668
		cmd  = f"{bcftools} mpileup -f {ref_genome} {bamfile}  --max-depth 10000 | "
		cmd += f"{bcftools} call -mv -Oz -o {pileup_dir}/{root}.vcf.gz  2> /dev/null"
		print(cmd)
		subprocess.call(["bash", "-c", cmd])
		# not sure why bcftools does not do it itself, bcs in bcftools view it is asking for index
		cmd  = f"{bcftools} index {pileup_dir}/{root}.vcf.gz  2> /dev/null"
		subprocess.call(["bash", "-c", cmd])


def output_pileup_cvg(cvg_dir, samtools, merged_regions, bamfiles):
	for chrom, intervals in merged_regions.items():
		for intv in intervals:
			fromto = f"{intv[0]}-{intv[1]}"
			print(f"chr{chrom}:{fromto}")
			# Each input file produces a separate group of pileup columns in the output.
			# -a : Output all positions, including those with zero depth.
			cmd  = f"{samtools} mpileup -a -r  chr{chrom}:{fromto} "
			cmd += f"-o {cvg_dir}/pileup_chr{chrom}_{fromto}.tsv {bamfiles[0]} {bamfiles[1]}   2> /dev/null"
			# print(cmd)
			subprocess.call(["bash", "-c", cmd])


def summarize_variant_pileup(bcftools, big_vcf_file, regions):
	# genotype deciphering scheme https://samtools.github.io/hts-specs/VCFv4.1.pdf, p 6
	#  / is unresolverd and | resolved phenotype
	# alleles are indexed from 0; the list order of genotype j/k is j+k*(k+1)/2) (upper triangular?)
	# eg for r biallelic sites the ordering is:  AA,AB,BB; for triallelic:  AA,AB,BB,AC,BC,CC, etc.
	phred_likelihoods = {}
	for chrom, intervals in regions.items():
		print(chrom, len(intervals))
		for intv in intervals:
			fromto = f"{intv[0]}-{intv[1]}"
			# -H: suppress header
			# -O: output format v (plain VCF)
			# QUAL gives an estimate of how likely it is to observe a call purely by chance
			# '%QUAL>=20' is some zeroth order check for the meaninfulness of the call
			# see https://samtools.github.io/bcftools/howtos/variant-calling.html
			cmd = f"{bcftools} view -H  -Ov -i '%QUAL>=20' -r chr{chrom}:{fromto} {big_vcf_file} "
			output = subprocess.check_output(["bash", "-c", cmd])
			if not output: continue
			for line in output.decode().split("\n"):
				if not line: continue
				fields = line.split("\t")
				position = fields[1]
				phred_likelihood = fields[-1]
				if chrom not in phred_likelihoods: phred_likelihoods[chrom] = {}
				if intv[0] not in phred_likelihoods[chrom]: phred_likelihoods[chrom][intv[0]] = []
				phred_likelihoods[chrom][intv[0]].append([position, phred_likelihood])
	return phred_likelihoods


def summarize_coverage_pileup(cvg_dir):
	pileups = [fnm for  fnm in os.listdir(cvg_dir) if "pileup" in fnm]
	coverage   = {}
	avg_depth  = {}
	region_end = {}
	max_depth = -1
	for fnm in pileups:
		fields = re.split("_|-", fnm.replace("pileup_", "").replace(".tsv", "").replace("chr", ""))
		[chrom, start, end] = [int(f) for f  in  fields]
		region_end[start] = end
		depth_total   = [0, 0]  # the first column refers to bamfiles[0], the second to  bamfiles[1]
		depth_nonzero = [0, 0]
		region_length = 0
		with open(f"{cvg_dir}/{fnm}") as inf:
			for line in inf:
				field = line.split("\t")
				# note this works only if there are no special flags in output_pileup mpileup command
				# (if there are, the column content is different)
				# see http://www.htslib.org/doc/samtools-mpileup.html
				depth = [int(field[3]), int(field[6])]
				region_length += 1
				for i in range(2):
					depth_total[i]   += depth[i]
					depth_nonzero[i] += 1 if depth[i]>10 else 0
					if max_depth<depth[i]: max_depth=depth[i]

		if region_length:
			if chrom not in coverage:
				coverage[chrom]  = {}
				avg_depth[chrom] = {}
			coverage[chrom][start]  = [0, 0]
			avg_depth[chrom][start] = [0, 0]
			for i in range(2):
				coverage[chrom][start][i]  = depth_nonzero[i]/region_length
				avg_depth[chrom][start][i] = depth_total[i]/region_length

	# print(max_depth)
	# fmt = "%3d  %10d      %.2f  %8.3f       %.2f  %8.3f  "
	# for chrom in sorted(coverage.keys()):
	# 	for start in sorted(coverage[chrom].keys()):
	# 		print(fmt%(chrom, start, coverage[chrom][start][0], avg_depth[chrom][start][0],
	# 				 coverage[chrom][start][1], avg_depth[chrom][start][1]))
	return [coverage, avg_depth, region_end]
