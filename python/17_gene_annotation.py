#! /usr/bin/python3

""" Gene-disease correlation table

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020


This script assumes the availability of compact uniprot
table with gene/disease correlation. Contact the author
to obtain it. Or try here:

https://www.uniprot.org/uniprot/?query=reviewed:yes#customize-columns

This this script assumes the directory tree of the format
.
├── alignments
├── annotation
├── clean_fastq
├── fastqc
├── pileup
├── python
└── task_dna

In particular, note the annotation directory, which should have
been populated in the previous step, 16_annotation2bed.bash.
Other dependencies shoud be set in config.py.
"""


import matplotlib.pyplot as plt
from pileup import *
from utils import *


def get_genes(ncbi_refseq):
	gene_names = set()
	with open(ncbi_refseq) as inf:
		for line in inf:
			fields  = line.strip().split()
			gene_name = fields[-1].strip("\";")
			gene_names.add(gene_name)
	return list(gene_names)


def shorten_disease_string(disease_str, split_on):
	# if there are more than 100 characters, take only the first 100
	# if there are more than three sentences, take only first three
	return ".".join(disease_str.split(split_on)[0].replace("Note=","")[:100] .strip().split(".")[:3])


def disease_parse(uniprot_disease_string):
	if not uniprot_disease_string: return "unknown"
	if "[MIM" in uniprot_disease_string:
		split_on = "[MIM"
	else:
		split_on = "."
	diseases  = [shorten_disease_string(substr, split_on) for substr in uniprot_disease_string.split("DISEASE:")[1:]]
	return "; ".join(diseases)


def find_description(uniprot_info, gene_symbols):
	# I happen to have this one on disk :)
	# fields
	# uniprot id, protein name, gene names (the first should be hgnc;
	# but then several genes can correspond to the same roteins, e.g. histones),
	# length, expression tissue, involvement in disease
	gene_symbols = set(gene_symbols)
	name = {}
	expression = {}
	uniprot = {}
	disease = {}
	with open(uniprot_info) as inf:
		for line in inf:
			field = line.strip().split("\t")
			if len(field)<6: continue
			if len(field[2].replace(" ", "")) == 0: continue
			uniprot_gene_names  = set([x.strip(";,") for x in field[2].split()])
			for hgnc in uniprot_gene_names.intersection(gene_symbols):
				uniprot[hgnc]    = field[0]
				name[hgnc]       = field[1].split("(")[0].strip()
				expression[hgnc] = field[4].replace("TISSUE SPECIFICITY: ","")
				disease[hgnc]    = disease_parse(field[5])
	return [name, expression, disease, uniprot]


def main():

	home_path = get_home_path()

	cvg_dir = f"{home_path}/pileup/coverage"
	annotation_dir = f"{home_path}/annotation"
	ncbi_refseq    = f"{annotation_dir}/target_regions_annotated.bed"
	uniprot_info   = Config.basic_uniprot
	check_exist([cvg_dir, annotation_dir, ncbi_refseq, uniprot_info])

	gene_symbols = get_genes(ncbi_refseq)

	[name, expression, disease, uniprot] = find_description(uniprot_info, gene_symbols)
	with open(f"{annotation_dir}/gene_annotation.tsv", "w") as outf:
		for gene_symbol in sorted(gene_symbols):
			if gene_symbol in name:
				fields = [gene_symbol, name[gene_symbol], uniprot[gene_symbol], expression[gene_symbol], disease[gene_symbol]]
			else:
				fields = [gene_symbol, "unknown", "unknown", "unknown", "unknown"]
			print("\t".join(fields), file=outf)


if __name__ == "__main__":
	main()
