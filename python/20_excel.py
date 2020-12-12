#! /usr/bin/python3

""" Compile the pipeline output into a spreadsheet table

Part of seqinspector toy NGS QC pipeline
Ivana Mihalek,  2020

This this script assumes the directory tree of the format
.
├── alignments
├── annotation
├── clean_fastq
├── fastqc
├── pileup
├── python
└── task_dna

where task_dna contains the fastq files (the name chosen to correspond
with the original task folder). Other dependencies shoud be set in config.py.
"""


import matplotlib.pyplot as plt
import xlsxwriter

from pileup import *
from utils import *
from styling import Styling


class Annotation:
	def __init__(self, region, gene):
		self.gene   = gene
		self.chrom  = None
		self.__region = region

	def region(self):
		if self.chrom and self.chrom in self.__region:
			return self.__region[self.chrom]
		else:
			return []


def read_region_annotation(ncbi_refseq):
	annotation = {}
	with open(ncbi_refseq) as inf:
		for line in inf:
			fields  = line.strip().split("\t")
			if len(fields) != 4: continue  # not interested
			[chrom, start, end, annotstr] = fields
			chrom_number = int(chrom.replace("chr",""))
			if chrom_number not in annotation: annotation[chrom_number] = {}
			annotation[chrom_number][f"{start} {end}"] = annotstr
	return annotation


def read_gene_annotation(gene_annot_file):
	column_names = ["gene_symbol", "name", "uniprot", "expression", "disease"]
	annotation = {}
	with open(gene_annot_file) as inf:
		for line in inf:
			columns  = line.strip().split("\t")
			annotation[columns[0]] = dict(zip(column_names[1:], columns[1:]))

	return annotation


def read_annotation(ncbi_refseq, gene_annot_file):
	return Annotation(read_region_annotation(ncbi_refseq), read_gene_annotation(gene_annot_file))


def find_annotation(annotation, start, end):
	for region, annotstr in annotation.region().items():
		[annot_start, annot_end] = [int(i) for i in region.split()]
		if annot_start<=start<=annot_end or annot_start<=end<=annot_end:
			return annotstr
	return "annot not found"


def get_gene_name(annotation, start, end):
	# paranoia to worry about two gene name fields in here?
	annotation_string = find_annotation(annotation, start, end)
	gene_name = list(filter(lambda field: "gene_name" in field, annotation_string.split(";")))[0]
	return gene_name.strip().split()[1].replace("\"", "")


def group_by_gene(stats, region_end, annotation):
	# there are several statistics on this list, they shoul
	# all have region start as their first argument
	gene = {}
	for start in sorted(stats[0].keys()):
		end = region_end[start]
		gene_name = get_gene_name(annotation, start, end)
		if gene_name not in gene: gene[gene_name] = []
		gene[gene_name].append(f"[hg19] {start}-{end}")
	return gene


def write_datasets(worksheet, style, stats, region, row, column):
	[coverage, avg_depth, calls_per_interval] = stats
	rows_spanned = 0
	# avg_depth is organized by region start, given as integer;
	# regions ids a string of the format f"{assembly} {start}-{end}"
	start = int(re.split(" |-", region)[1])
	row_span = len(avg_depth[start])
	if len(coverage[start]) != row_span:
		print(f"mismatch in the number of datasets at pos: {start}")
		exit()
	# I think this is the most legible format for the table
	for i in [1, 0]:
		worksheet.write_string(row, column, "%.0f" % avg_depth[start][i])
		column += 1
	for i in [1, 0]:
		worksheet.write_string(row, column, "%.2f" % coverage[start][i])
		column += 1

	worksheet.write_string(row, column, "%d" % calls_per_interval[start][1])
	column += 1

	# draw attention to calls in AH that fall in regions not covered by CH
	if calls_per_interval[start][0]>0 and avg_depth[start][0]>10 and avg_depth[start][0]>avg_depth[start][1]:
		worksheet.write_string(row, column, "%d" % calls_per_interval[start][0], style.xlsx_format["red_border"])
	else:
		worksheet.write_string(row, column, "%d" % calls_per_interval[start][0])
	rows_spanned += 1
	return rows_spanned


def write_region(worksheet, style, stats, regions, annotation, row, column):
	# TODO which exon(s) does the region correspond to
	rows_spanned = 0
	for region in regions:
		row_span = write_datasets(worksheet, style, stats, region, row, column+1)
		if not row_span: continue
		if row_span == 1:
			worksheet.write_string(row, column, region)
		else:
			worksheet.merge_range(row, column, row+row_span-1, column, region)
		row += row_span
		rows_spanned += row_span
	return rows_spanned


def write_chromosome(worksheet, style, stats, region_end, annotation, row, column):
	gene = group_by_gene(stats, region_end, annotation)
	rows_spanned = 0
	for gene_name, regions in gene.items():
		row_span = write_region(worksheet, style, stats, regions, annotation, row, column+3)
		if not row_span: continue
		disease = annotation.gene[gene_name]["disease"]
		uniprot = annotation.gene[gene_name]['uniprot']
		uniprot_hyperlink = f"https://www.uniprot.org/uniprot/{uniprot}"
		if row_span == 1:
			worksheet.write_string(row, column, gene_name)
			worksheet.write_url(row, column+1, uniprot_hyperlink, string=uniprot)
			worksheet.write_string(row, column+2, disease)
		else:  # merge range cannot have range of 1 (duh ...)
			worksheet.merge_range(row, column, row+row_span-1, column, gene_name)
			worksheet.merge_range(row, column+1, row+row_span-1, column+1, "", style.xlsx_format["hyperlink"])
			worksheet.write_url(row, column+1, uniprot_hyperlink, string=uniprot)
			worksheet.merge_range(row, column+2, row+row_span-1, column+2, disease)
		row += row_span
		rows_spanned += row_span
	return rows_spanned


def read_calls_per_interval(calls_file):
	calls_per_interval = {}
	with open(calls_file) as inf:
		for line in inf:
			[chrom, start, count1, count2] = [int_cast(x, calls_file) for x in line.strip().split("\t")]
			if chrom not in calls_per_interval: calls_per_interval[chrom] = {}
			calls_per_interval[chrom][start] = [count1, count2]
	return calls_per_interval


################
def column_string(idx):
	# setting width does not like column indices
	char = chr(ord('A')+idx)
	return f"{char}:{char}"


def set_column_widths(worksheet, style):
	# gene name, disease and region   columns
	worksheet.set_column(column_string(1), 12)
	worksheet.set_column(column_string(3), 50, style.xlsx_format["wordwrap"])
	worksheet.set_column(column_string(4), 30)


def write_header(worksheet, style):
	header = ["chrom", "gene", "uniprot", "disease", "region",
				"depth in CH", "depth in AH", "coverage in CH", "coverage in AH", "calls in CH", "calls in AH", ]
	worksheet.set_row(0, 40, style.xlsx_format["header"])
	for column in range(len(header)):
		worksheet.write_string(0, column, header[column])


def main():

	home_path = get_home_path()

	cvg_dir = f"{home_path}/pileup/coverage"
	vcf_dir = f"{home_path}/pileup/variants"
	annotation_dir       = f"{home_path}/annotation"
	calls_file           = f"{vcf_dir}/calls_per_interval.tsv"
	ncbi_refseq          = f"{annotation_dir}/target_regions_annotated.bed"
	gene_annotation_file = f"{annotation_dir}/gene_annotation.tsv"
	check_exist([cvg_dir, annotation_dir, calls_file, ncbi_refseq, gene_annotation_file])

	annotation = read_annotation(ncbi_refseq, gene_annotation_file)

	# summarize coverage/depth for each region
	[coverage, avg_depth, region_end] = summarize_coverage_pileup(cvg_dir)
	calls_per_interval = read_calls_per_interval(calls_file)

	# Create an new Excel file and add a worksheet.
	workbook  = xlsxwriter.Workbook('region_summary.xlsx')
	worksheet = workbook.add_worksheet("Summary for  AH and CH")
	style = Styling(workbook)
	set_column_widths(worksheet, style)
	write_header(worksheet, style)

	# write
	row    = 1
	column = 0
	for chrom in sorted(coverage.keys()):
		stats = [coverage[chrom], avg_depth[chrom], calls_per_interval[chrom]]
		annotation.chrom = chrom
		row_span = write_chromosome(worksheet, style, stats, region_end, annotation, row, column+1)
		if not row_span: continue
		# merge_range(first_row, first_col, last_row, last_col, data, cell_format=None)
		if row_span == 1:
			worksheet.write_string(row, column, f"chr{chrom}")
		else:
			worksheet.merge_range(row, column, row+row_span-1, column,  f"chr{chrom}")

		row += row_span


	workbook.close()


	return


if __name__ == "__main__":
	main()
